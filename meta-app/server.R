library(plotly)
library(SummarizedExperiment)
library(data.table)
library(PCAtools)
library(BiocSingular)
library(dplyr)
library(pool)
library(DT)


# global data
se <- readRDS("data/se.rds")
pathways <- readRDS("data/pathways.rds")
dge_se <- se[rowData(se)$feature_type == "gene", ]
dre_se <- se[rowData(se)$feature_type == "re", ]

# Connect to DE database
pool <- dbPool(
  drv = RSQLite::SQLite(),
  dbname = "data/db.sqlite"
)

# ------------------------------------------------------------------------------
server <- function(input, output, session) {
  data <- eventReactive(input$run, {
    showModal(modalDialog("Processing data. Please wait...", footer = NULL))
    
    # Specify SummarizedExperiment dataset
    if (input$dataset == "Protein Coding Genes") {
      se <- dge_se
    } else {
      se <- dre_se
    }

    # Filter SE based on user input
    se_filtered <- se[, se$contrast %chin% input$contrasts &
      se$experiment %chin% input$experiments &
      se$cell_line %chin% input$cell_lines &
      se$drug %chin% input$drugs &
      se$tissue %chin% input$tissues &
      se$disease %chin% input$diseases &
      se$epigenetic_class %chin% input$epigenetic]
    
    # If a gene set is selected subset the SE object
    if (input$pathway != "ALL GENES") {
      if (input$dataset == "Protein Coding Genes") {
        se_filtered <- se_filtered[rownames(se_filtered) %chin% pathways[[input$pathway]], ]
      }
    }

    # Extract matrix of user specified assay
    M <- assay(se_filtered, input$data)

    # Filter matrix if complete cases is specified
    if (input$complete) {
      M <- na.omit(M)
    }

    # Perform PCA
    algo <- switch (input$algorithm,
      "Exact" = ExactParam(), 
      "Irlba" = IrlbaParam(), 
      "Random" = RandomParam(), 
      "Auto" = FastAutoParam()
    )
    
    # Replace NA values with imputed values
    replace_val <- switch (input$data,
      "lfc" = 0L,
      "fdr" = 1L
    )
    M[is.na(M)] <- replace_val
    
    # Run PCA
    pca_res <- pca(
      mat = M,
      metadata = as.data.frame(colData(se_filtered)),
      scale = input$scale,
      center = input$center,
      removeVar = input$var,
      rank = input$components,
      BSPARAM = algo
    )

    removeModal()
    list("PCA" = pca_res, "ids" = colnames(M), "features" = rownames(M), "se" = se_filtered)
  })

  # Create single reactive data.table for use in plotting functions
  dt <- reactive({
    pca_dt <- setDT(data()[["PCA"]]$rotated, keep.rownames = "sample_name")
    meta_dt <- setDT(data()[["PCA"]]$metadata, keep.rownames = "sample_name")
    dt <- merge.data.table(x = pca_dt, y = meta_dt, by = "sample_name", all.x = TRUE)
    dt
  })

  output$biplot <- renderPlotly({
    dt <- dt()
    p <- plot_ly(
      dt,
      x = ~ get(input$x_axis),
      y = ~ get(input$y_axis),
      color = ~ get(input$col),
      colors = "Set1",
      hovertext = paste(
        "BioProject:", dt$experiment,
        "<br>Cell Line:", dt$cell_line,
        "<br>Disease:", dt$disease,
        "<br>Drug:", dt$drug,
        "<br>Dose", dt$dose,
        "<br>Time", dt$time_hr,
        "<br>Epigenetic Class", dt$epigenetic_class
      )
    ) |>
      add_markers(size = 180) |>
      layout(
        dragmode = "lasso",
        showlegend = FALSE,
        xaxis = list(title = paste(input$x_axis)),
        yaxis = list(title = paste(input$y_axis))
      )

    event_register(p, "plotly_selected")

    p
  })

  # Table
  output$table <- DT::renderDT({
    d <- event_data("plotly_selected")
    if (!is.null(d)) {
      xmin <- min(d$x) - 1e-6
      xmax <- max(d$x) + 1e-6
      ymin <- min(d$y) - 1e-6
      ymax <- max(d$y) + 1e-6
      dt <- dt()
      dt[
        data.table::between(get(input$x_axis), xmin, xmax, incbounds = TRUE) &
        data.table::between(get(input$y_axis), ymin, ymax, incbounds = TRUE),
        .(experiment, tissue, cell_line, disease, drug, dose, time_hr, epigenetic_class)
      ]
    }
  })

  # Meta Vote-Counting
  metavote_data <- reactive({
    keep_ids <- data()[["ids"]]
    keep_features <- data()[["features"]]

    # Filter the original SE object for only the selected features
    keep_se <- se[keep_features, keep_ids]
    lfc <- assay(keep_se, "lfc")
    fdr <- assay(keep_se, "fdr")

    # Calculate significance and direction based on both input matrices
    fdr_mat <- fdr < input$meta_pval
    lfc_mat <- matrix(
      data.table::fcase(
        lfc > 0 & abs(lfc) > input$meta_lfc, 1L,
        lfc < 0 & abs(lfc) > input$meta_lfc, -1L,
        default = 0L),
      nrow = nrow(keep_se)
      )
    
    exp_mat <- fdr_mat * lfc_mat

    # Calculate voting stats
    n_studies <- ncol(exp_mat)
    sign_consistency <- rowSums(exp_mat, na.rm = TRUE)
    n_de <- rowSums(abs(exp_mat), na.rm = TRUE)
    prop_de <- n_de / n_studies
    vote <- data.table::fcase(
      (abs(sign_consistency) >= input$meta_prop * n_studies) & (sign_consistency < 0), "down",
      (abs(sign_consistency) >= input$meta_prop * n_studies) & (sign_consistency > 0), "up",
      default = "unperturbed"
    )

    result_dt <- data.table::as.data.table(
      data.frame(
        "n_studies" = n_studies,
        "n_de" = n_de, 
        "prop_de" = round(prop_de, 2), 
        "sign_consistency" = sign_consistency, 
        "vote" = vote),
      keep.rownames = "feature_id"
    )
    result_dt
  })

  # Meta-Vote MetaVolcano plot
  output$metavolcano <- renderPlot({ coriell::plot_metavolcano(metavote_data()) })

  # Meta-Vote Count Table
  output$metatable <- DT::renderDT({
    metavote_data()[order(sign_consistency, decreasing = TRUE)]
  })
  
  # Differential Expression data
  d <- reactive({
    if (input$type == "both") {
      data <- pool |> 
        tbl("results") |> 
        filter(id == !!input$id) |> 
        collect()
    } else {
      data <- pool |> 
        tbl("results") |> 
        filter(id == !!input$id & feature_type == !!input$type) |> 
        collect()
    }
    return(data)
  })
  
  # Differential Expression volcano plot
  output$volcano <- renderPlot({ coriell::plot_volcano(d(), fdr = input$fdr, down_color = "blue", annotate_counts = FALSE) })
  
  # Differential Expression ma plot
  output$ma <- renderPlot({ coriell::plot_md(d(), x = "expr", fdr = input$fdr, annotate_counts = FALSE) })
  
  # Differential Expression Summary table
  output$de_tbl <- DT::renderDT({
    df <- coriell::summarize_dge(d(), fdr = input$fdr)  
    datatable(df) |> formatRound('Percent', 2)
    })
  
  # Experiment Ranking
  output$ranking <- DT::renderDT({
    se <- data()[["se"]]
    fdr <- assay(se, "fdr")
    lfc <- assay(se, "lfc")
    
    up <- (lfc > 0) & (fdr <= input$ranking_fdr)
    down <- (lfc < 0) & (fdr <= input$ranking_fdr)
    non <- fdr > 0.05
    
    up_counts <- colSums(up, na.rm = TRUE)
    down_counts <- colSums(down, na.rm = TRUE)
    non_counts <- colSums(non, na.rm = TRUE)
    
    df <- data.frame(up = up_counts, down = down_counts, non = non_counts)
    df <- cbind(df, colData(se)[, c("experiment", "cell_line", "drug", "dose", "time_hr", "epigenetic_class", "tissue", "disease")])
    df <- df[order(df$up, decreasing = TRUE), ]
    df
  })
  
  # Metadata table
  output$metadata <- DT::renderDT({
    df <- pool |> tbl("metadata") |> collect() |> dplyr::select(-id)
    datatable(
      df,
      colnames = c("Cell Line", "Drug", "BioProject", "Contrast", "Dose",
                   "Time (hr)", "Batch", "Mutation", "Comment", "Description",
                   "Epigenetic Class", "Tissue", "Disease"),
      options = list(columnDefs = list(list(
        targets = 4,
        render = JS(
          "function(data, type, row, meta) {",
          "return type === 'display' && data.length > 12 ?",
          "'<span title=\"' + data + '\">' + data.substr(0, 12) + '...</span>' : data;",
          "}")),
        list(
          targets = 10,
          width = "50%"
        ))),
      callback = JS('table.page(3).draw(false);')
      )
  })
}
