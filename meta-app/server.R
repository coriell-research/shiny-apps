suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(PCAtools))
suppressPackageStartupMessages(library(BiocSingular))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(pool))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(patchwork))


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
    algo <- switch(input$algorithm,
      "Exact" = ExactParam(),
      "Irlba" = IrlbaParam(),
      "Random" = RandomParam(),
      "Auto" = FastAutoParam()
    )

    # Replace NA values with imputed values
    replace_val <- switch(input$data,
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
        default = 0L
      ),
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
        "vote" = vote
      ),
      keep.rownames = "feature_id"
    )
    result_dt
  })

  # Meta-Vote MetaVolcano plot
  output$metavolcano <- renderPlot({
    coriell::plot_metavolcano(metavote_data())
  })

  # Meta-Vote Count Table
  output$metatable <- DT::renderDT({
    metavote_data()[order(sign_consistency, decreasing = TRUE)]
  })

  # Meta-Vote Results Download
  output$downloadMetaVote <- downloadHandler(
    filename = function() {
      paste("MetaVote-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      df <- metavote_data()[order(sign_consistency, decreasing = TRUE)]
      fwrite(df, file)
    }
  )

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
  output$volcano <- renderPlot({
    coriell::plot_volcano(d(), fdr = input$fdr, lfc = log2(1.5), annotate_counts = FALSE)
  })

  # Differential Expression ma plot
  output$ma <- renderPlot({
    coriell::plot_md(d(), x = "expr", fdr = input$fdr, lfc = log2(1.5), annotate_counts = FALSE)
  })

  # Differential Expression Summary table
  output$de_tbl <- DT::renderDT({
    df <- coriell::summarize_dge(d(), fdr = input$fdr)
    datatable(df) |> formatRound("Percent", 2)
  })

  output$downloadDE <- downloadHandler(
    filename = paste0("DifferentialExpression-", Sys.Date(), ".csv"),
    content = function(file) {
      fwrite(d(), file)
    }
  )
  
  output$downloadDEPlots <- downloadHandler(
    filename = paste0("DifferentialExpression-", Sys.Date(), ".png"),
    content = function(file) {
      png(file, width = 14, height = 6, units = "in", res = 300)
      ptitle <- gsub(".", " ", input$id, fixed = TRUE)
      if (input$type == "both") {
        ftype <- "All Features"
      } else if (input$type == "gene") {
        ftype <- "Protein Coding Genes"
      } else {
        ftype <- "Repetitive Elements"
      }
      v <- coriell::plot_volcano(d(), fdr = input$fdr)
      m <- coriell::plot_md(d(), x = "expr", fdr = input$fdr)
      print(
        (v | m) + plot_layout(guides = "collect") + 
        plot_annotation(title = ptitle,
                        subtitle = ftype) & 
        coriell::theme_coriell()
        )
      dev.off()
      },
    contentType = "image/png"
  )

  # Experiment Ranking
  ranking_dt <- reactive({
    se <- data()[["se"]]
    fdr <- assay(se, "fdr")
    lfc <- assay(se, "lfc")

    up <- (lfc > 0) & (fdr <= input$ranking_fdr)
    down <- (lfc < 0) & (fdr <= input$ranking_fdr)
    non <- fdr > input$ranking_fdr

    up_counts <- colSums(up, na.rm = TRUE)
    down_counts <- colSums(down, na.rm = TRUE)
    non_counts <- colSums(non, na.rm = TRUE)

    df <- data.frame(up = up_counts, down = down_counts, non = non_counts)
    df <- cbind(df, colData(se)[, c("experiment", "cell_line", "drug", "dose", "time_hr", "epigenetic_class", "tissue", "disease")])
    df <- df[order(df$up, decreasing = TRUE), ]
    df
  })

  output$ranking <- DT::renderDT({
    ranking_dt()
  })
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("ExpressionRanking-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      fwrite(ranking_dt(), file)
    }
  )

  # GSEA
  ranks <- reactive({
    df <- pool |>
      tbl("results") |>
      filter(id == !!input$gsea_id & feature_type == "gene") |>
      collect()

    ranks <- df$logFC * -1 * log10(df$PValue)
    names(ranks) <- df$feature_id
    ranks
  })

  output$gsea_results <- DT::renderDT({
    res <- fgsea(
      pathways = pathways[input$gsea_pathway],
      stats = ranks(),
      eps = 0.0,
      minSize = 10,
      maxSize = 500,
      nPermSimple = input$gsea_nPerm
    )
    res[, leadingEdge := NULL]
    datatable(res) |> formatRound(c("pval", "padj", "log2err", "ES", "NES"), digits = 2)
  })

  output$enrichment_plot <- renderPlot({
    plotEnrichment(
      pathways[[input$gsea_pathway]],
      ranks()
    ) +
      geom_line(color = "darkorange2", size = 2) +
      labs(title = input$gsea_pathway) +
      theme(plot.title = element_text(size = 14, face = "bold"))
  })

  # UpSet plot
  upset_data <- reactive({
    data <- if (input$upset_type == "both") {
      data <- pool |>
        tbl("results") |>
        filter(id %in% !!input$upset_ids) |>
        collect()
    } else {
      data <- pool |>
        tbl("results") |>
        filter(id %in% !!input$upset_ids & feature_type == !!input$upset_type) |>
        collect()
    }
    return(data)
  })

  output$upset <- renderPlot({
    dt <- setDT(upset_data())
    up <- dt[FDR < input$upset_fdr & logFC > 0]
    down <- dt[FDR < input$upset_fdr & logFC < 0]
    up_by_con <- split(up, by = "contrast")
    down_by_con <- split(down, by = "contrast")

    names(up_by_con) <- gsub("PRJNA[0-9]+\\.", "", names(up_by_con))
    names(up_by_con) <- gsub("PRJNA[0-9]+\\.", "", names(up_by_con))
    names(up_by_con) <- paste(names(up_by_con), "up", sep = ".")
    names(down_by_con) <- paste(names(down_by_con), "down", sep = ".")

    l_up <- lapply(up_by_con, function(x) unique(x$feature_id))
    l_down <- lapply(down_by_con, function(x) unique(x$feature_id))
    l <- c(l_up, l_down)
    m <- make_comb_mat(list_to_matrix(l), mode = input$upset_mode)
    m <- m[comb_degree(m) == 2]

    hm <- UpSet(m,
      comb_order = order(comb_size(m), decreasing = TRUE),
      top_annotation = upset_top_annotation(m, add_numbers = TRUE, numbers_gp = grid::gpar(fontsize = 16)),
      right_annotation = upset_right_annotation(m, add_numbers = TRUE, numbers_gp = grid::gpar(fontsize = 16))
    )
    draw(hm, padding = unit(c(12, 12, 12, 12), "mm"))
  })

  # Metadata table
  output$metadata <- DT::renderDT({
    df <- pool |>
      tbl("metadata") |>
      collect() |>
      dplyr::select(-id)
    datatable(
      df,
      colnames = c(
        "Cell Line", "Drug", "BioProject", "Contrast", "Dose",
        "Time (hr)", "Batch", "Mutation", "Comment", "Description",
        "Epigenetic Class", "Tissue", "Disease"
      ),
      options = list(columnDefs = list(
        list(
          targets = 4,
          render = JS(
            "function(data, type, row, meta) {",
            "return type === 'display' && data.length > 12 ?",
            "'<span title=\"' + data + '\">' + data.substr(0, 12) + '...</span>' : data;",
            "}"
          )
        ),
        list(
          targets = 10,
          width = "50%"
        )
      )),
      callback = JS("table.page(3).draw(false);")
    )
  })
}
