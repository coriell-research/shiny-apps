library(shiny)
library(shinyWidgets)
library(here)
library(data.table)
library(SummarizedExperiment)
library(plotly)
library(PCAtools)
library(magrittr)
library(patchwork)
library(coriell)


dge_se <- readRDS(here("data", "dge_se.rds"))
dre_se <- readRDS(here("data", "dre_se.rds"))
dge_res <- fread(here("data", "dge-all-experiments.csv.gz"))
dre_res <- fread(here("data", "dre-all-experiments.csv.gz"))
meta_dt <- as.data.table(x = colData(dge_se), keep.rownames = "id")
res_dt <- rbindlist(list("dge" = dge_res, "dre" = dre_res))
features <- sort(unique(res_dt$feature_id))

# ------------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("GEO Cancer Meta-Analysis"),
  sidebarLayout(
    sidebarPanel(
      width = 2,
      radioButtons("dataset",
        label = "Datasets",
        choices = c("Protein Coding Genes", "Repetitive Elements"),
        selected = "Protein Coding Genes"
      ),
      radioButtons("data",
        label = "Data",
        choices = c("Log2FC" = "lfc", "FDR" = "fdr", "Rank" = "rank"),
        selected = "lfc"
      ),
      pickerInput("contrasts",
        label = "Conditions",
        choices = sort(unique(colData(dge_se)$contrast)),
        selected = unique(colData(dge_se)$contrast),
        multiple = TRUE,
        options = list(`actions-box` = TRUE),
        width = "100%"
      ),
      pickerInput("experiments",
        label = "Experiments",
        choices = sort(unique(colData(dge_se)$experiment)),
        selected = unique(colData(dge_se)$experiment),
        multiple = TRUE,
        options = list(`actions-box` = TRUE),
        width = "100%"
      ),
      pickerInput("cell_lines",
        label = "Cell lines",
        choices = sort(unique(colData(dge_se)$cell_line)),
        selected = unique(colData(dge_se)$cell_line),
        multiple = TRUE,
        options = list(`actions-box` = TRUE),
        width = "100%"
      ),
      pickerInput("drugs",
        label = "Drugs",
        choices = sort(unique(colData(dge_se)$drug)),
        selected = unique(colData(dge_se)$drug),
        multiple = TRUE,
        options = list(`actions-box` = TRUE),
        width = "100%"
      ),
      pickerInput("tissues",
        label = "Tissues",
        choices = sort(unique(colData(dge_se)$tissue)),
        selected = unique(colData(dge_se)$tissue),
        multiple = TRUE,
        options = list(`actions-box` = TRUE),
        width = "100%"
      ),
      pickerInput("diseases",
        label = "Diseases",
        choices = sort(unique(colData(dge_se)$disease)),
        selected = unique(colData(dge_se)$disease),
        multiple = TRUE,
        options = list(`actions-box` = TRUE),
        width = "100%"
      ),
      pickerInput("epigenetic",
        label = "Epigenetic Class",
        choices = sort(unique(colData(dge_se)$epigenetic_class)),
        selected = unique(colData(dge_se)$epigenetic_class),
        multiple = TRUE,
        options = list(`actions-box` = TRUE),
        width = "100%"
        ),
      h4("PCA Params:"),
      numericInput("var",
        label = "Remove Variance",
        value = 0.1,
        min = 0,
        max = 0.95,
        step = 0.05
      ),
      checkboxInput("scale",
        label = "Scale Data",
        value = FALSE
      ),
      checkboxInput("center",
        label = "Center Data",
        value = TRUE
      ),
      checkboxInput("complete",
        label = "Use Complete Cases",
        value = FALSE
      ),
      actionBttn("run",
        label = "Run"
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "PCA Biplot",
          selectInput("col",
            label = "Color by:",
            choices = c("tissue", "cell_line", "drug", "disease", "experiment", "epigenetic_class"),
            selected = "tissue"
          ),
          selectInput("x_axis",
            label = "x-axis",
            choices = paste0("PC", 1:15),
            selected = "PC1"
          ),
          selectInput("y_axis",
            label = "y-axis",
            choices = paste0("PC", 1:15),
            selected = "PC2"
          ),
          plotlyOutput("biplot"),
          DT::dataTableOutput("table")
        ),
        tabPanel(
          "PCA Loadings",
          plotOutput("loadings",
            width = "1600px",
            height = "900px"
          )
        ),
        tabPanel(
          "Scree Plot",
          plotOutput("scree",
            width = "1600px",
            height = "900px"
          )
        ),
        tabPanel(
          "Diff Exp.",
          pickerInput("patch_con",
                      label = "Contrasts",
                      choices = sort(unique(colData(dge_se)$contrast)),
                      selected = "HH1.25uM_vs_DMSO_96hr",
                      multiple = FALSE
          ),
          numericInput("patch_fdr",
                       label = "FDR cutoff",
                       value = 0.1,
                       min = 0,
                       max = 1,
                       step = 0.05),
          numericInput("patch_fc",
                       label = "FC cutoff",
                       value = 1.5,
                       min = 0, 
                       max = Inf),
          plotOutput("patch")
        ),
        tabPanel(
          "Metadata",
          DT::dataTableOutput("metadata")
        ),
        tabPanel(
          "Features",
          selectizeInput("features",
            label = "Feature Name",
            choices = NULL,
            multiple = TRUE,
            options = list(maxOptions = length(features))
          ),
          DT::dataTableOutput("features")
        ),
        tabPanel(
          "Meta-Vote",
          numericInput("metav_lfc",
                       label = "LogFC cutoff",
                       value = 0,
                       min = 0,
                       max = 25),
          numericInput("metav_pval",
                       label = "FDR cutoff",
                       value = 0.1,
                       min = 0,
                       max = 1),
          numericInput("metav_prop",
                       label = "Proportion of experiments gene is DE",
                       value = 0.1,
                       min = 0,
                       max = 1),
          checkboxInput("metav_common",
                        label = "Use only genes common to all experiments",
                        value = FALSE),
          plotlyOutput("metavolcano"),
        ),
        tabPanel(
          "Meta-Pcombine",
          selectInput("metap_method",
                      label = "P-value combination method",
                      choices = c("fisher", "pearson", "tippet", "wilkinson"),
                      selected = "fisher",
                      multiple = FALSE),
          selectInput("metap_fun",
                      label = "Function to combine logFC values",
                      choices = c("mean", "median", "max", "min"),
                      selected = "mean",
                      multiple = FALSE),
          checkboxInput("metap_common",
                        label = "Use only genes common to all experiments",
                        value = FALSE),
          plotlyOutput("metap_volcano"),
          numericInput("metap_fdr",
                       label = "FDR cutoff",
                       value = 0.05,
                       min = 0,
                       max = 1),
          numericInput("metap_lfc",
                       label = "LogFC cutoff",
                       value = 0,
                       min = 0,
                       max = 100)
        )
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  updateSelectizeInput(session, "features", choices = features, selected = NULL, server = TRUE)
  
  data <- eventReactive(input$run, {

    # Select the dataset based on user input
    d <- dge_se
    d2 <- dge_res
    if (input$dataset == "Repetitive Elements") {
      d <- dre_se
      d2 <- dre_res
    }

    # filter dataset based on input
    d_filtered <- d[, d$contrast %chin% input$contrasts &
      d$experiment %chin% input$experiments &
      d$cell_line %chin% input$cell_lines &
      d$drug %chin% input$drugs &
      d$tissue %chin% input$tissues &
      d$disease %chin% input$diseases &
      d$epigenetic_class %chin% input$epigenetic]
    
    # get the ids for the contrasts in the filtered list
    ids <- colnames(d_filtered)
    dge_list <- split(d2[id %chin% ids], by = "id")
    
    M <- assay(d_filtered, input$data)
    if (input$complete) {
      M[M == 0] <- NA
      M <- na.omit(M)
    }

    # perform PCA using PCAtools
    pca_res <- pca(
      mat = M,
      metadata = as.data.frame(colData(d_filtered)),
      scale = input$scale,
      center = input$center,
      removeVar = input$var
    )

    # Return PCA, filtered matrix and dge list as reactive objects
    list("PCA" = pca_res, "M" = M, "dge_list" = dge_list)
  })

  # Create single reactive data.table for use in plotting functions
  dt <- reactive({
    pca_dt <- as.data.table(data()[["PCA"]]$rotated, keep.rownames = "sample_name")
    meta_dt <- as.data.table(data()[["PCA"]]$metadata, keep.rownames = "sample_name")
    dt <- merge.data.table(x = pca_dt, y = meta_dt, by = "sample_name", all.x = TRUE)
    dt
  })

  # PCA biplot
  output$biplot <- renderPlotly({
    dt <- dt()
    plot_ly(
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
    ) %>%
      add_markers(size = 180) %>%
      layout(
        dragmode = "lasso",
        showlegend = FALSE,
        xaxis = list(title = paste(input$x_axis)),
        yaxis = list(title = paste(input$y_axis))
      )
  })

  # Table
  output$table <- DT::renderDT({
    d <- event_data("plotly_selected")
    if (!is.null(d)) {
      xmin <- min(d$x) - 0.000001
      xmax <- max(d$x) + 0.000001
      ymin <- min(d$y) - 0.000001
      ymax <- max(d$y) + 0.000001
      dt <- dt()
      dt[
        between(get(input$x_axis), xmin, xmax, incbounds = TRUE) &
          between(get(input$y_axis), ymin, ymax, incbounds = TRUE),
        .(experiment, tissue, cell_line, disease, drug, dose, time_hr, epigenetic_class)
      ]
    }
  })

  # Loadings plot
  output$loadings <- renderPlot({
    plotloadings(data()[["PCA"]],
      components = 1:5,
      labSize = 6,
      drawConnectors = TRUE
    )
  })

  # Screeplot
  output$scree <- renderPlot({
    screeplot(data()[["PCA"]],
      components = 1:15,
      drawCumulativeSumLine = TRUE,
      drawCumulativeSumPoints = TRUE
    )
  })
  
  # Differential expression plots
  output$patch <- renderPlot({
    d <- dge_res
    if (input$dataset == "Repetitive Elements") {
      d <- dre_res
    }
    
    data <- d[contrast == input$patch_con]
    vol <- plot_volcano(data, fdr = input$patch_fdr, lfc = log2(input$patch_fc))
    md <- plot_md(data, fdr = input$patch_fdr, lfc = log2(input$patch_fc))
    (vol | md)
  })
  
  # Metadata table
  output$metadata <- DT::renderDataTable(meta_dt)
  
  # Differential expression data
  output$features <- DT::renderDataTable({res_dt[feature_id %chin% input$features]})
  
  # Meta-Analysis vote counting
  meta_votecount <- reactive({
    meta_vote(
      exp_list = data()[["dge_list"]],
      lfc = input$metav_lfc,
      pval = input$metav_pval,
      meta_prop = input$metav_prop,
      all_common = input$metav_common
      )
  })
  
  meta_pcomb <- reactive({
    LFC_FUN  <-  switch(input$metap_fun,
                        mean = mean,
                        median = median,
                        max = max,
                        min = min)
    meta_pcombine(
      exp_list = data()[["dge_list"]],
      method = input$metap_method,
      lfc_fun = LFC_FUN,
      all_common = input$metap_common,
      na.rm = TRUE
    )
  })
  
  # meta-volcano (vote) plot
  output$metavolcano <- renderPlotly({
    ggplotly(plot_metavolcano(meta_votecount()))
  })
  
  # meta-volcano pcombine plot
  output$metap_volcano <- renderPlotly({
    ggplotly(
      plot_volcano(df = meta_pcomb(), x = meta_lfc, y = meta_p, fdr = input$metap_fdr, lfc = input$metap_lfc)
    )
  })
}

# Run the application
shinyApp(ui = ui, server = server)
