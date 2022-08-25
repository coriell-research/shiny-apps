library(shiny)
library(shinythemes)
library(shinyWidgets)
library(SummarizedExperiment)
library(plotly)


se <- readRDS("data/se.rds")
pathways <- readRDS("data/pathways.rds")
id_choices <- sort(unique(colData(se)$id))
contrast_choices <- sort(unique(colData(se)$contrast))
project_choices <- sort(unique(colData(se)$experiment))
cell_choices <- sort(unique(colData(se)$cell_line))
drug_choices <- sort(unique(colData(se)$drug))
tissue_choices <- sort(unique(colData(se)$tissue))
disease_choices <- sort(unique(colData(se)$disease))
epiclass_choices <- sort(unique(colData(se)$epigenetic_class))
pathway_choices <- c("ALL GENES", sort(names(pathways)))

# ------------------------------------------------------------------------------
ui <- fluidPage(
  theme = shinytheme("yeti"),
  titlePanel("RNA-seq Meta-Analysis"),
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
        choices = c("LogFC" = "lfc", "FDR" = "fdr"),
        selected = "lfc"
      ),
      pickerInput("contrasts",
        label = "Conditions",
        choices = contrast_choices,
        selected = contrast_choices,
        multiple = TRUE,
        options = list(`actions-box` = TRUE),
        width = "100%"
      ),
      pickerInput("experiments",
        label = "Experiments",
        choices = project_choices,
        selected = project_choices,
        multiple = TRUE,
        options = list(`actions-box` = TRUE),
        width = "100%"
      ),
      pickerInput("cell_lines",
        label = "Cell lines",
        choices = cell_choices,
        selected = cell_choices,
        multiple = TRUE,
        options = list(`actions-box` = TRUE),
        width = "100%"
      ),
      pickerInput("drugs",
        label = "Drugs",
        choices = drug_choices,
        selected = drug_choices,
        multiple = TRUE,
        options = list(`actions-box` = TRUE),
        width = "100%"
      ),
      pickerInput("tissues",
        label = "Tissues",
        choices = tissue_choices,
        selected = tissue_choices,
        multiple = TRUE,
        options = list(`actions-box` = TRUE),
        width = "100%"
      ),
      pickerInput("diseases",
        label = "Diseases",
        choices = disease_choices,
        selected = disease_choices,
        multiple = TRUE,
        options = list(`actions-box` = TRUE),
        width = "100%"
      ),
      pickerInput("epigenetic",
        label = "Epigenetic Class",
        choices = epiclass_choices,
        selected = epiclass_choices,
        multiple = TRUE,
        options = list(`actions-box` = TRUE),
        width = "100%"
      ),
      pickerInput("pathway",
                  label = "Gene Set",
                  choices = pathway_choices,
                  selected = "ALL GENES",
                  multiple = FALSE,
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
      numericInput("components",
                   label = "Components",
                   value = 10,
                   min = 2,
                   max = 100
      ),
      pickerInput("algorithm",
                  label = "Algorithm",
                  choices = c("Exact", "Irlba", "Random", "Auto"),
                  selected = "Auto",
                  multiple = FALSE,
                  options = list(`actions-box` = TRUE),
                  width = "100%"
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
          fluidRow(
            column(
              4,
              selectInput(
                "col",
                label = "Color by:",
                choices = c(
                  "Tissue" = "tissue", "Cell Line" = "cell_line",
                  "Treatment" = "drug", "Disease" = "disease",
                  "Experiment" = "experiment",
                  "Epigenetic Class" = "epigenetic_class"
                ),
                selected = "tissue"
              )
            ),
            column(
              4,
              selectInput(
                "x_axis",
                label = "x-axis",
                choices = paste0("PC", 1:10),
                selected = "PC1"
              )
            ),
            column(
              4,
              selectInput(
                "y_axis",
                label = "y-axis",
                choices = paste0("PC", 1:10),
                selected = "PC2"
              )
            )
          ),
          plotlyOutput("biplot"),
          DT::dataTableOutput("table")
        ),
        tabPanel(
          "Meta-Vote",
          fluidRow(
            column(
              4,
              numericInput("meta_pval", label = "FDR Cutoff", value = 0.1, min = 0, max = 1, step = 0.1)
            ),
            column(
              4,
              numericInput("meta_lfc", label = "LogFC Cutoff", value = 0, min = 0, max = Inf, step = 1)
            ),
            column(
              4,
              numericInput("meta_prop", label = "Proportion Exp. DE", value = 0.1, min = 0, max = 1, step = 0.1)
            )
          ),
          plotOutput("metavolcano"),
          DT::dataTableOutput("metatable")
        ),
        tabPanel(
          "Differential Expression",
          fluidRow(
            column(6, plotOutput("volcano")),
            column(6, plotOutput("ma"))
          ),
          fluidRow(
            column(6, selectizeInput("id", label = "Contrast", choices = id_choices, multiple = FALSE, selected = NULL, width = "100%")),
            column(3, pickerInput("type", label = "Feature Type", choices = c("gene", "re", "both"), multiple = FALSE, selected = "gene")),
            column(3, numericInput("fdr", label = "FDR Cutoff", value = 0.1, min = 0, max = 1, step = 0.1))
          ),
          DT::dataTableOutput("de_tbl")
        ),
        tabPanel(
          "Expression Ranking",
          numericInput("ranking_fdr", label = "FDR Cutoff", value = 0.1, min = 0, max = 1, step = 0.1),
          DT::dataTableOutput("ranking")
        ),
        tabPanel(
          "Metadata",
          DT::dataTableOutput("metadata")
        )
      )
    )
  )
)
