library(shiny)
library(plotly)
library(shinythemes)
library(shinyWidgets)


# load choices for select inputs -----------------------------------------------
all_samples <- readRDS("data/all-samples.rds")
all_projects <- readRDS("data/all-projects.rds")
all_cells <- readRDS("data/all-cells.rds")
all_drugs <- readRDS("data/all-drugs.rds")


# Define UI for application ----------------------------------------------------
shinyUI(
    fluidPage(
        theme = shinytheme("cerulean"),
        titlePanel("Explore GEO RNA-Seq"),
        sidebarLayout(
            sidebarPanel(
                radioButtons("count_matrix", 
                             label = "Select Count Matrix",
                             choices = c("Batch Corrected Gene Counts", "Batch Corrected RE Counts", "Raw Gene Counts", "Raw RE Counts")),
                textInput("remove_regex",
                          label = "Remove features with regex",
                          value = "Simple_repeat|Unknown",
                          placeholder = "Simple_repeat|Unknown"),
                pickerInput("samples",
                            label = "Select Samples",
                            choices = all_samples,
                            selected = all_samples,
                            multiple = TRUE,
                            options = list(`actions-box` = TRUE)),
                pickerInput("projects",
                            label = "Select Experiments",
                            choices = all_projects,
                            selected = all_projects,
                            multiple = TRUE,
                            options = list(`actions-box` = TRUE)),
                pickerInput("cells",
                            label = "Select Cell Types",
                            choices = all_cells,
                            selected = all_cells,
                            multiple = TRUE,
                            options = list(`actions-box` = TRUE)),
                pickerInput("drugs",
                            label = "Select Drugs",
                            choices = all_drugs,
                            selected = all_drugs,
                            multiple = TRUE,
                            options = list(`actions-box` = TRUE)),
                checkboxInput("center_data",
                              label = "Center data in PCA?",
                              value = TRUE),
                checkboxInput("scale_data",
                              label = "Scale data in PCA?",
                              value = FALSE),
                numericInput("remove_var",
                             label = "Remove this proportion of variables based on low variance",
                             value = NULL,
                             min = 0,
                             max = 0.99,
                             step = 0.1),
                h3("Run PCA and Plot"),
                actionBttn("run",
                           label = "Calculate and Plot")
                ),
            mainPanel(
                tabsetPanel(
                    tabPanel("3D PCA",
                             plotlyOutput("pca", width = "100%", height = "700px"),
                             column(4,
                                    selectInput("pca_xvar",
                                                label = "x-axis",
                                                choices = paste0("PC", 1:25),
                                                selected = "PC1"),
                                    selectInput("pca_yvar",
                                                label = "y-axis",
                                                choices = paste0("PC", 1:25),
                                                selected = "PC2"),
                                    selectInput("pca_zvar",
                                                label = "z-axis",
                                                choices = paste0("PC", 1:25),
                                                selected = "PC3")
                                    ),
                             column(4,
                                    selectInput("pca_color",
                                                label = "Color points by:",
                                                choices = c("BioProject" = "bio_project", "Drug" = "drug", "Cell Type" = "cell_type"),                                                            ),
                                                selected = "Drug")
                             ),
                    tabPanel("Biplot",
                             plotOutput("biplot"),
                             column(4,
                                    selectInput("biplot_x",
                                                label = "x-axis",
                                                choices = paste0("PC", 1:25),
                                                selected = "PC1"),
                                    selectInput("biplot_y",
                                                label = "y-axis",
                                                choices = paste0("PC", 1:25),
                                                selected = "PC2"),
                                    checkboxInput("show_loading",
                                                  label = "Show loadings",
                                                  value = FALSE)
                                    ),
                             column(4,
                                    selectInput("biplot_color",
                                                label = "Color points by:",
                                                choices = c("None" = "None",
                                                            "Project" = "bio_project",
                                                            "Cell Type" = "cell_type",
                                                            "Drug" = "drug")),
                                    selectInput("biplot_shape",
                                                label = "Shape points by:",
                                                choices = c("None" = "None",
                                                            "Project" = "bio_project",
                                                            "Cell Type" = "cell_type",
                                                            "Drug" = "drug")),
                                    checkboxInput("biplot_labs",
                                                  label = "Show point labels",
                                                  value = FALSE)
                                    )
                             ),
                    tabPanel("Scree Plot",
                             plotOutput("screeplot")
                             ),
                    tabPanel("Pairs PLot",
                             plotOutput("pairsplots", width = "95%", height = "800px")
                             ),
                    tabPanel("Loadings Plot",
                             plotOutput("loadings", width = "95%", height = "700px")
                             )
                )
            )
        ))
)
