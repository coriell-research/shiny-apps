library(shiny)
library(shinythemes)


ui <- fluidPage(
  theme = shinytheme("cerulean"),
  titlePanel("Explore Normalization Methods"),
  sidebarLayout(
    sidebarPanel(
      width = 2,
      h3("Step 1. Import and Load Data"),
      p("The SummarizedExperiment object must contain an assay with raw 
              counts named 'counts' and colData with a column named 'group' 
              specifying the grouping factor for the samples. Other metadata 
              columns can be present but 'group' is necessary."),
      fileInput("se",
        label = "Import SummarizedExperiment",
        accept = c(".rds")
      ),
      actionButton("load",
        label = "Load Dataset",
        width = "100%"
      ),
      h3("Step 2. Filtering Parameters"),
      numericInput("min_count",
        label = "Min Count per Group",
        value = 10,
        min = 0,
        step = 1
      ),
      numericInput("min_total_count",
        label = "Minimum Total Count per Feature",
        value = 15,
        min = 0,
        step = 1
      ),
      numericInput("min_prop",
        label = "Minimum Proportion per Group",
        value = 0.7,
        min = 0, 
        max = 1,
        step = 0.1
      ),
      h3("Step 3. Normalization Parameters"),
      numericInput("percentile",
        label = "Percentile for UQ Normalization",
        value = 0.7,
        min = 0,
        max = 1,
        step = 0.1
      ),
      numericInput("pseudocount",
        label = "Pseudocount for logCPM calculations",
        value = 2,
        min = 0,
        step = 1
      ),
      selectInput("reference_col",
        label = "Enter reference column for TMM normalization (optional)",
        choices = NULL
      ),
      fileInput("control_genes",
        label = "Upload a list of control genes for RUVg",
        accept = c(".rds")
      ),
      numericInput("K",
        label = "Enter Number Factors (k) for RUVg",
        value = 2,
        min = 1,
        step = 1
      ),
      h3("Step 4. Run"),
      actionButton("run",
        label = "Run filtering and normalization",
        width = "100%"
      )
    ),

    mainPanel(
      tabsetPanel(
        tabPanel(
          "Filtering",
          selectInput("hist_sample",
            label = "Select sample",
            choices = NULL
          ),
          numericInput("hist_pseudocount",
            label = "Pseudocount after filtering",
            value = 0,
            min = 0,
            step = 1
          ),
          sliderInput("hist_bins",
            label = "Enter number of bins",
            value = 30, 
            min = 5, 
            max = 100
          ),
          fluidRow(
            column(
              width = 6,
              h3("Raw Data"),
              plotOutput("raw_hist")
            ),
            column(
              width = 6,
              h3("Filtered Data"),
              plotOutput("filt_hist")
            )
          ),
        ),
        tabPanel(
          "RLE",
          fluidRow(
            column(
              width = 6,
              selectInput("rle_norm_method",
                label = "Select normalization method",
                choices = c(list("TMM" = "logTMM", "RLE" = "logRLE", "UQ" = "logUQ", "RUVg" = "logRUVg", "QS" = "logQS", "LibrarySize" = "logLibrarySize")),
                selected = "TMM"
              ),
              selectInput("rle_fill_by",
                label = "Fill by",
                choices = NULL
              ),
              numericInput("rle_outlier_shape",
                label = "Outlier shape",
                value = NA,
                min = 0, 
                max = 25,
                step = 1
              ),
              numericInput("rle_outlier_alpha",
                label = "Outlier alpha",
                value = 0.5,
                min = 0, 
                max = 1,
                step = 0.1
              ),
              plotOutput("rle")
            ),
            column(
              width = 6,
              selectInput("rle_norm_method2",
                label = "Select normalization method",
                choices = c(list("TMM" = "logTMM", "RLE" = "logRLE", "UQ" = "logUQ", "RUVg" = "logRUVg", "QS" = "logQS", "LibrarySize" = "logLibrarySize")),
                selected = "TMM"
              ),
              selectInput("rle_fill_by2",
                label = "Fill by",
                choices = NULL
              ),
              numericInput("rle_outlier_shape2",
                label = "Outlier shape",
                value = NA,
                min = 0,
                max = 25,
                step = 1
              ),
              numericInput("rle_outlier_alpha2",
                label = "Outlier alpha",
                value = 0.5,
                min = 0,
                max = 1,
                step = 0.1
              ),
              plotOutput("rle2")
            )
          )
        ),
        tabPanel(
          "MA",
          fluidRow(
            column(
              width = 6,
              selectInput("ma_norm_method",
                label = "Select normalization method",
                choices = c("TMM", "RLE", "UQ", "RUVg", "QS", "LibrarySize"),
                selected = "TMM"
              ),
              selectInput("ma_sample1",
                label = "Sample 1",
                choices = NULL
              ),
              selectInput("ma_sample2",
                label = "Sample 2",
                choices = NULL
              ),
              checkboxInput("ma_smooth",
                label = "Smooth scatter",
                value = FALSE
              ),
              checkboxInput("ma_loess",
                label = "Plot lowess line",
                value = FALSE
              ),
              plotOutput("ma")
            ),
            column(
              width = 6,
              selectInput("ma_norm_method2",
                label = "Select normalization method",
                choices = c("TMM", "RLE", "UQ", "RUVg", "QS", "LibrarySize"),
                selected = "TMM"
              ),
              selectInput("ma_sample12",
                label = "Sample 1",
                choices = NULL
              ),
              selectInput("ma_sample22",
                label = "Sample 2",
                choices = NULL
              ),
              checkboxInput("ma_smooth2",
                label = "Smooth scatter",
                value = FALSE
              ),
              checkboxInput("ma_loess2",
                label = "Plot lowess line",
                value = FALSE
              ),
              plotOutput("ma2")
            )
          )
        ),
        tabPanel(
          "Scatter",
          fluidRow(
            column(
              width = 6,
              selectInput("scatter_norm_method",
                label = "Select normalization method",
                choices = c("TMM", "RLE", "UQ", "RUVg", "QS", "LibrarySize"),
                selected = "TMM"
              ),
              selectInput("scatter_sample1",
                label = "Sample 1",
                choices = NULL
              ),
              selectInput("scatter_sample2",
                label = "Sample 2",
                choices = NULL
              ),
              numericInput("scatter_pt_alpha",
                label = "Point Alpha",
                value = 0.5,
                min = 0,
                max = 1,
                step = 0.1
              ),
              checkboxInput("scatter_log",
                label = "Log scale axis",
                value = TRUE
              ),
              plotOutput("scatter")
            ),
            column(
              width = 6,
              selectInput("scatter_norm_method2",
                label = "Select normalization method",
                choices = c("TMM", "RLE", "UQ", "RUVg", "QS", "LibrarySize"),
                selected = "TMM"
              ),
              selectInput("scatter_sample12",
                label = "Sample 1",
                choices = NULL
              ),
              selectInput("scatter_sample22",
                label = "Sample 2",
                choices = NULL
              ),
              numericInput("scatter_pt_alpha2",
                label = "Point Alpha",
                value = 0.5,
                min = 0,
                max = 1,
                step = 0.1
              ),
              checkboxInput("scatter_log2",
                label = "Log scale axis",
                value = TRUE
              ),
              plotOutput("scatter2")
            )
          )
        ),
        tabPanel(
          "PCA",
          fluidRow(
            column(
              width = 6,
              selectInput("pca_norm_method",
                label = "Select normalization method",
                choices = c(list("TMM" = "logTMM", "RLE" = "logRLE", "UQ" = "logUQ", "RUVg" = "logRUVg", "QS" = "logQS", "LibrarySize" = "logLibrarySize")),
                selected = "TMM"
              ),
              checkboxInput("pca_scale",
                label = "Scale data",
                value = TRUE
              ),
              checkboxInput("pca_center",
                label = "Center data",
                value = TRUE
              ),
              selectInput("pca_component1",
                label = "X-axis",
                choices = c("PC1", "PC2", "PC3", "PC4", "PC5"),
                selected = "PC1"
              ),
              selectInput("pca_component2",
                label = "Y-axis",
                choices = c("PC1", "PC2", "PC3", "PC4", "PC5"),
                selected = "PC2"
              ),
              selectInput("pca_color_by",
                label = "Color by",
                choices = NULL
              ),
              selectInput("pca_shape_by",
                label = "Shape by",
                choices = NULL
              ),
              plotOutput("pca")
            ),
            column(
              width = 6,
              selectInput("pca_norm_method2",
                label = "Select normalization method",
                choices = c(list("TMM" = "logTMM", "RLE" = "logRLE", "UQ" = "logUQ", "RUVg" = "logRUVg", "QS" = "logQS", "LibrarySize" = "logLibrarySize")),
                selected = "TMM"
              ),
              checkboxInput("pca_scale2",
                label = "Scale data",
                value = TRUE
              ),
              checkboxInput("pca_center2",
                label = "Center data",
                value = TRUE
              ),
              selectInput("pca_component12",
                label = "X-axis",
                choices = c("PC1", "PC2", "PC3", "PC4", "PC5"),
                selected = "PC1"
              ),
              selectInput("pca_component22",
                label = "Y-axis",
                choices = c("PC1", "PC2", "PC3", "PC4", "PC5"),
                selected = "PC2"
              ),
              selectInput("pca_color_by2",
                label = "Color by",
                choices = NULL
              ),
              selectInput("pca_shape_by2",
                label = "Shape by",
                choices = NULL
              ),
              plotOutput("pca2")
            )
          )
        ),
        tabPanel(
          "Corr",
          column(
            width = 6,
            selectInput("cor_norm_method",
                        label = "Select normalization method",
                        choices = c("TMM", "RLE", "UQ", "RUVg", "QS", "LibrarySize"),
                        selected = "TMM"),
            selectizeInput("cor_samples",
                           label = "Select samples",
                           choices = NULL,
                           multiple = TRUE,
                           size = 10),
            selectInput("cor_viz_method",
                        label = "Visualization method",
                        choices = c("circle", "square", "ellipse", "number", "pie", "shade", "color")),
            plotOutput("corplot")
          ),
          column(
            width = 6,
            selectInput("cor_norm_method2",
                        label = "Select normalization method",
                        choices = c("TMM", "RLE", "UQ", "RUVg", "QS", "LibrarySize"),
                        selected = "TMM"),
            selectizeInput("cor_samples2",
                           label = "Select samples",
                           choices = NULL,
                           multiple = TRUE,
                           size = 10),
            selectInput("cor_viz_method2",
                        label = "Visualization method",
                        choices = c("circle", "square", "ellipse", "number", "pie", "shade", "color")),
            plotOutput("corplot2")
          ),
        ),
        tabPanel(
          "Heatmap",
          column(
            width = 6,
            selectInput("hm_norm_method",
                        label = "Select normalization method",
                        choices = c(list("TMM" = "logTMM", "RLE" = "logRLE", "UQ" = "logUQ", "RUVg" = "logRUVg", "QS" = "logQS", "LibrarySize" = "logLibrarySize")),
                        selected = "TMM"
            ),
            selectizeInput("hm_samples",
                           label = "Select samples",
                           choices = NULL,
                           multiple = TRUE,
                           size = 10),
            numericInput("hm_features",
                         label = "Top N variable features",
                         value = 1000,
                         min = 100, 
                         max = 10000),
            selectInput("hm_cd_rows",
                        label = "Cluster distance rows",
                        choices = c("euclidean", "correlation", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                        selected = "euclidean"),
            selectInput("hm_cd_cols",
                        label = "Cluster distance cols",
                        choices = c("euclidean", "correlation", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                        selected = "euclidean"),
            selectInput("hm_clust_method",
                        label = "Clustering method",
                        choices = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"),
                        selected = "complete"),
            plotOutput("hm")
          ),
          column(
            width = 6,
            selectInput("hm_norm_method2",
                        label = "Select normalization method",
                        choices = c(list("TMM" = "logTMM", "RLE" = "logRLE", "UQ" = "logUQ", "RUVg" = "logRUVg", "QS" = "logQS", "LibrarySize" = "logLibrarySize")),
                        selected = "TMM"
            ),
            selectizeInput("hm_samples2",
                           label = "Select samples",
                           choices = NULL,
                           multiple = TRUE,
                           size = 10),
            numericInput("hm_features2",
                         label = "Top N variable features",
                         value = 1000,
                         min = 100, 
                         max = 10000),
            selectInput("hm_cd_rows2",
                        label = "Cluster distance rows",
                        choices = c("euclidean", "correlation", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                        selected = "euclidean"),
            selectInput("hm_cd_cols2",
                        label = "Cluster distance cols",
                        choices = c("euclidean", "correlation", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                        selected = "euclidean"),
            selectInput("hm_clust_method2",
                        label = "Clustering method",
                        choices = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"),
                        selected = "complete"),
            plotOutput("hm2")
          )
        )
      )
    )
  )
)
