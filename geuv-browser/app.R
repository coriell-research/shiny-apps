library(shiny)
library(shinythemes)
library(shinyWidgets)
library(pool)

library(forcats)
library(ggplot2)
library(dplyr)
library(dbplyr)
library(here)


# set up global data -----------------------------------------------------------
# create database connection
pool <- dbPool(
  drv = RSQLite::SQLite(),
  dbname = here("data", "data.sqlite")
)

# Read in input selections
load(here("data", "input-selections.RData"))

# set up ui --------------------------------------------------------------------

ui <- fluidPage(
  navbarPage(
    title = "Coriell GEUVADIS Expression Browser",
    theme = shinytheme("flatly"),
    tabPanel(
      "Browser",
      pageWithSidebar(
        headerPanel("Select Data"),
        sidebarPanel(
          width = 2,
          selectizeInput(
            inputId = "gene",
            label = "Enter Gene Symbol",
            choices = NULL,
            multiple = FALSE
          ),
          pickerInput(
            inputId = "individuals",
            label = "Individual(s)",
            choices = individuals,
            selected = individuals,
            options = list(`actions-box` = TRUE, size = 10),
            multiple = TRUE
          ),
          checkboxGroupInput(
            inputId = "populations",
            label = "Populations(s)",
            choices = populations,
            selected = populations
          ),
          checkboxGroupInput(
            inputId = "sexes",
            label = "Sex",
            choices = sexes,
            selected = sexes
          ),
          selectInput(
            inputId = "metric",
            label = "Count Metric",
            choices = metrics,
            selected = "CPM",
            multiple = FALSE
          )
        ),
        mainPanel = mainPanel(
          tabsetPanel(
            id = "main",
            tabPanel(
              "Table", DT::dataTableOutput("table"),
              downloadButton("downloadData", "Download")
            ),
            tabPanel(
              "Expression Histogram",
              sliderInput("bins",
                label = "Select number of bins",
                value = 30,
                min = 1,
                max = 462
              ),
              selectizeInput("individual",
                "Show expression for:",
                choices = individuals,
                options = list(
                  placeholder = "Select an indivual sample",
                  onInitialize = I('function() { this.setValue(""); }')
                )
              ),
              plotOutput("hist"),
              h3("Distribution Summary"),
              verbatimTextOutput("summary")
            )
          )
        )
      )
    ),
    tabPanel(
      "About",
      includeMarkdown(here("doc", "about.md"))
    )
  )
)


# set up server functions ------------------------------------------------------

server <- function(input, output, session) {

  # create serverside selectize for gene list
  updateSelectizeInput(session, "gene", choices = genes, selected = "BRCA1", server = TRUE)

  # set up base reactive table used in most tabsets -----------------
  d <- reactive({
    pool %>%
      tbl("data") %>%
      filter(
        gene_name == !!input$gene,
        metric == !!input$metric,
        individual %in% !!input$individuals,
        population %in% !!input$populations,
        sex %in% !!input$sexes
      )
  })

  # datatable tabset ------------------------------------------------
  output$table <- DT::renderDataTable({
    d() %>%
      collect() %>%
      as.data.frame() %>% 
      mutate(individual = paste0('<a href="https://www.coriell.org/0/Sections/Search/Sample_Detail.aspx?Ref=', individual, '&Product=DNA">', individual, '</a>'),
             value = round(value, 2)) %>% 
      arrange(desc(value))
  }, 
  escape = c('sex', 'population', 'gene_name', 'gene_id', 'seqnames', 'start', 'end', 'width', 'strand', 'gene_type', 'metric', 'value')
  )

  output$downloadData <- downloadHandler(
    filename = function() {
      paste("Coriell-GEUV-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(d() %>% collect() %>% as.data.frame(), file)
    }
  )

  # histogram tabset ------------------------------------------------
  individual_expr <- reactive({
    d() %>%
      filter(individual == !!input$individual) %>%
      pull(value)
  })

  output$hist <- renderPlot({
    d() %>%
      ggplot(aes(x = value)) +
      geom_histogram(fill = "steelblue", color = "black", bins = input$bins) +
      geom_vline(xintercept = individual_expr(), linetype = 2, color = "red", size = 3) +
      theme_light() +
      labs(
        title = paste("Distribution of", input$gene, "Expression"),
        x = input$metric,
        y = "Count"
      )
  })

  output$summary <- renderPrint({
    summary(pull(.data = collect(d()), var = value))
  })
}


# run app ----------------------------------------------------------------------

shinyApp(ui = ui, server = server)
