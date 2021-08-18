library(shiny)
library(shinythemes)
library(tidyverse)
library(readxl)
library(lubridate)
library(openxlsx)

source("helpers.R")
ssrs_cols <- readRDS("data/ssrs-cols.rds")


# Define UI for application ----------------------------------------------------
ui <- fluidPage(
  theme = shinytheme("cosmo"),
  titlePanel("Process NHGRI Orders"),
  sidebarLayout(
    sidebarPanel(
      h3("Step 1:"),
      p("Upload the unprocessed NHGRI Shipping History Detail .xlsx file."),
      fileInput("ssrs",
        label = "Upload SSRS .xlsx file",
        accept = c(".xlsx", ".xls")
      ),
      tags$style(type = 'text/css', '#ssrs_messages {color: red; border-style: none;}'),
      p(textOutput("ssrs_messages")),
      h3("Step 2:"),
      p("Upload the lay summary .xlsx file."),
      fileInput("lay",
        label = "Upload lay summary .xlsx file",
        accept = c(".xlsx", "xls")
      ),
      tags$style(type = 'text/css', '#lay_messages {color: red; border-style: none;}'),
      p(textOutput("lay_messages")),
      h3("Step 3:"),
      p("Select the date range the covered in the Shipping History Detail file."),
      dateRangeInput("date_range",
        label = "Select covering range"
      ),
      h3("Step 4:"),
      p("Process the files. Processed results will appear as a table on the right. Check the Processing Messages tab for errors or messages."),
      actionButton("run",
        label = "Process data"
      ),
      h3("Step 5:"),
      p("Preview the processed results and download .zip folder of xlsx files."),
      downloadButton("download",
        label = "Download .zip"
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Preview",
                 selectInput("dataset",
                             label = "Select dataset to preview",
                             choices = NULL
                 ),
                 selectInput("population",
                             label = "Select population to preview",
                             choices = NULL
                 ),
                 p("NOTE: 'first' and 'last' columns will be formatted as a single 'Investigator' column in the final output"),
                 DT::dataTableOutput("preview")
        ),
        tabPanel("Processing Messages",
                 verbatimTextOutput("processing_messages"))
      )
    )
  )
)

# Define server logic ----------------------------------------------------------
server <- function(input, output, session) {
    ssrs <- reactive({
        req(input$ssrs)
        df <- read_excel(input$ssrs$datapath, skip = 1) %>% janitor::clean_names()
        
        validate(
            need(identical(colnames(df), ssrs_cols), "Error: Input columns do not match expected columns.")
        )
        
        df
    })
    
    # If the output is not an error don;t render anything
    output$ssrs_messages <- renderText( {if (is_tibble(ssrs())) cat("") else ssrs()} )
    
    lay <- reactive({
        req(input$lay)
        required_cols <- c("lay_summary", "customer", "population")
        df <- read_excel(input$lay$datapath) %>% janitor::clean_names()
        
        validate(
            need(all(required_cols %in% colnames(df)), "Error: Missing a required column (either Lay Summary, Customer, or Population) in the input csv file.")
        )
        
        df
    })
    
    # If the output is not an error don;t render anything
    output$lay_messages <- renderText( {if (is_tibble(lay())) cat("") else lay()} )
    
  observeEvent(input$run, {
    report_dfs <- clean_and_split_ssrs(ssrs())
    lay_unnested <- clean_and_unnest_lay(lay())
    lay_dfs <- join_lay_onto_reports(report_dfs, lay_unnested)

    # drop population column from report dfs after joins
    report_dfs <- report_dfs %>%
      map(select, -Population) %>%
      map(arrange, Investigator)

    # split the Investigator column into first last for both sets of data frames
    report_dfs <- split_investigator(report_dfs)
    lay_dfs <- split_investigator(lay_dfs)

    # after processing, update preview choices
    updateSelectInput(session, "dataset", choices = c("Research Intent", "Lay Summary"))
    updateSelectInput(session, "population", choices = names(report_dfs))

    d <- reactive({
      if (input$dataset == "Research Intent") {
        report_dfs[[input$population]]
      } else {
        lay_dfs[[input$population]]
      }
    })

    output$preview <- DT::renderDataTable(d())
    
    output$processing_messages <- renderText({
      msgs <- check_processing(report_dfs, lay_dfs)
      paste(msgs, collapse = "\n")
      })

    output$download <- downloadHandler(
      filename = function() {
        paste0(Sys.Date(), "-NHGRI-community-reports.zip")
      },
      content = function(fname) {
        fs <- c()
        tmpdir <- tempdir()
        setwd(tmpdir)
        for (pop in names(report_dfs)) {
          path <- paste0(pop, ".xlsx")
          fs <- c(fs, path)
          create_excel(
            report_df = report_dfs[[pop]],
            lay_df = lay_dfs[[pop]],
            title = pop,
            from_date = input$date_range[[1]],
            to_date = input$date_range[[2]],
            filename = paste0(pop, ".xlsx")
          )
        }
        zip(zipfile = fname, files = fs)
      },
      contentType = "application/zip"
    )
  })
}

# Run the application ----------------------------------------------------------
shinyApp(ui = ui, server = server)
