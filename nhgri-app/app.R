library(shiny)
library(shinythemes)
source("helpers.R")

# required columns for reading in data
ssrs_required <- c("Collection_Type_id", "Ref", "Exp", "Product", "From_Product", 
                   "Quantity_in_Container", "Quantity", "Quantity2", "Order_Date", 
                   "Order_Id", "Ship_Date", "Order_Type", "DiagDesc", "RIntentType", 
                   "RIntent", "Name", "Customer_Id", "Institution", "Institution_Type", 
                   "Country", "Order_Remark", "Invoice_Item_Price", "Queue_Price_11", 
                   "Queue_Price_21", "Queue_Price_31", "Invoice_Amt", "Price_Adjustment", 
                   "Ship_Charge", "Special_Handling_Charge", "Prepayment_Amt", "Invoice_Paid_Amt")

lay_required <- c("Lay_Summary1", "Customer", "Country", "Institution_Name", 
                  "Population", "Collection", "Addtl_Info")


# Define UI for application ----------------------------------------------------
ui <- fluidPage(
  theme = shinytheme("cosmo"),
  titlePanel("Process NHGRI Orders"),
  sidebarLayout(
    sidebarPanel(
      h3("Step 1:"),
      p("Upload the unprocessed NHGRI Shipping History Detail .csv file."),
      fileInput("ssrs",
        label = "Upload 'Shipping History' .csv file",
        accept = c(".csv")
      ),
      tags$style(type = 'text/css', '#ssrs_messages {color: red; border-style: none;}'),
      p(textOutput("ssrs_messages")),
      h3("Step 2:"),
      p("Upload the 'Lay Summary' .csv file."),
      fileInput("lay",
        label = "Upload 'Lay Summary' .csv file",
        accept = c(".csv")
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
        dt <- readInput(input$ssrs$datapath, skip_rows = 0)
        
        validate(
            need(identical(colnames(dt), ssrs_required), "Error: Input columns do not match expected columns.")
        )
        
        dt
    })
    
    # If the output is not an error don't render anything
    output$ssrs_messages <- renderText( {if (is.data.table(ssrs())) cat("") else ssrs()} )
    
    lay <- reactive({
        req(input$lay)
        dt <- readInput(input$lay$datapath, skip_rows = 3)
        
        validate(
            need(identical(colnames(dt), lay_required), "Error: Missing a required column (either Lay Summary, Customer, or Population) in the input csv file.")
        )
        
        dt
    })
    
    # If the output is not an error don't render anything
    output$lay_messages <- renderText( {if (is.data.table(lay())) cat("") else lay()} )
    
  observeEvent(input$run, {
    report_dfs <- cleanGroupSplitSSRS(ssrs())
    lay_unnested <- cleanAndUnestLay(lay())
    lay_dfs <- joinLayOntoSSRS(report_dfs, lay_unnested)

    # drop population column from report dfs after joins and reorder
    report_dfs <- lapply(report_dfs, function(dt) {dt[, Population := NULL][order(Investigator)]})

    # split the Investigator column into first last for both sets of data frames
    report_dfs <- lapply(report_dfs, splitInvestigator)
    lay_dfs <- lapply(lay_dfs, splitInvestigator)

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
      msgs <- checkProcessing(report_dfs, lay_dfs)
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
