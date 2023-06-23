# Shiny app to simulate a simple model of clonal expansion
# Gennaro Calendo
# 2023-06-03
#
# -----------------------------------------------------------------------------
library(shiny)
library(data.table)
library(ggplot2)
library(paletteer)
library(plotly)


# Define UI for application
ui <- fluidPage(
  titlePanel("Clonal Expansion Simulation"),
  sidebarPanel(
    numericInput("N", label="Generations", value=100, min=1, max=Inf),
    helpText("How many generations to simulate"),
    numericInput("phenotypes", label="Phenotypes", min=1, max=26, value=5),
    helpText("How many different cell types to simulate."),
    fluidRow(
      column(4, uiOutput("cells")),
      column(4, uiOutput("meth")),
      column(4, uiOutput("fitness"))
      ),
    helpText("'Cells' sets the initial number of cells for the given phenotype.",
             "'Methylation' sets the average methylation for the given phenotype.",
             "'Fitness' sets the relative probability for each cell of a given class to be selected at each generation.")
  ),
  mainPanel(
    plotlyOutput("evo"),
    plotlyOutput("avgmeth"),
    textOutput("debug")
    ),
  )


# Define server logic 
server <- function(input, output) {
  N_PHENOTYPES <- reactive({ input$phenotypes })
  PHENOTYPES <- reactive({ LETTERS[1:N_PHENOTYPES()] })
  
  # Generate ui inputs depending on user input ---
  output$cells <- renderUI({
    tagList(lapply(PHENOTYPES(), function(x) { 
      numericInput(paste0("cell", x), 
                   label=paste("Cells:", x), 
                   min=0, 
                   max=100, 
                   value=100,
                   step=1) })
      )
    })
  
  output$meth <- renderUI({
    tagList(lapply(PHENOTYPES(), function(x) { 
      numericInput(paste0("meth", x), 
                   label=paste("Methylation:", x), 
                   min=0, 
                   max=1, 
                   value=round(runif(1), 2),
                   step=0.1) })
    )
  })
  
  output$fitness <- renderUI({
    tagList(lapply(PHENOTYPES(), function(x) { 
      numericInput(paste0("fit", x), 
                   label=paste("Fitness:", x), 
                   min=0, 
                   max=Inf, 
                   value=1,
                   step=0.1) })
      )
  })
  
  output$debug <- renderText({
    N <- input$N
    phenotypes <- PHENOTYPES()
    cells <- sapply(phenotypes, function(x) input[[paste0("cell", x)]])
    fits <- sapply(phenotypes, function(x) input[[paste0("fit", x)]])
    
    paste("Simulations:", N, ";",
          "Phenotypes:", paste(phenotypes, collapse=","), ";",
          "Cells:", paste(cells, collapse=","), ";",
          "Fitness:", paste(fits, collapse=","))
  })
  
  # Simulation ---
  
  # Create the simulated data
  simdata <- reactive({
    N <- input$N
    phenotypes <- PHENOTYPES()
    cells <- sapply(phenotypes, function(x) input[[paste0("cell", x)]])
    fitness <- sapply(phenotypes, function(x) input[[paste0("fit", x)]])
    names(fitness) <- phenotypes

    # Generate the initial population pool
    population <- rep(phenotypes, times = cells)
    probs <- fitness[population]

    # Accumulate data for each generation
    generations <- vector("list", N)
    names(generations) <- 1:N
    for (i in 1:N) {
      generations[[i]] <- data.table(table(population, dnn = "Phenotype"))
      population <- sample(population, replace = TRUE, prob = probs)
      probs <- fitness[population]
      probs[is.na(probs)] <- 0L
    }

    # Collapse all into a single data.table
    freq <- rbindlist(generations, idcol = "Generation")
    freq[, Generation := as.numeric(Generation)]

    return(freq)
  })

  # Plot the population trajectory across generations
  output$evo <- renderPlotly({
    ggplotly(
      ggplot(simdata(), aes(x = Generation, y = N, color = Phenotype)) +
      geom_line(linewidth = 1) +
      scale_color_paletteer_d("ggsci::default_igv") +
      labs(title = paste("Population Structure After", input$N, "Generations"),
           x = "Generation",
           y = "Number of Cells") +
      theme_minimal() +
      theme(legend.position = "bottom",
            plot.title = element_text(face = "bold", size = 18)
        )
    )
  })

  # Join the methylation data onto the simulation data and summarize the
  #   percent methyation after N generations
  methdata <- reactive({
    meth <- sapply(PHENOTYPES(), function(x) input[[paste0("meth", x)]])
    mm <- data.table(
      Phenotype = PHENOTYPES(),
      avg_methylation = meth
    )

    # Calculate the overall methylation at each generation
    by_generation <- simdata()[mm, on = "Phenotype"][,
                              .(mean_methylation = weighted.mean(avg_methylation, w = N),
                                min_methylation = min(avg_methylation),
                                max_methylation = max(avg_methylation)),
                              by = Generation]

    return(by_generation)
  })

  # Plot the average percent methylation of the population over time
  output$avgmeth <- renderPlotly({
    ggplotly(
      ggplot(methdata(), aes(x = Generation, y = mean_methylation)) +
      geom_line(linewidth = 2, color = "midnightblue") +
      geom_line(aes(y = min_methylation), color = "red2", linewidth = 0.5, linetype = 3) +
      geom_line(aes(y = max_methylation), color = "red2", linewidth = 0.5, linetype = 3) +
      labs(title = paste("Population Average Methylation Across", input$N, "Generations"),
           x = "Generation",
           y = "Methylation"
     ) +
     scale_y_continuous(limits = c(0, 1)) +
     theme_minimal() +
     theme(plot.title = element_text(face = "bold", size = 18))
    )
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
