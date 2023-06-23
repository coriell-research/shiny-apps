library(shiny)
library(coriell)
library(ggplot2)
library(data.table)
library(genefilter)


# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("T-test power simulation"), 
  sidebarLayout(
    sidebarPanel(
      numericInput("n1", "N Group A", value = 10, min = 1, max = Inf),
      numericInput("n2", "N Group B", value = 10, min = 1, max = Inf),
      numericInput("u1", "Mean Group A", value = 0.5, min = 0, max = 1, step = 0.1),
      numericInput("sd1", "Stdev Group A", value = 0.1, min = 0, max = Inf, step = 0.1),
      numericInput("u2", "Mean Group B", value = 0.5, min = 0, max = 1, step = 0.1),
      numericInput("sd2", "Stdev Group B", value = 0.1, min = 0, max = Inf, step = 0.1),
      numericInput("sims", "Simulations", value = 1e3, min = 1, max = Inf, step = 100),
      numericInput("alpha", "Alpha level", value = 0.05, min = 0, max = 1, step = 0.01)
      ),
    mainPanel(
      tabsetPanel(
        tabPanel("Histogram",
          plotOutput("hist")
          ),
        tabPanel("t-tests",
          DT::dataTableOutput("results")
         ),
        tabPanel("Boxplots",
          sliderInput("simulation", 
            "Simulation #", value = 1, min = 1, max = Inf, ticks = TRUE, animate = TRUE),
          plotOutput("boxplot")
          ),
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  # Update the simulations available to the boxplot function
  observeEvent(input$sims, {
    updateSliderInput(inputId = "simulation", max = input$sims)
  }) 
  
  # Generate simulated data and perform t-test on each simulation
  data <- reactive({
    # Generate simulated data
    a <- replicate(n=input$sims, expr = clamp(rnorm(input$n1, mean=input$u1, sd=input$sd1), 0, 1))
    b <- replicate(n=input$sims, expr = clamp(rnorm(input$n2, mean=input$u2, sd=input$sd2), 0, 1))
    a <- t(as.matrix(a))
    b <- t(as.matrix(b))
    colnames(a) <- paste0("A", 1:input$n1)
    colnames(b) <- paste0("B", 1:input$n2)
    
    # Bind simulated data into a single matrix for testing
    m <- cbind(a, b)
    
    # Shape data into ggplot friendly format for plotting
    dt <- as.data.table(data.frame(m), keep.rownames = "Simulation")
    dt.m <- melt(dt, id.vars = "Simulation", variable.name = "Sample", value.name = "Methylation")
    dt.m[, Group := fifelse(Sample %like% "^A", "A", "B")]
    
    # Calculate the test results for all simulations
    Group <- factor(rep(c("A", "B"), times = c(input$n1, input$n2)))
    results <- rowttests(m, Group)
    results <- as.data.table(results, keep.rownames = "simulation")
    results[, `:=`(dm = round(dm, 1), 
                   statistic = round(statistic, 1), 
                   p.value = round(p.value, 2),
                   sig = fifelse(p.value < input$alpha, TRUE, FALSE))]
    
    # Return test results and simulated data
    list(results=results, simdata=dt.m)
  })
  
  # Render the results from the row t-tests
  output$results <- DT::renderDataTable({
    data()[["results"]]
  })
  
  # Render the p-value distribution plot
  output$hist <- renderPlot({
    df <- data()[["results"]]
    prop <- table(df$sig) / input$sims
    
    ggplot(df, aes(x = p.value)) + 
      geom_histogram() +
      geom_vline(xintercept = input$alpha, color = "red2", linetype = 3, linewidth = 2) +
      labs(
        title = "P.Value Distribution Across Simulations",
        subtitle = paste("Proportion Simulations Where P <", input$alpha, "(i.e. Power) :", prop[["TRUE"]]),
        x = "P.Value",
        y = "Count") +
      theme_coriell()
  })
  
  # Render the simulation boxplots
  output$boxplot <- renderPlot({
    df <- data()[["simdata"]]
    df <- df[Simulation == input$simulation]

    ggplot(df, aes(x = Group, y = Methylation, color = Group)) +
      geom_jitter(width = 0.1, size = 3, alpha = 0.8) +
      stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                   geom = "crossbar", width = 0.4, color = "black") +
      labs(title = paste("Simulated Experiment: ", input$simulation),
           x = "Group",
           y = "Methylation") +
      theme_coriell() +
      scale_color_brewer(palette = "Set1") +
      scale_y_continuous(limits = c(0, 1))
      
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
