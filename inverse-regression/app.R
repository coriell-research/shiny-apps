library(shiny)


# Inverse regression
meth2age <- function(x, model) {
  Intercept <- unname(coef(model)[1])
  Beta <- unname(coef(model)[2])
  (x - Intercept) / Beta
}

ui <- fluidPage(
  titlePanel("Inverse Regression Fit"),
  sidebarLayout(
    sidebarPanel(
      width = 2,
      numericInput("ntrain", "N Train", value = 100, min = 2, max = Inf, step = 10),
      helpText("Select the number of samples in the training set"),
      numericInput("ntest", "N Test", value = 20, min = 1, max = Inf, step = 10),
      helpText("Select the number of sample in the test set"),
      numericInput("ageLow", "Min Age", value = 0, min = 0, max = 150, step = 10),
      helpText("Specify the lowest possible age"),
      numericInput("ageHigh", "Max Age", value = 100, min = 0, max = 150, step = 10),
      helpText("Specify the highest possible age"),
      numericInput("B0", "Intercept", value = 25, min = 0, max = 100, step = 5),
      helpText("Specify the value of the intercept for the model"),
      numericInput("Slope", "Slope", value = 0.5, min = -Inf, max = Inf, step = 0.05),
      helpText("Specify the value of the slope for the model"),
      numericInput("Error", "Mean Error", value = 0, min = 0, max = Inf, step = 1),
      helpText("Specify the amount of error around the mean"),
      numericInput("SdError", "Stdev Error", value = 5, min = 0, max = Inf, step = 1),
      helpText("Specify the standard deviation of the error around the mean"),
      h2("Description:"),
      helpText("Ages for the training and testing sets are separately drawn from a uniform distribution.",
               "An initial model is made for the Training/Testing set: Methylation ~ Intercept + Slope * Age + error.",
               "Where 'error' is generated from a normal distribution with Mean = 'Error' and SD = 'Stdev Error'.", 
               "After fitting the model on the Training data, predicted ages of the Testing set are generated using the inverse regression formula: Predicted Age = Methylation - Intercept / Slope.",
               "The top left plot shows the model fit on the Training data with the predicted ages of the Testing data.",
               "The top right plot shows the relationship we wish to predict (Age ~ Metyhylation) of the Training data with the predicted ages of the Testing data projected onto it.",
               "The bottom left plot shows the Actual ages of the Testing data vs their Predicted ages from the inverse regression equation.",
               "The Bootom right plot shows the actual Testing data table."
      )
    ),
    mainPanel(
      fluidRow(
        column(6, plotOutput("scatterPlot")),
        column(6, plotOutput("invScatterPlot"))
      ),
      fluidRow(
        column(6, plotOutput("predPlot")),
        column(6, DT::dataTableOutput("df"))
      )
    )
  )
)

server <- function(input, output, session) {
  data <- reactive({
    # Generate training and testing data
    train_age <- floor(runif(input$ntrain, min = input$ageLow, max = input$ageHigh))
    train_meth <- coriell::clamp(input$B0 + (input$Slope * train_age) + rnorm(input$ntrain, input$Error, input$SdError), 0, 100)
    test_age <- floor(runif(input$ntest, min = input$ageLow, max = input$ageHigh))
    test_meth <- coriell::clamp(input$B0 + (input$Slope * test_age) + rnorm(input$ntest, input$Error, input$SdError), 0, 100)
    
    # Forward model -- Meth ~ Age
    model <- lm(train_meth ~ train_age)
    
    # Inverse regression on testing data
    predAge <- meth2age(test_meth, model)
    
    # Create data.frames from the training and testing data
    train_df <- data.frame(Age = train_age, Meth = train_meth)
    test_df <- data.frame(Age = test_age, AgeInv = round(predAge, 1), 
                          Meth = round(test_meth, 1), 
                          row.names = paste0("Test", 1:input$ntest))
    
    # Return the data.frames as reactive objects
    list(Train=train_df, Test=test_df)
  })
  
  output$scatterPlot <- renderPlot({
    Train <- data()$Train
    Test <- data()$Test
    
    plot(x=Train$Age, y=Train$Meth, xlab="Age", ylab="Methylation (%)", 
         pch = 16, col = "black", main = "Linear Regression\nMethylation ~ Age (blue)", 
         ylim = c(0, 100), xlim = c(0, input$ageHigh))
    abline(lm(Meth ~ Age, Train), col = "blue2", lty = 2, lwd = 2)
    points(x=Test$AgeInv, y=Test$Meth, col="red2", pch=16)
    legend("topleft", legend = c("Test", "Train"), col = c("red2", "black"), 
           pch=c(16, 16))
  })
  
  output$invScatterPlot <- renderPlot({
    Train <- data()$Train
    Test <- data()$Test
    
    plot(x=Train$Meth, y=Train$Age, xlab="Methylation (%)", ylab="Age",
         main="Inverse Regression\nAge ~ Methylation (blue)", xlim=c(0, 100), pch = 16,
         col = "black", ylim=c(0, input$ageHigh))
    abline(lm(Age ~ Meth, Train), col = "blue2", lty = 2, lwd = 2)
    points(x=Test$Meth, y=Test$AgeInv, col="red2", pch=16)
    legend("topleft", legend = c("Test", "Train"), col = c("red2", "black"), 
           pch=c(16, 16))
  })
  
  output$predPlot <- renderPlot({
    Test <- data()$Test
    MAE <- round(mean(abs(Test$AgeInv - Test$Age)), 1)
    MedAE <-  round(median(abs(Test$AgeInv - Test$Age)), 1)
    
    plot(x = Test$Age, y = Test$AgeInv, xlab = "Actual", ylab="Predicted", 
         main = paste0("Errors\nMAE: ", MAE, "; MedAE: ", MedAE))
    abline(a=0, b=1, lty=3, col="black", lwd=2)
  })
  
  output$df <- DT::renderDataTable({ 
    DT::datatable(
      data()$Test,
      colnames = c("Actual Age", "Age by Inv. Reg.", "Methylation"),
      caption = htmltools::tags$caption( style = 'caption-side: top; text-align: center; color:black;  font-size:200% ;','Testing Data') ) 
    })
}

shinyApp(ui, server)