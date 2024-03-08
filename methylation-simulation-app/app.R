library(shiny)
library(shinyWidgets)
library(data.table)
library(ggplot2)


# Function to count number of methylated loci per 4-mer (default)
patternCount <- function(x, k = 4) {
  counts <- vector("integer", length = k + 1)
  for (i in 1:(length(x) - k + 1)) {
    j <- i + k - 1
    count <- sum(x[i:j])
    counts[count + 1] <- counts[count + 1] + 1L
  }
  counts / sum(counts)
}


# UI ------------------------------------------------------------------------------------------


ui <- fluidPage(
  shinyFeedback::useShinyFeedback(),
  titlePanel("Simulate Methylation"),
  sidebarLayout(
    sidebarPanel(
      h2("Model generating function parameters:"),
      p("This tool assumes a linear change in methylation with age.
        Ages are fixed from 0-100. The parameters below are used to
        generate expected percent methylation values as a function of
        age. These percentage methylation values are then used to
        generate allele-level data by drawing random samples
        from a binomial distribution with p = expected percent
        methylation for each of the random replicates."),
      numericInput(
        "seed",
        "Random Seed",
        value = 123,
        min = 1,
        max = Inf,
        step = 1
      ),
      p("Random seed used to reproduce random sampling. Setting
        different seed values allows you to resample from the model
        generating function without changing other parameters."),
      numericInput(
        "intercept",
        "Intercept (Methylation %)",
        value = 1,
        min = -Inf,
        max = Inf,
        step = 1
      ),
      p("Percent methylation for Age=0"),
      numericInput(
        "slope",
        "Slope (Methylation %)",
        value = 0.3,
        min = -Inf,
        max = Inf,
        step = 0.1
      ),
      p("Percent change in methylation per year"),
      selectInput(
        "fun",
        "Error generating function",
        choices = c("none", "identity", "constant", "sqrt", "log1p", "arcsin"),
        selected = "none"
      ),
      p("The error generating function is used to vary the amount of
        sampling error as a function of age"),
      numericInput(
        "offset",
        "Offset",
        value = 1,
        min = 0,
        max = Inf,
        step = 0.1
      ),
      p("Offset to adjust the magnitude of the error function. Offset is
        multiplied to the result of the error function, i.e. f(age) * offset"),
      numericInput(
        "nrep",
        "N Replicates",
        value = 1,
        min = 1,
        max = Inf,
        step = 1
      ),
      p("Number of samples to simulate at each age"),
      numericInput(
        "nallele",
        "N Alleles",
        value = 100,
        min = 2,
        max = Inf
      ),
      p("Number of alleles to simulate for each sample"),
      h2("JSD Reference"),
      p("Select the proportion of methylated CpGs used as a reference for JSD calculation. 
        Proportions must sum to 1"),
      fluidRow(
        column(2, numericInput("p0", "p0", value = 0.2, min = 0, max = 1, step = 0.1)),
        column(2, numericInput("p1", "p0", value = 0.2, min = 0, max = 1, step = 0.1)),
        column(2, numericInput("p2", "p0", value = 0.2, min = 0, max = 1, step = 0.1)),
        column(2, numericInput("p3", "p0", value = 0.2, min = 0, max = 1, step = 0.1)),
        column(2, numericInput("p4", "p0", value = 0.2, min = 0, max = 1, step = 0.1)),
      ),
      h2("Generate Data"),
      p("Use the button below to generate data from the specified model.
        Be sure to regenerate data every time you make a desired change to
        the model."),
      actionButton("run", "Generate Data"),
      h2("Download data"),
      p("Data and plots are exported as a zip folder."),
      downloadButton("download")
    ),
    mainPanel(
      plotOutput("model"),
      plotOutput("simulated"),
      selectInput(
        "yaxis",
        "y-axis variable",
        choices = c("Avg. Methylation" = "AvgMeth", "JSD" = "JSD", "Residuals" = "Residual"),
        selected = "Avg. Methylation"
      )
    )
  )
)


# Server --------------------------------------------------------------------------------------


server <- function(input, output) {
  modfun <- function(slope, intercept, fun) {
    x <- 0:100
    m <- slope
    b <- intercept
    offset <- input$offset

    err <- switch(fun,
      "none" = 0,
      "identity" = x * offset,
      "constant" = 1 * offset,
      "sqrt" = sqrt(x) * offset,
      "log1p" = log1p(x) * offset,
      "arcsin" = asin(sqrt(x / 100)) * offset
    )

    eps <- rnorm(x, mean = 0, sd = err)
    y <- m * x + b + eps
    dt <- data.table(X = x, Y = y)

    return(dt)
  }

  # Results of drawing N replicates from model generating function
  simdata <- reactive({
    set.seed(input$seed)

    data <- replicate(
      input$nrep,
      modfun(input$slope, input$intercept, input$fun),
      simplify = FALSE
    )
    data <- rbindlist(data)

    # Ensures methylation values are in the range 0-100
    data[, Y := coriell::clamp(Y, min = 0, max = 100)]

    return(data)
  })

  modplot <- reactive({
    form <- paste0(
      "Methylation ~ ", input$intercept, " + ", input$slope,
      " * Age + N(0, ", input$fun, "(Age) * ", input$offset, ")"
    )

    ggplot(simdata(), aes(x = X, y = Y)) +
      geom_point(size = 1) +
      geom_smooth(se = FALSE, col = "red2", linetype = 2) +
      geom_smooth(method = "lm", se = FALSE, col = "blue2") +
      scale_y_continuous(limits = c(0, 100)) +
      scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5)) +
      labs(
        title = "Model Generating Function",
        subtitle = form,
        x = "Age",
        y = "Methylation (%)",
        caption = paste(
          "Random Seed:", input$seed,
          "; Replicates:", input$nrep,
          "; N Alleles:", input$nallele
        )
      ) +
      theme_light() +
      theme(
        plot.title = element_text(face = "bold", size = 20),
        plot.subtitle = element_text(face = "italic", size = 16)
      )
  })

  output$model <- renderPlot(modplot())

  data <- reactive({
    show_alert(
      title = "Generating Data",
      text = "Please Wait...",
      closeOnClickOutside = FALSE,
      btn_labels = NA,
    )

    # Random allele generation
    p <- simdata()[["Y"]] / 100
    n <- input$nallele

    isOne <- n == 1
    shinyFeedback::feedbackDanger("nallele", isOne, "Please select at least 2 alleles")
    req(!isOne, cancelOutput = TRUE)

    d <- t(sapply(p, \(x) rbinom(n, 1L, x)))
    cnames <- paste0("Allele", 1:n)
    colnames(d) <- cnames
    rownames(d) <- rep(0:100, input$nrep)

    # JSD calculation
    p0 <- input$p0
    p1 <- input$p1
    p2 <- input$p2
    p3 <- input$p3
    p4 <- input$p4
    total <- p0 + p1 + p2 + p3 + p4
    shinyFeedback::feedbackDanger("p4", total != 1, "Proportions must sum to 1.0")
    req(total == 1.0)
    
    k <- 4
    pattern_props <- do.call(rbind, apply(d, 1, patternCount, k = k, simplify = FALSE))
    colnames(pattern_props) <- paste0(seq(0, k, by = 1), "_Methylated")
    m <- rbind(c(p0, p1, p2, p3, p4), pattern_props)
    jsd <- sqrt(philentropy::JSD(m, test.na = FALSE))
    jsd0 <- jsd[1, 2:ncol(jsd)]  # Select every sample relative to the reference 

    d <- cbind(d, pattern_props)
    d <- as.data.table(d, keep.rownames = "Age")
    d[, Age := as.integer(Age)]
    d[, AvgMeth := Reduce(`+`, .SD) / n * 100, .SDcols = cnames]
    d[, JSD := jsd0]

    # Add a column for the expected methylation value from the linear fit
    E <- (input$slope * d$Age) + input$intercept
    d[, ExpectedMeth := E][, Residual := AvgMeth - ExpectedMeth]

    setcolorder(d, c("Age", "AvgMeth", "JSD", "ExpectedMeth", "Residual"))

    closeSweetAlert()
    return(d)
  }) |> bindEvent(input$run)


  simplot <- reactive({

    ylims <- c(0, 100)
    ylabel <- "Avg. Methylation (%)"
    if (input$yaxis == "JSD") {
      ylims <- c(0, 1)
      ylabel <- "JSD"
    } else if (input$yaxis == "Residual") {
      ylims = NULL
      ylabel = "Residual"
    }

    ggplot(data(), aes(x = Age, y = .data[[input$yaxis]])) +
      geom_point(size = 1) +
      geom_smooth(se = FALSE, col = "red2", linetype = 2) +
      geom_smooth(method = "lm", se = FALSE, col = "blue2") +
      scale_y_continuous(limits = ylims) +
      scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5)) +
      labs(
        title = "Simulated Data",
        x = "Age",
        y = ylabel
      ) +
      theme_light() +
      theme(
        plot.title = element_text(face = "bold", size = 20),
        plot.subtitle = element_text(face = "italic", size = 16)
      )
  })

  output$simulated <- renderPlot(simplot())

  output$download <- downloadHandler(
    filename = function() {
      paste0("methylation-simulator_", format(Sys.time(), "%Y-%m-%d"), ".zip")
    },
    content = function(file) {
      tmp <- tempdir()
      setwd(tmp)

      data_file <- fwrite(data(), "data.tsv", sep = "\t")
      model_plot <- ggsave("model-plot.pdf", plot = modplot(), device = "pdf", width = 9, height = 6)
      sim_plot <- ggsave("simdata-plot.pdf", plot = simplot(), device = "pdf", width = 9, height = 6)
      files <- c("data.tsv", "model-plot.pdf", "simdata-plot.pdf")

      zip(zipfile = file, files = files)
    },
    contentType = "application/zip"
  )
}

# Run the application
shinyApp(ui = ui, server = server)
