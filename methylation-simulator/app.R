library(shiny)
library(shinyMatrix)
library(philentropy)


# Helper functions ---
simulateReads <- function(depth, n, probs, sparsity = 0) {
  simulateRead <- function(...) {
    stopifnot("length(probs) must be equal to n" = length(probs) == n)
    read <- vector("numeric", length = n)
    for (i in 1:n) {
      u <- runif(1)
      if (u <= probs[i]) {
        read[i] <- 1L
      } else {
        read[i] <- 0L
      }
    }
    return(read)
  }
  
  # Simulate many reads and form into matrix
  reads <- replicate(depth, simulateRead(n, probs), simplify = FALSE)
  m <- do.call(rbind, reads)
  m[sample(seq_along(m), sparsity * length(m), replace = FALSE)] <- NA
  
  return(m)
}


methylationLevel <- function(m) {
  colMeans(m, na.rm = TRUE)
}


shannonEntropy <- function(x) {
  t1 <- x * ifelse(is.infinite(log2(x)), 0, log2(x))
  t2 <- (1 - x) * ifelse(is.infinite(log2(1 - x)), 0, log2(1-x))
  -1 * (t1 + t2)
}


normalizedEntropy <- function(x) {
  N <- length(x)
  (1 / (N + 1)) * -1 * sum(x * ifelse(is.infinite(log2(x)), 0, log2(x)), na.rm = TRUE)
}


patternProbs <- function(m) {
  counts <- rowSums(m, na.rm = TRUE)       # Number of methylated bases per read
  f <- factor(counts, levels = 0:ncol(m))  # Number of CpGs
  patterns <- tapply(counts, f, length)    # Sum counts per pattern
  patterns[is.na(patterns)] <- 0L          # Replace NAs for missing patterns
  probs <- patterns / nrow(m)              # Convert to probability
  names(probs) <- paste0(names(probs), "/", ncol(m))
  
  return(probs)
}

makeReference <- function(x) {
  x1 <- rep(FALSE, x)
  x2 <- rep(TRUE, x)
  reads <- vector("list", length = x - 1)
  for (i in 1:floor(x / 2)) {
    r1 <- rep(c(TRUE, FALSE), c(i, x-i))
    r2 <- !r1
    reads[[i]] <- r1
    reads[[x-i]] <- r2
  }
  xx <- do.call(rbind, reads)
  m <- do.call(rbind, list(x1, x2, xx))
  dimnames(m) <- list(paste0("read.", 1:(x+1)), paste0("CpG.", 1:x))
  m * 1L
}

# Set up default input values 
row1 <- c(100, 0.5, 0.5, 0.5, 0.5, 0.5)
rown <- matrix(data = 0L, nrow = 4, ncol = 6)
m <- rbind(row1, rown)
dimnames(m) <- NULL

# app ---
ui <- fluidPage(
  titlePanel("Methylation Comparison Simulator"),
  sidebarLayout(
    sidebarPanel(
      h3("Read Counts & CpG Probabilities"),
      p("Enter into the first column of the matrix the desired read depth for the 
        pattern. Enter into the subsequent columns the probability of being methylated
        for each individual CpG. For example, the default has been set to generate 
        100 reads from and amplicon with 5 CpGs each with a 50% probability of 
        being methytlated.", 
        br(), br(), 
        "New rows and columns can be added to the matrix dynamically to account 
        for more patterns and more CpGs.",
        br(), br(),
        "The sparsity input allows the user to replace a given percentage of cpg 
        values with NAs to simulate missing data. For example, if the total read 
        depth is 100 and there are 5 CpGs, if sparsity is set to 0.1 then 50 
        values in the matrix will be randomly replaced with NAs.",
        h4("Input Sample Parameters:")
        ),
      matrixInput(
        "sample",
        value = m,
        class = "numeric",
        rows = list(extend = TRUE),
        cols = list(extend = TRUE)
      ),
      numericInput(
        "sparsity",
        label = "Sparsity",
        value = 0,
        min = 0,
        max = 1,
        step = 0.1
      ),
      selectInput(
        "ref_type",
        label = "Reference Type",
        choices = c("Uniform", "Unmethylated", "Methylated", "FiftyPercent", "Random"),
        multiple = FALSE,
        selected = "Unmethylated"
        ),
      p(
        h4("Notes on references:"),
        br(),
        "Uniform reference creates a reference with the following properties. The proportion of reads with 0..n methylated reads is equal to 0.5. 
        The mean methylation at each CpG is 0.5. The entropy at each CpG is 1. Constructing a uniform reference should be done with
        only an odd number of CpGs. The depth of the uniform reference is arbitrary and not used in any downstream calculation.",
        br(), br(),
        "Methylated reference is a constructed with every CpG having a 100% probability of being methylated",
        br(), br(),
        "Unmethylated reference is constructed with every CpG having a 0% probabilty of being methylated",
        br(), br(),
        "FiftyPercent reference is constructed with every CpG having a 50% probability of being methylated",
        br(), br(),
        "Random reference is constructed with every CpG having a random probability of being methylated",
        br(), br(),
        h4("Description of summary stats:"),
        br(),
        "JSD: Jensen-shannon distance performed on the probability distributions from the number of methylated bases per cpg (i.e. distributions represented by the barplots)",
        br(), br(), 
        "NME: Normalized Methylation Entropy of the sample minus the reference.",
        br(), br(),
        "MML: Mean methylation level. The grand mean of the mean methylation at each CpG. 
        Delta MML is sample - reference",
        br(), br(),
        "H: Shannon entropy. Shannon's entropy at each base. Delta H is mean(H) sample - mean(H) reference",
        )
    ),

    # Show a plot of the generated distribution
    mainPanel(
      fluidRow(
        column(6,
          plotOutput("ref_hist")
          ),
        column(6,
          plotOutput("ref_hm")
          )
      ),
      fluidRow(
        column(6, 
          plotOutput("sample_hist")
          ),
        column(6,
         plotOutput("sample_hm")
         ),
        verbatimTextOutput("summary")
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # Get the input matrix
  data <- reactive({ 
    d <- input$sample
    depths <- d[, 1]               # Vector of depths
    prob_mat <- d[, 2:ncol(d)]     # Matrix of CpG probabilities
    cpgs <- ncol(d) - 1            # Total CpGs
    depth <- sum(d, na.rm = TRUE)  # Total depth
    
    list(depths = depths, prob_mat = prob_mat, cpgs = cpgs, depth = depth)
    })
  
  # Create the sample matrix from the user input
  SAMPLE <- reactive({
    depths <- data()$depths
    n <- data()$cpgs
    probs <- asplit(data()$prob_mat, 1)
    
    # Create list of matrices from user input
    l <- mapply(
      simulateReads, 
      depth = depths, 
      n = n, 
      probs = probs, 
      sparsity = input$sparsity, 
      SIMPLIFY = FALSE
      )
    
    # Remove the all 0 depths from list
    l <- l[sapply(l, function(x) length(x) > 0)]
    
    # Bind into single matrix
    m <- do.call(rbind, l)
    dimnames(m) <- list(paste0("read", 1:nrow(m)), paste0("CpG.", 1:ncol(m)))
    m
  })
  
  # Create the reference based on user input
  REF <- reactive({
    depth <- data()$depth
    n <- data()$cpgs
    
    switch(
      input$ref_type,
      Uniform = makeReference(n),
      FiftyPercent = simulateReads(depth, n, rep(0.5, n)),
      Unmethylated = simulateReads(depth, n, rep(0, n)),
      Methylated = simulateReads(depth, n, rep(1, n)),
      Random = simulateReads(depth, n, runif(n)),
      )
  })
  
  # Calculate statistics ---
  jsd <- reactive({
    P <- patternProbs(REF())
    Q <- patternProbs(SAMPLE())
    x <- rbind(P, Q)
    round(JSD(x, test.na = TRUE, unit = "log2"), 3)
  })
  
  nme <- reactive({
    nme_ref <- normalizedEntropy(methylationLevel(REF()))
    nme_sample <- normalizedEntropy(methylationLevel(SAMPLE()))
    round(nme_sample - nme_ref, 3)
  })
  
  mml <- reactive({
    mml_ref <- mean(methylationLevel(REF()))
    mml_sample <- mean(methylationLevel(SAMPLE()))
    round(mml_sample - mml_ref, 3)
  })
  
  H <- reactive({
    h_ref <- mean(shannonEntropy(methylationLevel(REF())))
    h_sample <- mean(shannonEntropy(methylationLevel(SAMPLE())))
    round(h_sample - h_ref, 3)
  })
  
  # Output plots ---
  output$ref_hist <- renderPlot({
    barplot(
      patternProbs(REF()), 
      ylim = c(0, 1),
      main = "Reference",
      xlab = "Number Methylated CpGs",
      ylab = "Proportion Reads") 
  })
  
  output$ref_hm <- renderPlot({
    heatmap(
      REF(), 
      Rowv = NA, 
      Colv = NA, 
      scale = "none", 
      col = topo.colors(4),
      main = "Reference Reads",
      ylab = "Read",
      labCol = paste0("CpG.", 1:(data()$cpgs))
    )
  })
  
  output$sample_hist <- renderPlot({ 
    barplot(
      patternProbs(SAMPLE()), 
      ylim = c(0, 1),
      main = paste0("Sample\nJSD: ", jsd(), "; dNME:", nme(), "; dMML: ", mml(), "; dH: ", H()),
      xlab = "Number Methylated CpGs",
      ylab = "Proportion Reads") 
    })
  
  output$sample_hm <- renderPlot({
    heatmap(
      SAMPLE(), 
      Rowv = NA, 
      Colv = NA, 
      scale = "none", 
      col = topo.colors(4),
      main = "Simulated Reads",
      ylab = "Read",
      labCol = paste0("CpG.", 1:(data()$cpgs))
    )
    
    output$summary <- renderText({
      paste("SUMMARY:\n",
        "JSD:", jsd(), "\n",
        "Delta Mean Methylation:", mml(), "\n",
        "Delta Normalized Methylation Entropy:", nme(), "\n",
        "Delta Mean Shannon Entropy:", H(), "\n\n",
        "Sample Mean Methylation per CpG:\n", paste(round(methylationLevel(SAMPLE()), 3), collapse = ","),
        "\n\n Sample Entropy of each CpG:\n", paste(round(shannonEntropy(methylationLevel(SAMPLE())), 3), collapse = ","),
        "\n\n Reference Mean Methylation per CpG:\n", paste(round(methylationLevel(REF()), 3), collapse = ","),
        "\n\n Reference Entropy of each CpG:\n", paste(round(shannonEntropy(methylationLevel(REF())), 3), collapse = ",")
        )
    })
  })
}

# Run the application
shinyApp(ui = ui, server = server)
