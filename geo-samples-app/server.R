library(shiny)
library(data.table)
library(PCAtools)
library(plotly)
library(here)


# load data --------------------------------------------------------------------
# batch corrected count matrices
gene_mat.c <- readRDS(here("data", "batchCorrected-gene-lcpms.rds"))
re_mat.c <- readRDS(here("data", "batchCorrected-re-lcpms.rds"))

# raw count matrices
gene_mat <- readRDS(here("data", "gene-lcpms.rds"))
re_mat <- readRDS(here("data", "re-lcpms.rds"))

# data.tables of gene counts in long format
gene_dt <- readRDS(here("data", "gene-lcpm-dt.rds"))
re_dt <- readRDS(here("data", "re-lcpm-dt.rds"))

# metadata for the annotated samples
metadata <- readRDS(here("data", "sample-metadata.rds"))
setkeyv(metadata, cols = c("sample_name", "bio_project", "cell_type", "drug"))


# Define server logic required to draw a histogram -----------------------------
shinyServer(function(input, output) {
    
    M <- reactive({
        switch(input$count_matrix,
               "Batch Corrected Gene Counts" = gene_mat.c,
               "Batch Corrected RE Counts" = re_mat.c,
               "Raw Gene Counts" = gene_mat,
               "Raw RE Counts" = re_mat
        )}
    )
    
    observeEvent(input$run, {
        selected_metadata <- metadata[sample_name %in% input$samples &
                                      bio_project %in% input$projects &
                                      cell_type %in% input$cells &
                                      drug %in% input$drugs]
        # convert to dataframe
        selected_metadata <- data.frame(
            as.matrix(selected_metadata, rownames = "sample_name")
        )
        
        input_matrix <- subset(M(),
                               subset = !grepl(input$remove_regex, rownames(M())),
                               select = rownames(selected_metadata)
                               )
        stopifnot(all(colnames(input_matrix) == rownames(selected_metadata)))
        
        # perform PCA on selected input
        p <- pca(input_matrix, 
                 selected_metadata, 
                 center = input$center_data,
                 scale = input$scale_data, 
                 removeVar = input$remove_var)
        
        # make long data from PCA output for plotly
        pca.dt <- as.data.table(p$rotated, keep.rownames = "sample_name")
        setkey(pca.dt, "sample_name")
        pca.meta <- as.data.table(p$metadata, keep.rownames = "sample_name")
        setkey(pca.meta, "sample_name")
        
        # join metadata onto PC data
        pca.dt <- pca.dt[pca.meta]
        
        # create plots from PCA ------------------------------------------------
        output$pca <- renderPlotly({
            plot_ly(pca.dt,
                    x = ~get(input$pca_xvar),
                    y = ~get(input$pca_yvar), 
                    z = ~get(input$pca_zvar), 
                    color = ~get(input$pca_color),
                    colors = "Set1",
                    hovertext = paste("BioProject:", pca.dt$bio_project,
                                      "<br> Cell Type:", pca.dt$cell_type,
                                      "<br> Drug:", pca.dt$drug,
                                      "<br> Dose", pca.dt$dose,
                                      "<br> Time", pca.dt$time_hr)) %>%
                add_markers()
        })
        
        output$biplot <- renderPlot({ biplot(p, 
                                             x = input$biplot_x,
                                             y = input$biplot_y,
                                             showLoadings = input$show_loading,
                                             colby = if (input$biplot_color == "None") NULL else input$biplot_color,
                                             shape = if (input$biplot_shape == "None") NULL else input$biplot_shape,
                                             hline = 0,
                                             vline = 0,
                                             hlineType = 2,
                                             vlineType = 2,
                                             lab = if (input$biplot_labs) rownames(p$metadata) else NULL) 
            })
        output$pairsplots <- renderPlot({ pairsplot(p) })
        output$loadings <- renderPlot({ plotloadings(p) })
        output$screeplot <- renderPlot({ screeplot(p) })
    
    })
})
