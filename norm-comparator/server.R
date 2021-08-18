library(shiny)

source(here::here("helpers.R"))


server <- function(input, output, session) {
  observeEvent(input$load, {
    inFile <- input$se

    if (is.null(inFile)) {
      return(NULL)
    }

    raw <- readRDS(inFile$datapath)

    if (class(raw) != "SummarizedExperiment") {
      return(NULL)
    }

    shiny::validate(
      need("counts" %in% assayNames(raw), "'counts' must be one of the assayNames"),
      need("group" %in% colnames(colData(raw)), "'group' must be present as column in colData")
    )
    
    # update user choices after loading in raw dataset
    samples <- colnames(raw)
    metadata_cols <- colnames(colData(raw))
    updateSelectInput(session, "hist_sample", choices = samples)
    updateSelectInput(session, "reference_col", choices = c(samples, "NULL"), selected = "NULL")
    updateSelectInput(session, "rle_fill_by", choices = metadata_cols, selected = "group")
    updateSelectInput(session, "rle_fill_by2", choices = metadata_cols, selected = "group")
    updateSelectInput(session, "ma_sample1", choices = samples, selected = samples[[1]])
    updateSelectInput(session, "ma_sample2", choices = samples, selected = samples[[2]])
    updateSelectInput(session, "ma_sample12", choices = samples, selected = samples[[1]])
    updateSelectInput(session, "ma_sample22", choices = samples, selected = samples[[2]])
    updateSelectInput(session, "scatter_sample1", choices = samples, selected = samples[[1]])
    updateSelectInput(session, "scatter_sample2", choices = samples, selected = samples[[2]])
    updateSelectInput(session, "scatter_sample12", choices = samples, selected = samples[[1]])
    updateSelectInput(session, "scatter_sample22", choices = samples, selected = samples[[2]])
    updateSelectInput(session, "pca_color_by", choices = metadata_cols, selected = "group")
    updateSelectInput(session, "pca_shape_by", choices = c(metadata_cols, "NULL"), selected = "NULL")
    updateSelectInput(session, "pca_color_by2", choices = metadata_cols, selected = "group")
    updateSelectInput(session, "pca_shape_by2", choices = c(metadata_cols, "NULL"), selected = "NULL")
    updateSelectizeInput(session, "cor_samples", choices = samples, selected = samples[1:10])
    updateSelectizeInput(session, "cor_samples2", choices = samples, selected = samples[1:10])
    updateSelectizeInput(session, "hm_samples", choices = samples, selected = samples[1:10])
    updateSelectizeInput(session, "hm_samples2", choices = samples, selected = samples[1:10])

    # filter/normalize/plot on RUN
    observeEvent(input$run, {
      keep <- edgeR::filterByExpr(raw,
        group = raw$group,
        min.count = input$min_count,
        min.total.count = input$min_total_count,
        min.prop = input$min_prop
      )
      filtered <- raw[keep, ]

      # hard coded here (modify later to accept input$control_genes)
      control_genes <- c(
        "ACTB", "B2M", "C1orf43", "CHMP2A", "GAPDH", "GPI", "GUSB",
        "HPRT1", "HSP90AB1", "PPIA", "PSMB2", "PSMB4", "RAB7A", "REEP5",
        "RPL13A", "RPLP0", "SNRPD3", "TFRC", "UBC", "VCP", "VPS29"
      )

      normed <- normalize_data(filtered, input$pseudocount, input$percentile, input$reference_col, control_genes, input$K)
      output$raw_hist <- renderPlot({
        plot_distribution(raw, input$hist_sample, input$hist_bins, input$hist_pseudocount)
      })
      output$filt_hist <- renderPlot({
        plot_distribution(normed, input$hist_sample, input$hist_bins, input$hist_pseudocount)
      })
      output$rle <- renderPlot({
        plot_rle(normed, input$rle_norm_method, input$rle_fill_by, input$rle_outlier_shape, input$rle_outlier_alpha)
      })
      output$rle2 <- renderPlot({
        plot_rle(normed, input$rle_norm_method2, input$rle_fill_by2, input$rle_outlier_shape2, input$rle_outlier_alpha2)
      })
      output$ma <- renderPlot({
        plot_smear(normed, input$ma_norm_method, input$ma_sample1, input$ma_sample2, input$ma_smooth, input$ma_loess)
      })
      output$ma2 <- renderPlot({
        plot_smear(normed, input$ma_norm_method2, input$ma_sample12, input$ma_sample22, input$ma_smooth2, input$ma_loess2)
      })
      output$scatter <- renderPlot({
        plot_scatter(normed, input$scatter_norm_method, input$scatter_sample1, input$scatter_sample2, input$scatter_pt_alpha, input$scatter_log)
      })
      output$scatter2 <- renderPlot({
        plot_scatter(normed, input$scatter_norm_method2, input$scatter_sample12, input$scatter_sample22, input$scatter_pt_alpha2, input$scatter_log2)
      })
      output$pca <- renderPlot({
        plot_pca(normed, input$pca_norm_method, input$pca_scale, input$pca_center, input$pca_component1, input$pca_component2, input$pca_color_by, input$pca_shape_by)
      })
      output$pca2 <- renderPlot({
        plot_pca(normed, input$pca_norm_method2, input$pca_scale2, input$pca_center2, input$pca_component12, input$pca_component22, input$pca_color_by2, input$pca_shape_by2)
      })
      output$corplot <- renderPlot({
        plot_cor(normed, input$cor_norm_method, input$cor_samples, input$cor_viz_method)
        })
      output$corplot2 <- renderPlot({
        plot_cor(normed, input$cor_norm_method2, input$cor_samples2, input$cor_viz_method2)
      })
      output$hm <- renderPlot({
        plot_heatmap(normed, input$hm_norm_method, input$hm_samples, input$hm_features, input$hm_cd_rows, input$hm_cd_cols, input$hm_clust_method)
      })
      output$hm2 <- renderPlot({
        plot_heatmap(normed, input$hm_norm_method2, input$hm_samples2, input$hm_features2, input$hm_cd_rows2, input$hm_cd_cols2, input$hm_clust_method2)
      })
    })
  })
}
