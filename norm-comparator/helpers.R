library(shiny)
library(tidyverse)
library(SummarizedExperiment)


# Function for normalizing filtered dataset
normalize_data <- function(se, pseudocount, percentile, reference_col, control_genes, K) {
  # Total count ------------------------------------------------------------------
  tc <- edgeR::cpm(assays(se)[["counts"]])
  tc_lcpm <- edgeR::cpm(assays(se)[["counts"]], log = TRUE, prior.count = pseudocount)
  
  assays(se)[["LibrarySize"]] <- tc
  assays(se)[["logLibrarySize"]] <- tc_lcpm
  
  # edgeR ------------------------------------------------------------------------
  tmm <- edgeR::calcNormFactors(se, method = "TMM", refColumn = reference_col)
  uq <- edgeR::calcNormFactors(se, method = "upperquartile", p = percentile)
  rle <- edgeR::calcNormFactors(se, method = "RLE")
  
  tmm_cpm <- edgeR::cpm(tmm)
  uq_cpm <- edgeR::cpm(uq)
  rle_cpm <- edgeR::cpm(rle)
  tmm_lcpm <- edgeR::cpm(tmm, log = TRUE, prior.count = pseudocount)
  uq_lcpm <- edgeR::cpm(uq, log = TRUE, prior.count = pseudocount)
  rle_lcpm <- edgeR::cpm(rle, log = TRUE, prior.count = pseudocount)
  
  assays(se)[["TMM"]] <- tmm_cpm
  assays(se)[["UQ"]] <- uq_cpm
  assays(se)[["RLE"]] <- rle_cpm
  assays(se)[["logTMM"]] <- tmm_lcpm
  assays(se)[["logUQ"]] <- uq_lcpm
  assays(se)[["logRLE"]] <- rle_lcpm
  
  # QSmooth ----------------------------------------------------------------------
  qs <- qsmooth::qsmooth(se, group_factor = se$group)
  qs_cpm <- edgeR::cpm(qsmooth::qsmoothData(qs))
  qs_lcpm <- edgeR::cpm(qsmooth::qsmoothData(qs), log = TRUE, prior.count = pseudocount)
  
  assays(se)[["QS"]] <- qs_cpm
  assays(se)[["logQS"]] <- qs_lcpm
  
  # RUVg -------------------------------------------------------------------------
  ruv_set <- RUVSeq::RUVg(as.matrix(assays(se)[["counts"]]), control_genes, k = K)
  ruv_cpm <- edgeR::cpm(ruv_set$normalizedCounts)
  ruv_lcpm <- edgeR::cpm(ruv_set$normalizedCounts, log = TRUE, prior.count = pseudocount)
  
  assays(se)[["RUVg"]] <- ruv_cpm
  assays(se)[["logRUVg"]] <- ruv_lcpm
  
  se
}

# Function for plotting the count distribution of a given sample ---------------
plot_distribution <- function(se, sample_name, n_bins, pseudocount) {
  assays(se)[["counts"]][, sample_name] %>%
    enframe(name = NULL, value = "x") %>%
    mutate(x = x + pseudocount) %>%
    ggplot(aes(x)) +
    geom_histogram(bins = n_bins, fill = "steelblue", color = "black") +
    theme_light() +
    scale_x_log10() +
    labs(x = "CPM")
}

# Function for plotting RLE boxplots following normalization -------------------
plot_rle <- function(se, assay_name, fill_by, outlier_shape, outlier_alpha) {
  median_expr <- apply(assays(se)[[assay_name]], 1, median)
  rle_data <- assays(se)[[assay_name]] - median_expr

  plot_df <- rle_data %>%
    as_tibble(rownames = "feature_id") %>%
    pivot_longer(-feature_id, names_to = "sample_name", values_to = "RLE") %>%
    left_join(as_tibble(colData(se), rownames = "sample_name"),
      by = "sample_name"
    )

  plot_df %>%
    ggplot(aes(x = sample_name, y = RLE, fill = .data[[fill_by]])) +
    geom_boxplot(outlier.shape = outlier_shape, outlier.alpha = outlier_alpha) +
    theme_light() +
    geom_hline(
      yintercept = median(plot_df$RLE, na.rm = TRUE),
      linetype = 2,
      color = "red"
    ) +
    labs(
      x = NULL,
      y = "RLE"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Function for plotting sample vs. sample scatter plots ------------------------
plot_scatter <- function(se, assay_name, sample1, sample2, pt_alpha, log_scale) {
  plot_df <- assays(se)[[assay_name]] %>%
    as_tibble(rownames = "feature_id") %>%
    select(feature_id, sample1, sample2)

  if (log_scale) {
    plot_df %>%
      ggplot(aes(x = .data[[sample1]], y = .data[[sample2]])) +
      geom_point(alpha = pt_alpha) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
      scale_x_log10() +
      scale_y_log10() +
      theme_light() +
      labs(title = paste(assay_name))
  } else {
    plot_df %>%
      ggplot(aes(x = .data[[sample1]], y = .data[[sample2]])) +
      geom_point(alpha = pt_alpha) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
      theme_light() +
      labs(title = paste(assay_name))
  }
}

# Function to create a smear plot of two samples -------------------------------
plot_smear <- function(se, assay_name, sample1, sample2, smooth, loess) {
  X <- assays(se)[[assay_name]][, c(sample1)]
  Y <- assays(se)[[assay_name]][, c(sample2)]
  edgeR::maPlot(X,
    Y,
    normalize = FALSE,
    smooth.scatter = smooth,
    lowess = loess,
    main = paste(assay_name, "-", sample1, "vs.", sample2),
    ylab = paste0("log2(", sample2, " / ", sample1, ")")
  )
  abline(a = 0, b = 0, lty = 2)
}

# Function for plotting the PCA of normalized data -----------------------------
plot_pca <- function(se,
                     assay_name,
                     scale_data,
                     center_data,
                     component1,
                     component2,
                     color_by,
                     shape_by) {
  pca <- prcomp(t(assays(se)[[assay_name]]), center = center_data, scale. = scale_data)
  pca_df <- as_tibble(pca$x, rownames = "sample_name") %>%
    left_join(as_tibble(colData(se), rownames = "sample_name"), by = "sample_name")

  if (shape_by == "NULL") {
    pca_df %>%
      ggplot(aes(x = .data[[component1]], y = .data[[component2]], color = .data[[color_by]])) +
      geom_point(size = 3) +
      geom_hline(yintercept = 0, linetype = 2) +
      geom_vline(xintercept = 0, linetype = 2) +
      theme_light()
  } else {
    pca_df %>%
      ggplot(aes(x = .data[[component1]], y = .data[[component2]], color = .data[[color_by]], shape = .data[[shape_by]])) +
      geom_point(size = 3) +
      geom_hline(yintercept = 0, linetype = 2) +
      geom_vline(xintercept = 0, linetype = 2) +
      theme_light()
  }
}

# Function for creating correlation plots --------------------------------------
plot_cor <- function(se, assay_name, cor_samples, viz_method) {
  cor_mat <- cor(assays(se)[[assay_name]][, cor_samples])
  corrplot::corrplot(
    cor_mat,
    method = viz_method,
    type = "lower", 
    col = viridis::viridis(100), 
    outline = TRUE,
    tl.col = "black",
    tl.srt = 45,
    diag = TRUE)
}

# Function for plotting heatmaps -----------------------------------------------
plot_heatmap <- function(se, assay_name, samples, n_features, cd_rows, cd_cols, cd_method) {
  mat <- assays(se)[[assay_name]]
  row_vars <- matrixStats::rowVars(mat)
  var_order <- order(row_vars, decreasing = TRUE)
  var_mat <- mat[head(var_order, n_features), samples]
  var_mat <- as.matrix(var_mat)
  var_mat <- na.omit(var_mat)
  pheatmap::pheatmap(var_mat, 
                     scale = "row",
                     color = colorRampPalette(c("dodgerblue3", "grey99", "firebrick3"))(50),
                     show_rownames = FALSE,
                     treeheight_row = 0,
                     clustering_distance_rows = cd_rows,
                     clustering_distance_cols = cd_cols,
                     clustering_method = cd_method)
}
