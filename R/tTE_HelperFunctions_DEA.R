TargetedApproach <- function(DESeq2_run, verbose = TRUE, fc_threshold, pval_threshold, sig_axis, target_factor) {

  if (verbose) message("- Executing a targeted analysis.\n", "- In this analysis no heatmap or PCA plot will be produced.\n", "- A volcano plot will be generated instead.")

  SummarizedExperiment::colData(DESeq2_run)[[target_factor]] <- factor(SummarizedExperiment::colData(DESeq2_run)[[target_factor]])
  levels <- levels(SummarizedExperiment::colData(DESeq2_run)[[target_factor]])

  res_names <- DESeq2::resultsNames(DESeq2_run) # Get all results names
  combs <- utils::combn(levels, 2, simplify = FALSE) # Generate all pairwise combinations

  results.list <- list()

  if (verbose) {
    message("    Levels found in target factor: ", paste(levels, collapse = ", "))
    message("    Number of pairwise combinations: ", length(combs))
  }

  for (cmb in combs) {
    # Construct possible result names
    contrast_name_forward  <- paste0(target_factor, "_", cmb[2], "_vs_", cmb[1])
    contrast_name_backward <- paste0(target_factor, "_", cmb[1], "_vs_", cmb[2])

    # Check which one exists in resultsNames()
    if (contrast_name_forward %in% res_names) {
      contrast <- c(target_factor, cmb[2], cmb[1])
      contrast_name <- paste(cmb[2], "vs", cmb[1], sep = "_")
    } else if (contrast_name_backward %in% res_names) {
      contrast <- c(target_factor, cmb[1], cmb[2])
      contrast_name <- paste(cmb[1], "vs", cmb[2], sep = "_")
    } else {
      warning("Contrast ", cmb[1], " vs ", cmb[2], " not found in resultsNames(); skipping.")
      next
    }

    if (verbose) message("Processing contrast: ", contrast_name)

    target_DESeq2_results <- as.data.frame(DESeq2::results(DESeq2_run, contrast = contrast))

    target_results_tibble <- tibble::tibble(
      feature = rownames(target_DESeq2_results),
      log2_FC = target_DESeq2_results$log2FoldChange,
      neg_log10_adjusted_p_value = -log10(target_DESeq2_results$padj)
    ) %>%
      dplyr::filter(!is.na(.data$log2_FC) & !is.na(.data$neg_log10_adjusted_p_value))

    results.list[[paste0("results_", contrast_name)]] <- target_results_tibble

    if (verbose) message("- Generating volcano plot.")
    volcano_plot <- GenerateVolcanoPlot( data = target_results_tibble, fc_threshold = fc_threshold, pval_threshold = pval_threshold, sig_axis = sig_axis) +
      ggplot2::ggtitle(contrast_name)

    results.list[[paste0("volcano_", contrast_name)]] <- volcano_plot
  }

  return(results.list)
}

ExploratoryApproach <- function(metadata, vst, heatmap = TRUE,
                                dim.reduct = c("PCA", "UMAP", "tSNE"),
                                numPC = 5, color_factor, shape_factor = NULL,
                                label_factor = NULL, color_palette = NULL, verbose = TRUE) {

  results.list <- list() # Store all outputs

  # 1. Heatmap
  if (heatmap) {
    if (verbose) message("- Generating heatmap with Euclidean distances.")
    heatmap_plot <- ProduceHeatmapDiffExp(data = vst)
    heatmap_plot <- grDevices::recordPlot()
    results.list[["heatmap"]] <- heatmap_plot
  }

  # 2. Dimensionality reduction plots
  if (!missing(dim.reduct)) {
    dim.reduct <- match.arg(dim.reduct)

    if (dim.reduct == "PCA" & numPC < 2) {
      message("- Invalid `numPC` given, using default `numPC = 2`.")
      numPC <- 2
    }

    if (is.null(color_factor) || !(color_factor %in% colnames(metadata))) stop("Please specify a `color_factor` present in the `metadata`.")
    if (!is.null(shape_factor) && !(shape_factor %in% colnames(metadata))) stop("Please specify a `shape_factor` present in the `metadata`.")
    if (!is.null(label_factor) && !(label_factor %in% colnames(metadata))) stop("Please specify a `label_factor` present in the `metadata`.")

    if (verbose) message("- Running dimensionality reduction: ", dim.reduct)

    dim_plots <- RunDimReduct(vst = vst, metadata = metadata, color_factor = color_factor,
                              shape_factor = shape_factor,  label_factor = label_factor,
                              color_palette = color_palette, dim.reduct = dim.reduct, numPC = numPC)

    results.list[["dim_reduction"]] <- dim_plots
  }

  return(results.list)
}

MakeScatterPlot <- function(df, x_col, y_col, color_factor, shape_factor = NULL, label_factor = NULL, title = NULL, color_palette = NULL) {

  # Helper function to generate scatter plots
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[color_factor]]))
  if (!is.null(shape_factor)) p <- p + ggplot2::aes(shape = .data[[shape_factor]])
  if (!is.null(label_factor)) p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = .data[[label_factor]]),
                                                                size = 2, show.legend = FALSE)
  p <- p + ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::labs(title = title, color = color_factor)
  if (!is.null(color_palette)) p <- p + ggplot2::scale_color_manual(values = color_palette)

  return(p)
}

RunDimReduct <- function(vst, metadata, color_factor, shape_factor = NULL, label_factor = NULL,
                         dim.reduct = c("PCA", "UMAP", "tSNE"), numPC = 5, color_palette = NULL) {

  # Main dimensionality reduction function
  dim.reduct <- match.arg(dim.reduct)
  plots <- list()

  # 1. Compute coordinates
  if (dim.reduct == "PCA") {
    # pca_res <- stats::prcomp(t(vst), scale. = TRUE, center = TRUE)
    # coords <- pca_res$x[, 1:numPC]
    # colnames(coords) <- paste0("PC", 1:numPC)
    # variance <- 100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)) # Variance explained
    # plots[["ElbowPlot"]] <- ProduceElbowPlot(data = pca_res, variance = variance) # Elbow plot

    # Transpose vst to samples x features
    mat <- t(vst)

    # Remove zero or near-zero variance columns (features)
    low_var_cols <- which(matrixStats::colVars(mat, useNames = TRUE) < 1e-5)
    if (length(low_var_cols) > 0) {
      mat <- mat[, -low_var_cols, drop = FALSE]
      message("- Removed ", length(low_var_cols), " low-variance features before PCA.")
    }

    pca_res <- stats::prcomp(mat, scale. = TRUE, center = TRUE)
    coords <- pca_res$x[, 1:min(numPC, ncol(pca_res$x))]
    colnames(coords) <- paste0("PC", 1:ncol(coords))

    variance <- 100 * (pca_res$sdev^2 / sum(pca_res$sdev^2))
    plots[["ElbowPlot"]] <- ProduceElbowPlot(data = pca_res, variance = variance)

  } else if (dim.reduct == "UMAP") {
    coords <- uwot::umap(t(vst))
    colnames(coords) <- c("Dim1", "Dim2")
  } else if (dim.reduct == "tSNE") {
    coords <- Rtsne::Rtsne(t(vst))$Y
    colnames(coords) <- c("Dim1", "Dim2")
  }

  plot_df <- cbind(as.data.frame(coords), metadata)

  # 2. Generate plots
  if (dim.reduct == "PCA") {
    for (i in 2:numPC) {
      pc_y <- colnames(coords)[i]
      title <- sprintf("PC1 vs %s (%.1f%% vs %.1f%%)", pc_y, variance[1], variance[i])
      plots[[paste0("PC1_vs_", pc_y)]] <- MakeScatterPlot(plot_df, "PC1", pc_y, color_factor, shape_factor, label_factor, title, color_palette)
    }
  } else {
    plots[[dim.reduct]] <- MakeScatterPlot(plot_df, colnames(coords)[1], colnames(coords)[2], color_factor, shape_factor, label_factor, dim.reduct, color_palette)
  }

  return(plots)
}

ComputeDESeq2 <- function(data, metadata, corr_factor, reduce){

  ###
  # CALL: ComputeSizeCorrection() and DESeq2runner()
  # DESCRIPTION: This function can be executed as an exploratory or targeted approach.
  # The targeted approach requires the parameter `targets` to be a table with two columns:
  # (i) specify the conditions to analyze and (ii) specify the labels to the selected conditions.
  # Returns the DESeq2 object.
  ###

  # R has a limit in the length of an integer that can be analyzed.
  if (any(data > .Machine$integer.max)){
    data <- round(data / reduce) # If the limit is overcome the `reduce` parameter will be used to divide all the values by the same factor.
    message("- Some values of the input matrix exceed the maximum value accepted by R.\n",
            paste("- The matrix has been divided by the factor, `reduce` = ", reduce,"."))
  }

  colData <- tibble::as_tibble(metadata) # Force the metadata to be interpreted as a tibble
  formula <- stats::as.formula(paste("~", corr_factor)) # Define the formula to use by DESeq2 to correct the `data` (`corr_factor`)

  # Run DESeq2: DESeqDataSetFromMatrix() and DESeq()
  if (ncol(data) != nrow(colData)) stop("Inconsistent features between data and metadata.")
  DESeq2_run <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(countData = data, colData = colData, design = formula))
  DESeq2_run <- DESeq2::estimateSizeFactors(DESeq2_run, type = "poscounts")
  DESeq2_run <- suppressMessages(DESeq2::DESeq(DESeq2_run, quiet = TRUE))

  return(DESeq2_run)
}

GenerateVolcanoPlot <- function(data, fc_threshold, pval_threshold, sig_axis){

  ###
  # CALL: DESeq2runner()
  # DESCRIPTION: This function generates a volcano plot based on the DESeq2 results when the user inputs a `targets` parameter (specific for a targeted analysis).
  ###

  volcano_plot <- ggplot2::ggplot() + ggplot2::theme_bw() +
    ggplot2::geom_point(data = data, mapping = ggplot2::aes(x = .data$log2_FC, y = .data$neg_log10_adjusted_p_value), color = "black", alpha = 0.5) +
    ggplot2::geom_point(data = dplyr::filter(data, abs(.data$log2_FC) > fc_threshold, .data$neg_log10_adjusted_p_value > -log10(pval_threshold)),
                        mapping = ggplot2::aes(x = .data$log2_FC, y = .data$neg_log10_adjusted_p_value), size = 3) # +

  if (isTRUE(sig_axis)){ # Highlight only the significant data points
    volcano_plot <- volcano_plot +
      ggrepel::geom_text_repel(data = dplyr::filter(data, abs(.data$log2_FC) > fc_threshold, .data$neg_log10_adjusted_p_value > -log10(pval_threshold)),
                               mapping = ggplot2::aes(x = .data$log2_FC, y = .data$neg_log10_adjusted_p_value, label = .data$feature), size = 2) # +
  } else { # Limit the axis based on the min and max values of the input data
    x_range <- range(data$log2_FC, na.rm = TRUE)
    y_range <- range(data$neg_log10_adjusted_p_value, na.rm = TRUE)
    volcano_plot <- volcano_plot + ggplot2::coord_cartesian(xlim = x_range, ylim = y_range)
  }


  volcano_plot <- volcano_plot +
    ggrepel::geom_text_repel(data = data, mapping = ggplot2::aes(x = .data$log2_FC, y = .data$neg_log10_adjusted_p_value, label = .data$feature), size = 2) +
    ggplot2::geom_abline(slope = 0, intercept = -log10(pval_threshold), linetype = "dotted") +
    ggplot2::geom_vline(xintercept = -1 * fc_threshold, linetype = "dotted") +
    ggplot2::geom_vline(xintercept = fc_threshold, linetype = "dotted") +
    ggplot2::labs(x = expression(log[2]~fold~change), y = expression(-log[10]~adjusted~p~value), title = "target group vs. others") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20))

  return(volcano_plot)
}

HeatmapFontSize <- function(data, base_size = 10, shrink_factor = 1.5, min_size = 2) {

  ###
  # CALL: ProduceHeatmapDiffExp()
  # DESCRIPTION: This function takes into consideration the number of features in `data` to automatically adjust the size of the labels.
  # The main objective is to avoid a messy display in the plot. Returns the size of the labels in the rows and columns
  ###

  fontsize_row <- max(min_size, base_size - shrink_factor * log10(nrow(data)))
  fontsize_col <- max(min_size, base_size - shrink_factor * log10(ncol(data)))

  return(list(fontsize_row = fontsize_row, fontsize_col = fontsize_col))
}

ProduceHeatmapDiffExp <- function(data){

  ###
  # CALL: DESeq2runner()
  # DESCRIPTION: This function generates a heatmap to represent the hierarchical relation of the data.
  ###

  Euc_Dists <- stats::dist(t(data)) # Compute Euclidean distances

  # Sort nodes in hierarchical clustering
  sort_hclust <- function(...) stats::as.hclust(dendsort::dendsort(stats::as.dendrogram(...)))
  mat_cluster_names <- sort_hclust(stats::hclust(Euc_Dists))

  font_size <- HeatmapFontSize(data = data) # Control the size of the labels in the heatmap

  # Produce heatmap that has been hierarchically clustered
  Dist_heatmap <- pheatmap::pheatmap(Euc_Dists, cluster_cols = mat_cluster_names, cluster_rows = mat_cluster_names, color = viridis::inferno(1000),
                                     border_color = NA, labels_col = mat_cluster_names$labels, labels_row = mat_cluster_names$labels,
                                     fontsize_row = font_size$fontsize_row, fontsize_col = font_size$fontsize_col) # Assigned the computed sizes for the labels
  return(Dist_heatmap)
}

ProduceElbowPlot <- function(data, variance){

  ###
  # CALL: RunPCA()
  # DESCRIPTION: This function generates an elbow plot, to select the number of principal components (PCs) worth to further evaluate.
  ###

  # Generates a table with each PC and the variance they explain
  elbow_tibble <- tibble::tibble( PC = 1:length(data$sdev), variance_explained = variance)

  # Generate plot - reorder() to add PCs in the x-axis and variance in the y-axis
  elbow_plot <- ggplot2::ggplot(data = elbow_tibble) + ggplot2::theme_bw() +
    ggplot2::geom_bar(mapping = ggplot2::aes(x = stats::reorder(.data$PC, -.data$variance_explained),  y = .data$variance_explained), stat = "identity") +
    ggplot2::labs(title = "Percentage of variance explained per PC", x = "PC", y = "% variance explained") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

  return(elbow_plot) # Outputs the elbow plot
}
