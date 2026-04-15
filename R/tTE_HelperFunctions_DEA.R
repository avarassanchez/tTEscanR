ComputeSizeCorrection <- function(data, metadata, batch = NULL, reduce = 100, verbose = TRUE){

  ###
  # CALL: Multiple functions inside tTEscanR
  # DESCRIPTION: This function uses the library DESeq2 to size correct the data to account for sequencing depth.
  ###

  CheckDataFrame(data = data) # Evaluate the input parameter: 'data'
  CheckDataFrame(data = metadata) # Evaluate the input parameters: 'metadata' and 'batch'
  if (!(batch %in% colnames(metadata))) stop("The correction factor ('batch') was not found in the 'metadata'.")
  if (verbose) message("- The 'data', 'metadata' and 'batch' parameters have been properly loaded.\n",
                       "- Size-correcting the counts using DESeq2.\n", "- Filtering (if necessary) the 'data' and 'metadata' for common features.")

  filtered <- FilterByMetadata(data = data, metadata = metadata, verbose = verbose) # Check consistency in conditions included in the data and metadata
  DESeq2_run <- suppressWarnings(ComputeDESeq2(data = filtered$data, metadata = filtered$metadata, batch = batch, reduce = reduce))

  if (verbose) message("- Extracting the size-corrected counts.")
  size_corrected_output_matrix <- suppressMessages(DESeq2::counts(DESeq2_run, normalized = TRUE)) # Extracting the normalized data from the DESeq2 object
  CheckDataFrame(data = size_corrected_output_matrix)
  return(size_corrected_output_matrix) # Return the size-corrected matrix
}

ComputeAllPairwiseComp <- function(dds, factor_name, padj_threshold = 0.05) {

  ###
  # CALL: ComputeDEResults()
  # DESCRIPTION: This function takes the fitted model and extracts the computed statistics.
  ###

  if (!factor_name %in% colnames(SummarizedExperiment::colData(dds))) stop("Factor '", factor_name, "' not found in colData.")
  factor_data <- SummarizedExperiment::colData(dds)[[factor_name]]
  if (!is.factor(factor_data)) stop("The column '", factor_name, "' must be a factor in the DESeqDataSet.")
  levels_factor <- levels(factor_data)
  if (length(levels_factor) < 2) stop("Factor must have at least 2 levels to compute contrasts.")

  combn_matrix <- utils::combn(levels_factor, 2) # Generate all pairwise combinations
  results_list <- list()
  for (i in seq_len(ncol(combn_matrix))) {
    level1 <- combn_matrix[1, i]
    level2 <- combn_matrix[2, i]
    contrast_name <- paste0(level1, "_vs_", level2)

    message("Extracting contrasts: ", contrast_name)
    res <- DESeq2::results(dds, contrast = c(factor_name, level1, level2), alpha = padj_threshold) # Alpha marks the significance threshold
    results_list[[contrast_name]] <- res
  }
  return(results_list)
}

ComputeDESeq2 <- function(data, metadata, condition = NULL, batch = NULL, reference = NULL, reduce = 100) {

  ###
  # CALL: RunDEAnalysis()
  # DESCRIPTION: Prepares the matrix, validates experimental factors, safely handles empty samples, builds a dynamic design formula, and runs the DESeq2 model.
  ###

  data <- as.matrix(data)
  if (any(data > .Machine$integer.max)) { # R has a limit in the length of an integer that can be analyzed.
    data <- round(data / reduce)
    message("- Some values exceed R's integer limit. Matrix divided by ", reduce, ".")
  }
  storage.mode(data) <- "integer" # Ensure the data is integer even if it came in as double

  sample_sums <- Matrix::colSums(data)
  valid_samples <- sample_sums > 0 # Remove empty samples

  if (sum(!valid_samples) > 0) {
    message("- Removing ", sum(!valid_samples), " samples with zero total counts to prevent DESeq2 crash.")
    data <- data[, valid_samples, drop = FALSE]
    metadata <- metadata[valid_samples, , drop = FALSE]
  }

  colData <- as.data.frame(metadata) # Ensures a standard data.frame
  if (ncol(data) != nrow(colData)) stop("Dimension mismatch: The count matrix has ", ncol(data), " samples, but the metadata has ", nrow(colData), " rows.")

  if (!is.null(condition)) {
    if (!condition %in% colnames(colData)) stop("The condition factor '", condition, "' was not found in the metadata.")
    colData[[condition]] <- factor(colData[[condition]])

    if (!is.null(reference)) { # Apply reference level to the primary condition if provided
      if (!(reference %in% levels(colData[[condition]]))) stop("The reference level '", reference, "' was not found in the 'condition' column.")
      colData[[condition]] <- stats::relevel(colData[[condition]], ref = reference)
      message("- Set '", reference, "' as the reference level for '", condition, "'.")
    }
  }

  if (!is.null(batch)) {
    if (!batch %in% colnames(colData)) stop("The batch factor '", batch, "' was not found in the metadata.")
    colData[[batch]] <- factor(colData[[batch]])
  }

  if (!is.null(batch) && !is.null(condition)) { # Dynamic definition of the design formula
    if (batch == condition) {
      design_formula <- stats::as.formula(paste("~", condition))
      message("- Batch and Condition are identical. Building single-factor model: ~ ", condition)
    } else {
      design_formula <- stats::as.formula(paste("~", batch, "+", condition))
      message("- Building model with batch correction: ~ ", batch, " + ", condition)
    }
  } else if (!is.null(condition)) {
    design_formula <- stats::as.formula(paste("~", condition))
    message("- Building model: ~ ", condition)
  } else if (!is.null(batch)) {
    design_formula <- stats::as.formula(paste("~", batch))
    message("- Building exploratory model based on batch: ~ ", batch)
  } else {
    design_formula <- stats::as.formula("~ 1")
    message("- No factors provided. Building intercept-only exploratory model: ~ 1")
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = data, colData = colData, design = design_formula)

  message("- Running DESeq2 pipeline...")
  dds <- DESeq2::DESeq(dds)

  return(dds)
}

TargetedApproach <- function(DESeq2_run, verbose = TRUE, fc_threshold, padj_threshold, label_significant, condition) {

  if (verbose) message("- Executing a targeted analysis based on '", condition, "'.")
  col_data <- SummarizedExperiment::colData(DESeq2_run)
  levels_factor <- levels(col_data[[condition]])
  combs <- utils::combn(levels_factor, 2, simplify = FALSE)
  results.list <- list()

  if (verbose) message("    Levels found in 'condition' factor: ", paste(levels_factor, collapse = ", "), "\n    Number of pairwise combinations: ", length(combs))
  for (cmb in combs) {
    level1 <- cmb[1]
    level2 <- cmb[2]
    contrast_name <- paste(level1, "vs", level2, sep = "_")

    if (verbose) message("Processing contrast: ", contrast_name)
    res <- DESeq2::results(DESeq2_run, contrast = c(condition, level1, level2))
    res_df <- as.data.frame(res)

    target_results_tibble <- tibble::tibble(feature = rownames(res_df), log2_FC = res_df$log2FoldChange, p_value = res_df$pvalue, neg_log10_adjusted_p_value = -log10(res_df$padj)) %>%
      dplyr::filter(!is.na(.data$log2_FC) & !is.na(.data$neg_log10_adjusted_p_value))

    results.list[[paste0("results_", contrast_name)]] <- target_results_tibble
    if (verbose) message("- Generating volcano plot.")
    volcano_plot <- GenerateVolcanoPlot(data = target_results_tibble, fc_threshold = fc_threshold, padj_threshold = padj_threshold, label_significant = label_significant) +
      ggplot2::ggtitle(contrast_name)
    results.list[[paste0("volcano_", contrast_name)]] <- volcano_plot
  }

  return(results.list)
}

MakeScatterPlot <- function(df, x_col, y_col, color_factor, shape_factor = NULL, label_factor = NULL, title = NULL, color_palette = NULL,
                            highlight_median = FALSE, point_alpha = 0.5, median_size = 6, x_axis = NULL, y_axis = NULL, show_legend = "none") {

  plot <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]]))

  if (!is.null(shape_factor)) {
    plot <- plot + ggplot2::geom_point(ggplot2::aes(color = .data[[color_factor]], shape = .data[[shape_factor]]), alpha = if (highlight_median) point_alpha else 1, size = 3)
  } else {
    plot <- plot + ggplot2::geom_point(ggplot2::aes(color = .data[[color_factor]]), alpha = if (highlight_median) point_alpha else 1, size = 3)
  }

  if (highlight_median) { # Median points and labels
    medians <- df %>% dplyr::group_by(.data[[color_factor]]) %>%
      dplyr::summarise(dplyr::across(.cols = dplyr::all_of(c(x_col, y_col)), .fns  = ~ stats::median(.x, na.rm = TRUE)), .groups = "drop")

    plot <- plot +
      ggplot2::geom_point(data = medians, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[color_factor]]), size = median_size, inherit.aes = FALSE) +
      ggrepel::geom_text_repel(data = medians, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], label = .data[[color_factor]]), fontface = "bold", size = 3, show.legend = FALSE, inherit.aes = FALSE)
  }

  if (!is.null(label_factor)) plot <- plot + ggrepel::geom_text_repel(ggplot2::aes(label = .data[[label_factor]]), size = 2.8, show.legend = FALSE)
  if (!is.null(color_palette)) plot <- plot + ggplot2::scale_color_manual(values = color_palette) else plot <- plot + ggplot2::scale_color_viridis_d(option = "viridis")

  plot <- plot + ggplot2::guides(color = ggplot2::guide_legend(), shape = ggplot2::guide_legend()) + ggplot2::theme_bw()
  plot <- plot + ggplot2::labs(title = title, color = color_factor, shape = shape_factor,
                               x = if (!is.null(x_axis)) sprintf("%s (%.1f%%)", x_col, x_axis) else x_col,
                               y = if (!is.null(y_axis)) sprintf("%s (%.1f%%)", y_col, y_axis) else y_col)
  plot <- plot + ggplot2::theme(legend.position = show_legend)
  return(plot)
}

RunDimReduct <- function(vst, metadata, condition = NULL, color_factor = condition, shape_factor = NULL, label_factor = NULL, show_legend = "none",
                         dim_reduct, highlight_median, numPC = 5, color_palette = NULL, scale_pca = FALSE) {

  plots <- list()
  mat <- t(vst) # Transpose vst to samples x features

  if (dim_reduct == "PCA") {
    # Remove zero or near-zero variance columns (features)
    low_var_cols <- which(matrixStats::colVars(mat, useNames = TRUE) < 1e-5)
    if (length(low_var_cols) > 0) {
      mat <- mat[, -low_var_cols, drop = FALSE]
      message("- Removed ", length(low_var_cols), " low-variance features before PCA.")
    }

    pca_res <- stats::prcomp(mat, scale. = scale_pca, center = TRUE) # Scale set to FALSE as the vst is already normalized
    coords <- pca_res$x[, 1:min(numPC, ncol(pca_res$x))]
    colnames(coords) <- paste0("PC", 1:ncol(coords))

    variance <- 100 * (pca_res$sdev^2 / sum(pca_res$sdev^2))
    plots[["ElbowPlot"]] <- ProduceElbowPlot(data = pca_res, variance = variance)
  }

  if (dim_reduct == "UMAP") {
    coords <- tryCatch({
      uwot::umap(mat)
    }, error = function(e) {
      warning("UMAP failed: ", e$message, "\nUsing PCA instead.")
      stats::prcomp(mat, scale. = TRUE, center = TRUE)$x[, 1:2]
    })
    colnames(coords) <- c("Dim1", "Dim2")
  }

  if (dim_reduct == "tSNE") {
    n_samples <- nrow(mat)
    dynamic_perp <- floor((n_samples - 1) / 3) # Dynamic calculation of the perplexity
    final_perp <- min(30, max(1, dynamic_perp))
    message("- Running t-SNE with perplexity = ", final_perp)

    coords <- Rtsne::Rtsne(mat, perplexity = final_perp, check_duplicates = FALSE)$Y
    colnames(coords) <- c("Dim1", "Dim2")
  }

  rownames(coords) <- rownames(mat) # Restore sample names before merging
  plot_df <- cbind(as.data.frame(coords), metadata)

  # Generate plots
  if (dim_reduct == "PCA") {
    pc_names <- colnames(coords)[seq_len(numPC)]
    pc_pairs <- utils::combn(pc_names, 2, simplify = FALSE)

    for (pair in pc_pairs) {
      x_pc <- pair[1]
      y_pc <- pair[2]
      plots[[paste0(x_pc, "_vs_", y_pc)]] <- MakeScatterPlot(df = plot_df, x_col = x_pc, y_col = y_pc, color_factor = color_factor, shape_factor = shape_factor,
                                                             label_factor = label_factor, title = sprintf("%s vs %s", x_pc, y_pc), x_axis = variance[which(pc_names == x_pc)],
                                                             y_axis = variance[which(pc_names == y_pc)], color_palette = color_palette, highlight_median = highlight_median, show_legend = show_legend)
    }

  } else {
    plots[[dim_reduct]] <- MakeScatterPlot(df = plot_df, x_col = colnames(coords)[1], y_col = colnames(coords)[2], color_factor = color_factor,
                                           shape_factor = shape_factor, label_factor = label_factor, title = dim_reduct,
                                           color_palette = color_palette, highlight_median = highlight_median, show_legend = show_legend)
  }

  return(plots)
}

GenerateVolcanoPlot <- function(data, fc_threshold, padj_threshold, label_significant){

  ###
  # CALL: DESeq2runner()
  # DESCRIPTION: This function generates a volcano plot based on the DESeq2 results when the user inputs a `targets` parameter (specific for a targeted analysis).
  ###

  sig_data <- data %>% dplyr::filter(abs(.data$log2_FC) > fc_threshold, .data$neg_log10_adjusted_p_value > -log10(padj_threshold)) # Extract the statistical significant data

  volcano_plot <- ggplot2::ggplot() + ggplot2::theme_bw() +
    ggplot2::geom_point(data = data, mapping = ggplot2::aes(x = .data$log2_FC, y = .data$neg_log10_adjusted_p_value), color = "black", alpha = 0.5) +
    ggplot2::geom_point(data = sig_data, mapping = ggplot2::aes(x = .data$log2_FC, y = .data$neg_log10_adjusted_p_value), color = "darkred", size = 3)

  if (isTRUE(label_significant)) {
    volcano_plot <- volcano_plot + ggrepel::geom_text_repel(data = sig_data, ggplot2::aes(x = .data$log2_FC, y = .data$neg_log10_adjusted_p_value, label = .data$feature), size = 2, max.overlaps = 20)
  }

  volcano_plot <- volcano_plot +
    ggplot2::geom_hline(yintercept = -log10(padj_threshold), linetype = "dotted") +
    ggplot2::geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dotted") +
    ggplot2::labs(x = expression(log[2]~fold~change), y = expression(-log[10]~adjusted~p~value)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20))

  return(volcano_plot)
}

HeatmapFontSize <- function(data, base_size = 10, shrink_factor = 1.5, min_size = 2) {

  ###
  # CALL: ProduceHeatmapDiffExp()
  # DESCRIPTION: Avoids messy displays by adjusting the size of the labels (output) based on the number of features.
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
  dist_matrix <- as.matrix(Euc_Dists)

  dend <- dendsort::dendsort(stats::as.dendrogram(stats::hclust(Euc_Dists))) # Sort nodes in hierarchical clustering
  font_size <- HeatmapFontSize(data = dist_matrix) # Control the size of the labels in the heatmap

  # Produce heatmap that has been hierarchically clustered
  heatmap_grop <- grid::grid.grabExpr({
    Dist_heatmap <- ComplexHeatmap::Heatmap(matrix = dist_matrix, name = "Euclidean\nDistance", col = viridis::inferno(100), cluster_columns = dend, cluster_rows = dend, show_row_dend = TRUE,
                                            show_column_dend = TRUE, border = FALSE, row_names_gp = grid::gpar(fontsize = font_size$fontsize_row), column_names_gp = grid::gpar(fontsize = font_size$fontsize_col))
    ComplexHeatmap::draw(Dist_heatmap)
  })

  return(heatmap_grop)
}

ProduceElbowPlot <- function(data, variance){

  ###
  # CALL: RunPCA()
  # DESCRIPTION: This function generates an elbow plot, to select the number of principal components (PCs) worth to further evaluate.
  ###

  numPC <- min(length(variance), 50) # In order to keep the output cleaner the number of PCs will be limited to 50
  elbow_tibble <- tibble::tibble(PC = 1:numPC, variance_explained = variance[1:numPC]) # Generates table with each PC and their variance

  # Generate plot - reorder() to add PCs in the x-axis and variance in the y-axis
  elbow_plot <- ggplot2::ggplot(data = elbow_tibble, ggplot2::aes(x = .data$PC, y = .data$variance_explained)) +
    ggplot2::theme_bw() + ggplot2::geom_col(fill = "steelblue", alpha = 0.7) +
    ggplot2::geom_line(mapping = ggplot2::aes(group = 1), color = "black", linetype = "dashed", linewidth = 0.8) +
    ggplot2::geom_point(color = "black", size = 2) +
    ggplot2::labs(title = "Elbow Plot: Variance Explained per PC", x = "Principal Component", y = "Variance Explained (%)") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14), axis.text.x = ggplot2::element_text(angle = 0))

  return(elbow_plot) # Outputs the elbow plot
}
