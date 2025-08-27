DESeq2runner <- function(data, metadata, corr_factor, targets = NULL, target_factor = NULL, fc_threshold = 1, pval_threshold = 0.05, reduce = 100,
                         heatmap = TRUE, PCA = TRUE, numPC = 2, color_factor = NULL, shape_factor = NULL, color_palette = NULL, labels = FALSE, verbose = TRUE){

  ###
  # CALL: ExecuteDESeq2runner()
  # DESCRIPTION: This function runs the DESeq2 analysis over an individual data input.
  # It will be called as many times as datasets included in the `list_data` parameter in ExecuteDESeq2runner().
  # The input parameters are checked in the main function ExecuteDESeq2runner().
  # Returns a list with the corresponding outcomes, depends on the analysis approach (exploratory or targeted).
  ###

  # Check the consistency in the conditions included in the data and metadata - If needed filters them accordingly
  if (verbose) message("- Filtering (if necessary) the `data` and `metadata` for common features.")
  filtered <- FilterByMetadata(data = data, metadata = metadata, verbose = verbose)
  data <- filtered[[1]]
  metadata_filtered <- filtered[[2]]

  results.list <- list() # Create empty list to store the results

  if (verbose) message("- Running differential expression analysis.")
  DESeq2_run <- ComputeDESeq2(data = data, metadata = metadata, corr_factor = corr_factor, targets = targets, target_factor = target_factor, reduce = reduce, verbose = verbose)
  if (verbose) message("- Differential expression analysis completed.")

  if (is.null(targets)){ # EXPLORATORY ANALYSIS (NO TARGETS)
    if (verbose) message("- Computing variance stabilization.")
    vst <- DESeq2::varianceStabilizingTransformation(DESeq2_run)
    vst <- SummarizedExperiment::assay(vst)
    CheckDataFrame(data = vst)

    if (heatmap == TRUE){
      if (verbose) message("- Generating heatmap with Euclidean distances.")
      Dist_heatmap <- ProduceHeatmapDiffExp(data = vst)
      Dist_heatmap <- grDevices::recordPlot()
      results.list <- append(results.list, Dist_heatmap)
    }

    if (PCA == TRUE){
      if (numPC < 2) message("- Invalid `numPC` given, the default `numPC = 2` will be used instead.")

      if (is.null(color_factor)){ # Select the color_factor to use in the PCA plot
        color_factor <- corr_factor
        message("- No `color_factor` has been given as input, the `corr_factor` will be used instead.")
      }

      # Check points:
      if (!is.character(color_factor)) color_factor <- as.character(color_factor)
      if (!(color_factor %in% colnames(metadata))) stop("Please specify a `color_factor` parameter present in the `metadata` used as input.")
      if(!is.null(shape_factor) && !(shape_factor %in% colnames(metadata))) stop("Please specify a `shape_factor` parameter present in the `metadata` used as input.")

      if (verbose) message("- Running principal component analysis with ", numPC, " principal components.")
      PCA_results <- RunPCA(data = vst, metadata = metadata, numPC = numPC, color_factor = color_factor, shape_factor = shape_factor,
                            color_palette = color_palette, labels = labels, verbose = verbose)

      results.list <- append(results.list, PCA_results)
    }

    if (verbose) message("- Computing size correction to account for sequencing depth.")
    size_corrected_output_matrix <- suppressMessages(DESeq2::counts(DESeq2_run, normalized = TRUE))
    CheckDataFrame(data = size_corrected_output_matrix)

    results.list <- append(results.list, size_corrected_output_matrix)

  } else{ # TARGETED ANALYSIS

    if (verbose) message("- Executing a targeted analysis.\n", "- In this analysis no heatmap or PCA plot will be produced.\n", "- A volcano plot will be generated instead.")

    # Converts the format of the ComputeDESeq2() output
    target_DESeq2_results <- as.data.frame(DESeq2::results(DESeq2_run[[1]], contrast = DESeq2_run[[2]])) # Retrieve the different information parts from ComputeDESeq2()

    # Generate a summary table with the statistical measures
    target_results_tibble <- tibble::tibble(feature = rownames(data), log2_FC = target_DESeq2_results$log2FoldChange,
                                            neg_log10_adjusted_p_value = -log10(target_DESeq2_results$padj)) %>%
      dplyr::filter(!is.na(.data$log2_FC) & !is.na(.data$neg_log10_adjusted_p_value))  # Remove NA values

    results.list <- append(results.list, target_results_tibble)

    if (verbose) message("- Generating volcano plot.")
    volcano_plot <- GenerateVolcanoPlot(data = target_results_tibble, fc_threshold = fc_threshold, pval_threshold = pval_threshold)
    results.list <- append(results.list, volcano_plot)
  }

  return(results.list) # Return list with the resuts based on the analysis approach
}

ComputeDESeq2 <- function(data, metadata, corr_factor, targets = NULL, target_factor = NULL, reduce = 100, verbose = TRUE){

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

  if (!is.null(targets)){ # TARGETED APPROACH
    colData$class <- "other" # Add to the colData a new column "class"
    contrast_vect <- c("class", "other")

    for (i in 1:nrow(targets)){ # Iterates over each entry (row) in the `targets` table
      colData$class[grep(targets[i,1], colData[[target_factor]])] <- targets[i,2] # Update "class" based on `targets`
      if (targets[i,2] %in% contrast_vect) contrast_vect <- contrast_vect else contrast_vect <- append(contrast_vect, targets[i,2])
    }

    # It is necessary to have a target to compare against the rest of the conditions - The rest of the conditions will be labeled as "other"
    if (length(unique(colData$class)) < 2) stop("A single class was identified.\n", paste("Current `target_factor` = ", target_factor, ".\n"),
                                                paste("Class distribution:", as.character(length(unique(colData$class)))))
    if (verbose) message(paste("- Class distribution:", as.character(length(unique(colData$class)))))
    corr_factor <- "class"
  }

  formula <- stats::as.formula(paste("~", corr_factor)) # Define the formula to use by DESeq2 to correct the `data` (`corr_factor`)

  # Run DESeq2: DESeqDataSetFromMatrix() and DESeq()
  DESeq2_run <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(countData = data, colData = colData, design = formula))
  if (!is.null(targets)) message("- Executing a targeted analysis.\n", paste("- The class taken as reference for DESeq2 is:", levels(DESeq2_run$class)[[1]]))
  DESeq2_run <- suppressMessages(DESeq2::DESeq(DESeq2_run, quiet = TRUE))

  if (is.null(targets)) return(DESeq2_run) else return(list(DESeq2_run, contrast_vect)) # When performing a targeted approach, report the table with the amount of entries per condition
}

GenerateVolcanoPlot <- function(data, fc_threshold, pval_threshold){

  ###
  # CALL: DESeq2runner()
  # DESCRIPTION: This function generates a volcano plot based on the DESeq2 results.
  # It is executed when the user inputs a `targets` parameter, therefore it is specific for a targeted analysis.
  ###

  volcano_plot <- ggplot2::ggplot() + ggplot2::theme_classic() +
    ggplot2::geom_point(data = data, mapping = ggplot2::aes(x = .data$log2_FC, y = .data$neg_log10_adjusted_p_value), color = "black", alpha = 0.5) +
    ggplot2::geom_point(data = dplyr::filter(data, abs(.data$log2_FC) > fc_threshold, .data$neg_log10_adjusted_p_value > -log10(pval_threshold)),
                        mapping = ggplot2::aes(x = .data$log2_FC, y = .data$neg_log10_adjusted_p_value), size = 3) +
    ggrepel::geom_text_repel(data = dplyr::filter(data, abs(.data$log2_FC) > fc_threshold, .data$neg_log10_adjusted_p_value > -log10(pval_threshold)),
                             mapping = ggplot2::aes(x = .data$log2_FC, y = .data$neg_log10_adjusted_p_value, label = .data$feature), size = 2) +
    ggplot2::geom_abline(slope = 0, intercept = -log10(pval_threshold), linetype = "dotted") +
    ggplot2::geom_vline(xintercept = -1 * fc_threshold, linetype = "dotted") +
    ggplot2::geom_vline(xintercept = fc_threshold, linetype = "dotted") +
    ggplot2::labs(x = expression(log[2]~fold~change), y = expression(-log[10]~adjusted~p~value), title = "brain neurons vs. other cell types") +
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

ProducePCAPlot <- function(data, variance, metadata, numPC, color_factor, shape_factor, color_palette, labels){

  ###
  # CALL: RunPCA()
  # DESCRIPTION: This function generates a PCA plot for each PC given (returns a list).
  # The plots are customized based on the user's input parameters.
  ###

  outputs_PCA <- list() # Define empty list to store the outputs (PCA plots)

  # Based on the metadata generates a table including the PCs that will be interrogated
  PC_columns <- lapply(1:numPC, function(i) data$x[, i])
  names(PC_columns) <- paste0("PC_", 1:numPC)

  PCA_plot_tibble <- tibble::tibble(metadata, !!!PC_columns )
  CheckDataFrame(data = PCA_plot_tibble)

  # Selects the variable that will be used to color the data points
  index <- which(colnames(PCA_plot_tibble) == color_factor)
  PCA_plot_tibble$factor <- factor(PCA_plot_tibble[[index]])

  id_factor <- colnames(PCA_plot_tibble)

  for (i in 2:numPC){ # Iterates over all the PCs - starts in 2 as PC = 1 is always included (comparison in pairs)
    if(is.null(shape_factor)){ # Shape is NOT considered
      plot <- ggplot2::ggplot(data = PCA_plot_tibble, mapping =  ggplot2::aes_string(x = "PC_1", y = paste("PC", i, sep = "_"), label = "id_factor", color = color_factor))
    } else { # Shape is considered
      plot <- ggplot2::ggplot(data = PCA_plot_tibble, mapping =  ggplot2::aes_string(x = "PC_1", y = paste("PC", i, sep = "_"), label = "id_factor", color = color_factor, shape = shape_factor))
    }

    # Customization:
    if (!is.null(color_palette)) plot <- plot + ggplot2::scale_color_manual(values = color_palette) # Usage of user-defined color palette
    plot <- plot + ggplot2::theme_bw() + ggplot2::geom_point() + ggplot2::labs(x = sprintf("PC1 (%.1f%%)", variance[1]), y = sprintf("PC%d (%.1f%%)", i, variance[i])) # Add the points and the titles

    # Update the output list by adding the generated plot
    if (isFALSE(labels)){
      outputs_PCA[[paste0("PCA_plot_PC", i)]] <- plot
    } else{ # Include the names of the data points
      plot <- plot + ggrepel::geom_text_repel(ggplot2::aes(label = id_factor), size = 2, show.legend = FALSE)
      outputs_PCA[[paste0("PCA_plot_labels_PC", i)]] <- plot
    }
  }

  return(outputs_PCA) # Returns list of plots
}

RunPCA <- function(data, metadata, numPC, color_factor, shape_factor, color_palette, labels, verbose){

  ###
  # CALL: DESeq2runner()
  # DESCRIPTION: This function runs a principal component analysis (PCA) and generates an elbow plot and PCA plots.
  # Returns a nested list with the two outputs mentioned above as independent lists.
  ###

  # Check that the input data has a suitable format
  vst_matrix <- as.matrix(data)
  CheckDataFrame(data = vst_matrix)

  # Generate transpose of feature by condition matrix
  condition_by_feature_matrix <- t(vst_matrix)
  CheckDataFrame(data = condition_by_feature_matrix)

  # Remove features with too low variances, if any
  too_low_variance_features <- which(matrixStats::colVars(condition_by_feature_matrix, useNames = FALSE) < 1E-5)

  if (length(too_low_variance_features) > 1) {
    condition_by_feature_matrix <- condition_by_feature_matrix[ , -c(too_low_variance_features)] # Filtering the matrix to remove low variance features
    CheckDataFrame(data = condition_by_feature_matrix)
  }

  if (verbose) message("- Running `prcomp()` with `scale = TRUE` and `center = TRUE`.")
  PCA_res <- stats::prcomp(condition_by_feature_matrix, scale = TRUE, center = TRUE) # Run PCA

  PC_variance_explained = 100 * (PCA_res$sdev ^ 2) / sum(PCA_res$sdev ^ 2) # Examine amount of variance explained

  if (verbose) message("- Producing elbow plot.")
  elbow_plot <- ProduceElbowPlot(data = PCA_res, variance = PC_variance_explained)

  if (verbose) message("- Producing PCA plots.")
  pca_plots <- ProducePCAPlot(data = PCA_res, variance = PC_variance_explained, metadata = metadata, numPC = numPC,
                              color_factor = color_factor, shape_factor = shape_factor, color_palette = color_palette, labels = labels)

  return(list(elbow_plot = elbow_plot, pca_plots = pca_plots)) # Returns a list with: (i) elbow plot, (ii) list with a PCA plot per each PC considered
}
