#' Perform Differential Expression Analysis Using DESeq2
#' @description
#' This function applies differential expression analysis using the DESeq2 framework on a matrix of expression values.
#' It supports both exploratory visualizations (heatmap and PCA) and targeted comparisons using a custom contrast table.
#'
#' @param list_data A \code{list} of matrices with features (i.e. genes, codons, anticodons or amino acids) as rows and samples or conditions as columns. The list can be named or unnamed.
#' @param metadata A \code{data.frame} with the data associated to the conditions in the matrices of \code{list_data}. There has to be one column with the same labels as the column names.
#' @param dim_reduct Either \code{"PCA"}, \code{"UMAP"} or \code{"tSNE"} to specify the dimensionality reduction approach to be executed. Defaults to \code{"PCA"}.
#' @param batch Optional; name of the categorical variable in \code{metadata} to correct the data.
#' @param reference Optional; factor from the \code{batch} variable to use as reference for the corrections. If not specified, the 1st factor that appears will be used instead.
#' @param condition Optional; a factor based on \code{metadata} columns to define the comparisons to perform in a targeted analysis. Defaults to \code{NULL}.
#' @param fc_threshold Numeric; fold change threshold used for highlighting significant features in the volcano plot (if \code{targets} is specified). Defaults to 1.
#' @param padj_threshold Numeric; p-value threshold used for highlighting significant features in the volcano plot. Defaults to 0.05.
#' @param label_significant Logical; if \code{TRUE} displays the axis of the plots based on \code{fc_threshold} and \code{padj_threshold}. Defaults to \code{TRUE}.
#' @param reduce Numeric; a scaling factor used to normalize large expression values that exceed R's handling capacity. Defaults to 100.
#' @param heatmap Logical; if \code{TRUE}, generates a heatmap for exploratory analysis. Not applicable if \code{targets} is specified. Defaults to \code{TRUE}.
#' @param highlight_median Logical; if \code{TRUE} the data points of each cluster will be summarized into the median. Defaults to \code{FALSE}.
#' @param numPC Numeric; number of principal components to include in the PCA analysis. Required if \code{dim_reduct} is \code{"PCA"}. Defaults to 2.
#' @param show_legend Either \code{"none"} (default), \code{"top"}, \code{"bottom"}, \code{"right"} or \code{"left"} to specify the position of the legend in the plot.
#' @param color_factor Optional; name of the categorical variable in \code{metadata} to group by color the data points. Used if \code{dim_reduct} specified.
#' @param shape_factor Optional; name of the categorical variable in \code{metadata} to group by shape the data points. Used if \code{dim_reduct} specified.
#' @param color_palette Optional; a \code{vector} of color codes to customize the plot appearance.
#' @param label_factor Optional; name of the categorical variable to label the data points. Used if \code{dim_reduct} specified.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return A \code{list} of outputs per each matrix in \code{list_data}, based on the enabled parameters: (i) exploratory plots (heatmap and/or PCA), (ii) targeted plot (volcano), and (iii) size-corrected \code{data}.
#' @export
#'
#' @examples
#' data(default_tTEscanR_tRNA_data, default_tTEscanR_metadata)
#' DE_analysis <- RunDEAnalysis(list_data = list(mRNA = default_tTEscanR_tRNA_data),
#'                              metadata = default_tTEscanR_metadata, batch = "tissue",
#'                              color_factor = "tissue")

RunDEAnalysis <- function(list_data, metadata, condition = NULL, batch = NULL, reference = NULL, reduce = 100, dim_reduct = NULL, color_factor = NULL, heatmap = TRUE,
                          shape_factor = NULL, label_factor = NULL, highlight_median = FALSE, numPC = 2, color_palette = NULL, fc_threshold = 1,
                          padj_threshold = 0.05, label_significant = TRUE, verbose = TRUE, show_legend = "none") {

  # Compute results
  results <- ComputeDEResults(list_data = list_data, metadata = metadata, condition = condition, batch = batch, reference = reference, reduce = reduce, verbose = verbose)

  # Plot all datasets
  plots <- list()
  for (dataset_name in names(results)) {
    plots[[dataset_name]] <- PlotDEResults(DE_results_list = results, dataset_name = dataset_name, dim_reduct = dim_reduct,
                                           color_factor = color_factor, shape_factor = shape_factor, label_factor = label_factor,
                                           highlight_median = highlight_median, numPC = numPC, color_palette = color_palette,
                                           condition = condition, fc_threshold = fc_threshold, padj_threshold = padj_threshold,
                                           heatmap = heatmap, label_significant = label_significant, verbose = verbose, show_legend = show_legend)
  }

  return(list(results = results, plots = plots))
}


#' Compute the DESeq2 Analysis
#'
#' @param list_data A \code{list} of matrices with features (i.e. genes, codons, anticodons or amino acids) as rows and samples or conditions as columns. The list can be named or unnamed.
#' @param metadata A \code{data.frame} with the data associated to the conditions in the matrices of \code{list_data}. There has to be one column with the same labels as the column names.
#' @param condition Optional; a factor based on \code{metadata} columns to define the comparisons to perform in a targeted analysis. Defaults to \code{NULL}.
#' @param batch Optional; name of the categorical variable in \code{metadata} to correct the data. Required if \code{dim_reduct} specified.
#' @param reference Optional; factor from the \code{batch} variable to use as reference for the corrections. If not specified, the 1st factor that appears will be used instead.
#' @param padj_threshold Numeric; p-value threshold used for highlighting significant features in the volcano plot. Defaults to 0.05.
#' @param reduce Numeric; a scaling factor used to normalize large expression values that exceed R's handling capacity. Defaults to 100.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @returns A DESeq2 object with the normalized and vst counts.
#' @export
#'
#' @examples
#' data(default_tTEscanR_tRNA_data, default_tTEscanR_metadata)
#' DE_analysis <- ComputeDEResults(list_data = list(tRNA = default_tTEscanR_tRNA_data),
#'                                 metadata = default_tTEscanR_metadata, batch = "tissue")

ComputeDEResults <- function(list_data, metadata, condition = NULL, batch = NULL, reference = NULL, reduce = 100, padj_threshold = 0.05, verbose = TRUE){

  DE_results_list <- list()

  # Determine the correction factor
  if (is.null(condition) && is.null(batch)) stop("Please specify the 'condition' or the 'batch'.")
  cf <- if (!is.null(batch)) batch else condition

  # Check the names of the datasets in list_data
  if (is.null(names(list_data))){
    names(list_data) <- paste0("dataset_", seq_along(list_data))
  } else {
    empty_names <- names(list_data) == ""
    if (any(empty_names))names(list_data)[empty_names] <- paste0("dataset_", which(empty_names))
  }

  for (i in seq_along(list_data)) {
    dataset_name <- names(list_data)[i]
    message(paste("- Processing dataset:", names(list_data)[i]))

    filtered <- FilterByMetadata(data = list_data[[i]], metadata = metadata, verbose = verbose)

    DESeq2_run <- suppressWarnings(ComputeDESeq2(data = filtered$data, metadata = filtered$metadata, batch = cf, condition = condition, reference = reference, reduce = reduce))
    pairwise_results <- ComputeAllPairwiseComp(dds = DESeq2_run, factor_name = cf, padj_threshold = padj_threshold) # Run all pairwise comparisons

    size_corrected <- suppressMessages(DESeq2::counts(DESeq2_run, normalized = TRUE)) # Extract normalized counts and vst
    vst <- SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(DESeq2_run))

    DE_results_list[[names(list_data)[i]]] <- list(DESeq2_run = DESeq2_run, pairwise_results = pairwise_results, size_corrected = size_corrected, vst = vst, metadata = filtered$metadata) # Store results and metadata
  }

  return(DE_results_list)
}

#' Generates Visualizations from the DEA data
#'
#' @param DE_results_list A DESeq2 object with the normalized and vst counts. Can be obtained by running \code{\link{ComputeDEResults}}
#' @param dataset_name String to specify the assay in \code{DE_results_list} to extract.
#' @param dim_reduct Either \code{"PCA"}, \code{"UMAP"} or \code{"tSNE"} to specify the dimensionality reduction approach to be executed. Defaults to \code{"PCA"}.
#' @param color_factor Optional; name of the categorical variable in \code{metadata} to group by color the data points. Used if \code{dim_reduct} specified.
#' @param shape_factor Optional; name of the categorical variable to label the data points. Used if \code{dim_reduct} specified.
#' @param label_factor Optional; name of the categorical variable to label the data points. Used if \code{dim_reduct} specified.
#' @param highlight_median  Logical; if \code{TRUE} the data points of each cluster will be summarized into the median. Defaults to \code{FALSE}.
#' @param numPC Numeric; number of principal components to include in the PCA analysis. Required if \code{dim_reduct} is \code{"PCA"}. Defaults to 2.
#' @param scale_pca Logical; if \code{TRUE}, scales the data if \code{dim_reduct} is \code{"PCA"}. Defaults to \code{FALSE}.
#' @param heatmap Logical; if \code{TRUE}, generates a heatmap for exploratory analysis. Defaults to \code{TRUE}.
#' @param color_palette Optional; a \code{vector} of color codes to customize the plot appearance.
#' @param show_legend Either \code{"none"} (default), \code{"top"}, \code{"bottom"}, \code{"right"} or \code{"left"} to specify the position of the legend in the plot.
#' @param condition Optional; a factor based on \code{metadata} columns to define the comparisons to perform in a targeted analysis. Defaults to \code{NULL}.
#' @param fc_threshold Numeric; fold change threshold used for highlighting significant features in the volcano plot (if \code{targets} is specified). Defaults to 1.
#' @param padj_threshold Numeric; p-value threshold used for highlighting significant features in the volcano plot. Defaults to 0.05.
#' @param label_significant Logical; if \code{TRUE} displays the axis of the plots based on \code{fc_threshold} and \code{padj_threshold}. Defaults to \code{TRUE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @returns Visualization of the DEA results. The "heatmap" element is a \code{gtable}/\code{grob}
#' produced by \code{ComplexHeatmap} and captured via \code{grid::grid.grabExpr}.
#' @export
#'
#' @examples
#' data(default_tTEscanR_tRNA_data, default_tTEscanR_metadata)
#' DE_analysis <- ComputeDEResults(list_data = list(tRNA = default_tTEscanR_tRNA_data),
#'                                 metadata = default_tTEscanR_metadata, batch = "tissue")
#' DE_plots <- PlotDEResults(DE_results_list = DE_analysis, dataset_name = "tRNA",
#'                           dim_reduct = "PCA", condition = "tissue", heatmap = FALSE)

PlotDEResults <- function(DE_results_list, dataset_name = NULL, dim_reduct = NULL, condition = NULL, color_factor = condition, shape_factor = NULL,
                          label_factor = NULL, highlight_median = FALSE, numPC = 2, scale_pca = FALSE, color_palette = NULL, heatmap = TRUE,
                          fc_threshold = 1, padj_threshold = 0.05, label_significant = TRUE, show_legend = "none", verbose = TRUE) {

  if (is.null(dataset_name)) dataset_name <- names(DE_results_list)[1]
  if (!dataset_name %in% names(DE_results_list)) stop("dataset_name '", dataset_name, "' not found in 'DE_results_list'.")

  vst <- DE_results_list[[dataset_name]]$vst # Extract precomputed data
  results_list <- list()

  if (isTRUE(heatmap)) {
    if (verbose) message("- Generating heatmap with Euclidean distances.")
    heatmap_plot <- ProduceHeatmapDiffExp(data = vst)
    results_list[["heatmap"]] <- heatmap_plot
  }

  # Exploratory Approach
  if (!is.null(dim_reduct)){
    dim_reduct <- match.arg(dim_reduct, choices = c("PCA", "UMAP", "tSNE"))

    if (dim_reduct != "PCA"){
      numPC <- NULL
    } else if (is.null(numPC) || numPC < 2) {
        message("- Invalid 'numPC' given, using default 'numPC = 2'.")
        numPC <- 2
    }

    metadata <- DE_results_list[[dataset_name]]$metadata

    if (is.null(color_factor)) stop("Please specify a suitable 'color_factor'.")
    if (length(color_factor) > 1 || length(shape_factor) > 1 || length(label_factor) > 1) stop("The parameters 'color_factor', 'shape_factor' and 'label_factor' are expected to be the name of a column in the `metadata`.")
    if (any(c(length(color_factor), length(shape_factor), length(label_factor)) > 1)) stop("The parameters 'color_factor', 'shape_factor' and 'label_factor' must be a single column name, not a vector")

    provided_factors <- c(color_factor, shape_factor, label_factor) # If a factor is set to NULL it will not be added to the vector
    missing_factors <- setdiff(provided_factors, colnames(metadata))
    if (length(missing_factors) > 0) stop("The following factors are not found in the metadata: ", paste(missing_factors, collapse = ","))

    if (verbose) message("- Generating dimensionality reduction plots.")
    exploratory_plots <- RunDimReduct(metadata = metadata, vst = vst, highlight_median = highlight_median, show_legend = show_legend,
                                      dim_reduct = dim_reduct, numPC = numPC, color_factor = color_factor, shape_factor = shape_factor,
                                      label_factor = label_factor, color_palette = color_palette, scale_pca = scale_pca)
    results_list[["exploratory"]] <- exploratory_plots
  }

  # Targeted Approach
  if (!is.null(condition)) {
      DESeq2_run <- DE_results_list[[dataset_name]]$DESeq2_run
      if (verbose) message("- Generating targeted analysis (volcano plots).")
      targeted_results <- TargetedApproach(DESeq2_run = DESeq2_run, condition = condition, fc_threshold = fc_threshold,
                                           padj_threshold = padj_threshold, label_significant = label_significant, verbose = verbose)
      results_list[["targeted"]] <- targeted_results
  }

  return(results_list)
}

