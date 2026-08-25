#' Perform Differential Expression Analysis Using DESeq2
#' @description
#' This function applies differential expression analysis using the DESeq2
#' framework on a matrix of expression values. It supports both exploratory
#' visualizations (heatmap and PCA) and targeted comparisons using a custom
#' contrast table.
#'
#' @param list_data A \code{list} of matrices with features (i.e. genes, codons,
#'     anticodons or amino acids) as rows and samples or conditions as columns.
#'     The list can be named or unnamed.
#' @param metadata A \code{data.frame} with the data associated to the
#'     conditions in the matrices of \code{list_data}. There has to be one
#'     column with the same labels as the column names.
#' @param dim_reduct Either \code{"PCA"}, \code{"UMAP"} or \code{"tSNE"} to
#'     specify the dimensionality reduction approach to be executed. Defaults
#'     to \code{"PCA"}.
#' @param batch Optional; name of the categorical variable in \code{metadata} to
#'     correct the data.
#' @param reference Optional; factor from the \code{batch} variable to use as
#'     reference for the corrections. If not specified, the 1st factor that
#'     appears will be used instead.
#' @param target Optional; a factor based on \code{metadata} columns to
#'     define the comparisons to perform in a targeted analysis. Defaults to
#'      \code{NULL}.
#' @param fc_threshold Numeric; fold change threshold used for highlighting
#'     significant features in the volcano plot. Defaults to 1.
#' @param padj_threshold Numeric; p-value threshold used for highlighting
#'     significant features in the volcano plot. Defaults to 0.05.
#' @param label_significant Logical; if \code{TRUE} displays the axis of the
#'     plots based on \code{fc_threshold} and \code{padj_threshold}. Defaults
#'     to \code{TRUE}.
#' @param compute_pairwise Logical; if \code{TRUE}, computes all the pairwise
#'     comparisons based on the conditions included in the input data.
#' @param reduce Numeric; a scaling factor used to normalize large expression
#'     values that exceed R's handling capacity. Defaults to 100.
#' @param heatmap Logical; if \code{TRUE}, generates a heatmap for exploratory
#'     analysis. Defaults to \code{TRUE}.
#' @param highlight_median Logical; if \code{TRUE} the data points of each
#'     cluster will be summarized into the median. Defaults to \code{FALSE}.
#' @param numPC Numeric; number of principal components to include in the PCA
#'     analysis. Required if \code{dim_reduct} is \code{"PCA"}. Defaults to 2.
#' @param show_legend Either \code{"none"} (default), \code{"top"},
#'     \code{"bottom"}, \code{"right"} or \code{"left"} to specify the
#'     position of the legend in the plot.
#' @param color_factor Optional; name of the categorical variable in
#'     \code{metadata} to group by color the data points. Used if
#'     \code{dim_reduct} specified.
#' @param shape_factor Optional; name of the categorical variable in
#'     \code{metadata} to group by shape the data points. Used if
#'     \code{dim_reduct} specified.
#' @param color_palette Optional; a \code{vector} of color codes to customize
#'     the plot appearance.
#' @param label_factor Optional; name of the categorical variable to label the
#'     data points. Used if \code{dim_reduct} specified.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#'     Defaults to \code{TRUE}.
#'
#' @return A \code{list} of outputs per each matrix in \code{list_data},
#'     based on the enabled parameters: (i) exploratory plots (heatmap and/or
#'     PCA), (ii) targeted plot (volcano), and (iii) size-corrected
#'     \code{data}.
#' @export
#'
#' @examples
#' data("default_tTEscanR_tRNA_data", package = "tTEscanR")
#' data("default_tTEscanR_metadata", package = "tTEscanR")
#' DE_analysis <- runDEAnalysis(
#'     list_data = list(tRNA = default_tTEscanR_tRNA_data),
#'     metadata = default_tTEscanR_metadata, batch = "tissue",
#'     color_factor = "tissue", compute_pairwise = FALSE
#' )
runDEAnalysis <- function(list_data, metadata, batch = NULL,
    reference = NULL, reduce = 100, dim_reduct = NULL, color_factor = batch,
    heatmap = TRUE, shape_factor = NULL, label_factor = NULL, target = NULL,
    highlight_median = FALSE, numPC = 2, color_palette = NULL, fc_threshold = 1,
    show_legend = "none", padj_threshold = 0.05, label_significant = TRUE,
    compute_pairwise = TRUE, verbose = TRUE) {
    if (verbose) message("--- Running the differential expression analysis ---")
    ## Compute results
    results <- computeDEResults(
        list_data = list_data, metadata = metadata, target = target,
        batch = batch, reference = reference, reduce = reduce,
        verbose = verbose, compute_pairwise = compute_pairwise
    )

    ## Plot all datasets
    plots <- list()
    for (dataset_name in names(results)) {
        plots[[dataset_name]] <- plotDEResults(
            DE_results_list = results, dataset_name = dataset_name,
            dim_reduct = dim_reduct, color_factor = color_factor,
            shape_factor = shape_factor, label_factor = label_factor,
            highlight_median = highlight_median, numPC = numPC,
            color_palette = color_palette, target = target,
            fc_threshold = fc_threshold, padj_threshold = padj_threshold,
            heatmap = heatmap, label_significant = label_significant,
            verbose = verbose, show_legend = show_legend
        )
    }
    if (verbose) {
        message(
            "--- The differential expression analysis was completed ",
            "successfully ---"
        )
    }
    return(list(results = results, plots = plots))
}

#' Compute the DESeq2 Analysis
#'
#' @param list_data A \code{list} of matrices with features (i.e. genes, codons,
#'     anticodons or amino acids) as rows and samples or conditions as columns.
#'     The list can be named or unnamed.
#' @param metadata A \code{data.frame} with the data associated to the
#'     conditions in the matrices of \code{list_data}. There has to be one
#'     column with the same labels as the column names.
#' @param target Optional; a factor based on \code{metadata} columns to
#'     define the comparisons to perform in a targeted analysis. Defaults to
#'     \code{NULL}.
#' @param batch Optional; name of the categorical variable in \code{metadata}
#'     to correct the data. Required if \code{dim_reduct} specified.
#' @param reference Optional; factor from the \code{batch} variable to use as
#'     reference for the corrections. If not specified, the 1st factor that
#'     appears will be used instead.
#' @param compute_pairwise Logical; if \code{TRUE}, computes all the pairwise
#'     comparisons based on the conditions included in the input data.
#' @param padj_threshold Numeric; p-value threshold used for highlighting
#'     significant features in the volcano plot. Defaults to 0.05.
#' @param reduce Numeric; a scaling factor used to normalize large expression
#'     values that exceed R's handling capacity. Defaults to 100.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#'     Defaults to \code{TRUE}.
#'
#' @returns A DESeq2 object with the normalized and vst counts.
#' @export
#'
#' @examples
#' data("default_tTEscanR_tRNA_data", package = "tTEscanR")
#' data("default_tTEscanR_metadata", package = "tTEscanR")
#'
#' DE_analysis <- computeDEResults(
#'     list_data = list(tRNA = default_tTEscanR_tRNA_data),
#'     metadata = default_tTEscanR_metadata, batch = "tissue",
#'     compute_pairwise = FALSE
#' )
computeDEResults <- function(list_data, metadata, target = NULL, batch = NULL,
    reference = NULL, reduce = 100, padj_threshold = 0.05, verbose = TRUE,
    compute_pairwise = TRUE) {
    DE_results_list <- list()
    if (is.null(target) && is.null(batch)) { # Set the correction factor
        stop("Please specify the 'target' or the 'batch'.")
    }
    cf <- if (!is.null(batch)) batch else target
    ## Standardize dataset names
    ds_names <- names(list_data)
    if (is.null(ds_names)) { # Check names of the datasets in list_data
        names(list_data) <- paste0("dataset_", seq_along(list_data))
    } else {
        empty <- ds_names == ""
        if (any(empty)) {
            names(list_data)[empty] <- paste0("dataset_", which(empty))
        }
    }
    ## Pre-allocate output list vector for efficiency
    DE_results_list <- vector("list", length(list_data))
    names(DE_results_list) <- names(list_data)
    for (i in seq_along(list_data)) {
        dataset_name <- names(list_data)[i]
        if (verbose) message("- Processisng dataset: ", dataset_name)
        filter <- filterByMetadata(
            data = list_data[[i]], metadata = metadata, verbose = verbose
        )
        DESeq2_run <- computeDESeq2(
            data = filter$data, condition = target, reduce = reduce, batch = cf,
            metadata = filter$metadata, reference = reference, verbose = verbose
        )
        pairwise <- if (isTRUE(compute_pairwise)) {
            computeAllPairwiseComp( # Run all pairwise comparisons
                dds = DESeq2_run, factor_name = cf,
                padj_threshold = padj_threshold, verbose = verbose
            )
        } else {
            NULL
        }
        size_corrected <- DESeq2::counts(DESeq2_run, normalized = TRUE)
        vst_mat <- SummarizedExperiment::assay(
            DESeq2::varianceStabilizingTransformation(DESeq2_run)
        )
        DE_results_list[[dataset_name]] <- list(
            DESeq2_run = DESeq2_run, pairwise_results = pairwise, vst = vst_mat,
            size_corrected = size_corrected, metadata = filter$metadata
        ) # Assign to pre-allocated list index
    }
    return(DE_results_list)
}

#' Generates Visualizations from the DEA data
#'
#' @param DE_results_list A DESeq2 object with the normalized and vst counts.
#'     Can be obtained by running \code{\link{computeDEResults}}
#' @param dataset_name String to specify the assay in \code{DE_results_list}
#'     to extract.
#' @param dim_reduct Either \code{"PCA"}, \code{"UMAP"} or \code{"tSNE"} to
#'     specify the dimensionality reduction approach to be executed. Defaults
#'     to \code{"PCA"}.
#' @param color_factor Optional; name of the categorical variable in
#'     \code{metadata} to group by color the data points. Used if
#'     \code{dim_reduct} specified.
#' @param shape_factor Optional; name of the categorical variable to label the
#'     data points. Used if \code{dim_reduct} specified.
#' @param label_factor Optional; name of the categorical variable to label the
#'     data points. Used if \code{dim_reduct} specified.
#' @param highlight_median  Logical; if \code{TRUE} the data points of each
#'     cluster will be summarized into the median. Defaults to \code{FALSE}.
#' @param numPC Numeric; number of principal components to include in the PCA
#'     analysis. Required if \code{dim_reduct} is \code{"PCA"}. Defaults to 2.
#' @param scale_pca Logical; if \code{TRUE}, scales the data if
#'     \code{dim_reduct} is \code{"PCA"}. Defaults to \code{TRUE}.
#' @param heatmap Logical; if \code{TRUE}, generates a heatmap for exploratory
#'     analysis. Defaults to \code{TRUE}.
#' @param color_palette Optional; a \code{vector} of color codes to customize
#'     the plot appearance.
#' @param show_legend Either \code{"none"} (default), \code{"top"},
#'     \code{"bottom"}, \code{"right"} or \code{"left"} to specify the
#'     position of the legend in the plot.
#' @param target Optional; a factor based on \code{metadata} columns to
#'     define the comparisons to perform in a targeted analysis. Defaults to
#'     \code{NULL}.
#' @param fc_threshold Numeric; fold change threshold used for highlighting
#'     significant features in the volcano plot (if \code{targets} is
#'     specified). Defaults to 1.
#' @param padj_threshold Numeric; p-value threshold used for highlighting
#'     significant features in the volcano plot. Defaults to 0.05.
#' @param label_significant Logical; if \code{TRUE} displays the axis of the
#'     plots based on \code{fc_threshold} and \code{padj_threshold}. Defaults
#'     to \code{TRUE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#'     Defaults to \code{TRUE}.
#'
#' @returns Visualization of the DEA results.
#' @export
#'
#' @examples
#' data("default_tTEscanR_tRNA_data", package = "tTEscanR")
#' data("default_tTEscanR_metadata", package = "tTEscanR")
#'
#' DE_analysis <- computeDEResults(
#'     list_data = list(tRNA = default_tTEscanR_tRNA_data),
#'     metadata = default_tTEscanR_metadata, batch = "tissue",
#'     compute_pairwise = FALSE
#' )
#' DE_plots <- plotDEResults(
#'     DE_results_list = DE_analysis, dataset_name = "tRNA",
#'     dim_reduct = "PCA", color_factor = "tissue", heatmap = FALSE
#' )
plotDEResults <- function(DE_results_list, dataset_name = NULL, heatmap = TRUE,
    dim_reduct = NULL, numPC = 2, target = NULL, verbose = TRUE,
    color_factor = NULL, shape_factor = NULL, label_factor = NULL,
    highlight_median = FALSE, scale_pca = TRUE, color_palette = NULL,
    fc_threshold = 1, padj_threshold = 0.05, label_significant = TRUE,
    show_legend = "none") {
    if (is.null(dataset_name)) dataset_name <- names(DE_results_list)[1]
    if (!dataset_name %in% names(DE_results_list)) {
        stop("dataset_name '", dataset_name, "' not in 'DE_results_list'.")
    }
    vst <- DE_results_list[[dataset_name]]$vst # Extract precomputed data
    results_list <- list()
    if (isTRUE(heatmap)) {
        if (verbose) message("- Generating heatmap with Euclidean distances.")
        heatmap_plot <- produceHeatmapDiffExp(data = vst)
        results_list[["heatmap"]] <- heatmap_plot
    }
    if (!is.null(dim_reduct)) { # Exploratory Approach
        dim_reduct <- match.arg(dim_reduct, choices = c("PCA", "UMAP", "tSNE"))
        if (dim_reduct != "PCA") {
            numPC <- NULL
        } else if (is.null(numPC) || numPC < 2) {
            warning("- Invalid 'numPC' given, using default 'numPC = 2'.")
            numPC <- 2
        }
        meta <- DE_results_list[[dataset_name]]$metadata
        CheckSpecificParameters(
            color = color_factor, meta = meta,
            shape = shape_factor, label = label_factor
        )
        if (verbose) message("- Generating dimensionality reduction plots.")
        exploratory_plots <- runDimReduct(
            metadata = meta, vst = vst, median = highlight_median,
            show_legend = show_legend, dim_reduct = dim_reduct, numPC = numPC,
            color = color_factor, shape = shape_factor, label = label_factor,
            palette = color_palette, scale = scale_pca, verbose = verbose
        ) # Scale set to FALSE if the vst was already normalized
        results_list[["exploratory"]] <- exploratory_plots
    }
    if (!is.null(target)) { # Targeted Approach
        run <- DE_results_list[[dataset_name]]$DESeq2_run
        targeted_results <- targetedApproach(
            DESeq2 = run, sig = label_significant, condition = target,
            fc = fc_threshold, padj = padj_threshold, verbose = verbose
        )
        results_list[["targeted"]] <- targeted_results
    }
    return(results_list)
}

CheckSpecificParameters <- function(color, shape, label, meta) {
    if (is.null(color)) stop("Please specify a suitable 'color_factor'.")
    if (length(color) > 1 || length(shape) > 1 || length(label) > 1) {
        stop(
            "The parameters 'color_factor', 'shape_factor' and 'label_factor'",
            " are expected to be the name of a column in the 'metadata'."
        )
    }
    if (any(c(length(color), length(shape), length(label)) > 1)) {
        stop(
            "The parameters 'color_factor', 'shape_factor' and 'label_factor'",
            " must be a single column name, not a vector"
        )
    }

    provided_factors <- c(color, shape, label)

    ## If a factor is set to NULL it will not be added to the vector
    missing_factors <- setdiff(provided_factors, colnames(meta))
    if (length(missing_factors) > 0) {
        stop(
            "The following factors are not found in the metadata: ",
            missing_factors,
            collapse = ","
        )
    }
}
