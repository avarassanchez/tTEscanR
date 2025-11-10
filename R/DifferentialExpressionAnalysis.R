#' Compute Size Correction to Account for Sequencing Depth
#' @description
#' This function normalizes the data to correct for differences in sequencing depth across samples, ensuring comparability.
#' It is particularly useful prior to comparative analysis across conditions where bias may be introduced.
#'
#' @param data A \code{matrix} with features (i.e. genes, codons, anticodons or amino acids) as rows and samples or conditions as columns.
#' @param metadata A \code{data.frame} with the meta-information related with the conditions in \code{data}. There has to be one column with the same labels as the column names.
#' @param corr_factor A factor based on \code{metadata} columns to define the variable to correct for. Required if \code{runDESeq2} is \code{TRUE}.
#' @param reduce Numeric; a scaling factor used to normalize large expression values that exceed R's handling capacity. Defaults to 100.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return A size corrected \code{data} matrix.
#' @export
#'
#' @examples
#' data(subset_mRNA_data, metadata)
#' size_correction <- ComputeSizeCorrection(data = subset_mRNA_data, metadata = metadata,
#'                                          corr_factor = "tissue")

ComputeSizeCorrection <- function(data, metadata, corr_factor = NULL, reduce = 100, verbose = TRUE){

  ###
  # CALL: User and multiple functions inside tTEscanR
  # DESCRIPTION: This function uses the library DESeq2 to size correct the data to account for sequencing depth.
  ###

  message("1 . Checking the format of the input data.")
  CheckDataFrame(data = data) # Evaluate the input parameter: `data`
  if (verbose) message("- The data parameter has been properly loaded.")
  CheckDataFrame(data = metadata) # Evaluate the input parameters: `metadata`and `corr_factor`
  if (!(corr_factor %in% colnames(metadata))) stop("The correction factor was not found in the metadata.")
  if (verbose) message("- The metadata and corr_factor parameters have been properly loaded.")

  message("  1 . COMPLETED\n", "2 . Executing DESeq2.")
  if (verbose) message("- Filtering (if necessary) the `data` and `metadata` for common features.")
  filtered <- FilterByMetadata(data = data, metadata = metadata, verbose = verbose) # Check consistency in conditions included in the data and metadata
  DESeq2_run <- suppressWarnings(ComputeDESeq2(data = filtered[[1]], metadata = filtered[[2]], corr_factor = corr_factor, reduce = reduce))

  message("  2 . COMPLETED\n", "3 . Extracting the size-corrected counts.")
  size_corrected_output_matrix <- suppressMessages(DESeq2::counts(DESeq2_run, normalized = TRUE)) # Extracting the normalized data from the DESeq2 object
  CheckDataFrame(data = size_corrected_output_matrix)
  message("  3 . COMPLETED")
  return(size_corrected_output_matrix) # Return the size-corrected matrix
}

#' Perform Differential Expression Analysis Using DESeq2
#' @description
#' This function applies differential expression analysis using the DESeq2 framework on a matrix of expression values.
#' It supports both exploratory visualizations (heatmap and PCA) and targeted comparisons using a custom contrast table.
#'
#' @param list_data A \code{list} of matrices with features (i.e. genes, codons, anticodons or amino acids) as rows and samples or conditions as columns. The list can be named or unnamed.
#' @param metadata A \code{data.frame} with the data associated to the conditions in the matrices of \code{list_data}. There has to be one column with the same labels as the column names.
#' @param target_factor Optional; a factor based on \code{metadata} columns to define the comparisons to perform in a targeted analysis. Defaults to \code{NULL}.
#' @param fc_threshold Numeric; fold change threshold used for highlighting significant features in the volcano plot (if \code{targets} is specified). Defaults to 1.
#' @param pval_threshold Numeric; p-value threshold used for highlighting significant features in the volcano plot. Required if \code{targets} is specified. Defaults to 0.05.
#' @param sig_axis Logical; if \code{TRUE} displays the axis of the plots based on \code{fc_threshold} and \code{pval_threshold}. Defaults to \code{TRUE}.
#' @param reduce Numeric; a scaling factor used to normalize large expression values that exceed R's handling capacity. Defaults to 100.
#' @param heatmap Logical; if \code{TRUE}, generates a heatmap for exploratory analysis. Not applicable if \code{targets} is specified. Defaults to \code{TRUE}.
#' @param dim.reduct Either \code{"PCA"}, \code{"UMAP"} or \code{"tSNE"} to specify the dimensionality reduction approach to be executed. Defaults to \code{"PCA"}.
#' @param numPC Numeric; number of principal components to include in the PCA analysis. Required if \code{dim.reduct} is \code{"PCA"}. Defaults to 2.
#' @param corr_factor Optional; name of the categorical variable to correct the data. Required if \code{dim.reduct} specified.
#' @param color_factor Optional; name of the categorical variable to group by color the data points. Used if \code{dim.reduct} specified.
#' @param shape_factor Optional; name of the categorical variable to group by shape the data points. Used if \code{dim.reduct} specified.
#' @param color_palette Optional; a \code{vector} of color codes to customize the plot appearance.
#' @param label_factor Optional; name of the categorical variable to label the data points. Used if \code{dim.reduct} specified.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return A \code{list} of outputs per each matrix in \code{list_data}, based on the enabled parameters: (i) exploratory plots (heatmap and/or PCA), (ii) targeted plot (volcano), and (iii) size-corrected \code{data}.
#' @export
#'
#' @examples
#' data(subset_mRNA_data, metadata)
#' DE_analysis <- RunDEAnalysis(list_data = list(mRNA = subset_mRNA_data),
#'                              metadata = metadata, corr_factor = "tissue")

RunDEAnalysis <- function(list_data, metadata, target_factor = NULL, corr_factor = NULL, fc_threshold = 1, pval_threshold = 0.05, reduce = 100, heatmap = TRUE,
                          dim.reduct = "PCA", numPC = 2, color_factor = NULL, shape_factor = NULL, color_palette = NULL, label_factor = NULL, sig_axis = TRUE, verbose = TRUE){

  ###
  # CALL: User and Run_tTEscanR_pipeline()
  # DESCRIPTION: This function executes the differential expression analysis.
  ###

  if (verbose) message("1 . Checking the input parameters.")
  if (is.null(names(list_data)) || all(names(list_data) == "")) names(list_data) <- paste("list_data", seq_along(list_data)) # If list_data is an unnamed list the name is assigned based on appearance order
  if (!is.null(dim.reduct) && (is.null(corr_factor) || !(corr_factor %in% colnames(metadata)))) stop("The correction factor was not found in the metadata.") # Checking the variable that will be used to correct the data when running DESeq2
  CheckDataFrame(data = metadata)

  if (!is.null(dim.reduct)) {
    if (!dim.reduct %in% c("PCA", "UMAP", "tSNE")) stop("Please indicate a suitable `dim.reduct` value.\n", "Supported formats: PCA, UMAP or tSNE")
  }

  # Checking the variable that will be used to classify the entries in a targeted approach
  if (!(is.null(target_factor)) && !(target_factor %in% colnames(metadata))) stop("Please indicate a suitable `target_factor`.\n", paste("Current `target_factor` =", target_factor, "not in `metadata`.\n"),
                                                                          paste("Consider as potential `target_factor`:", paste(colnames(metadata), collapse = ",")))
  if(!(is.null(label_factor))){
    if (!(label_factor %in% colnames(metadata))){
      message("The label_factor was not found in the metadata.\n", "The plot(s) will not display any label.")
      labels <- FALSE
    } else {
      labels <- TRUE
    }
  } else {
    labels <- FALSE
  }

  message("  1 . COMPLETED\n", "2 . Running the differential expression analysis.")
  DESeq2_results_list <- list() # Create an empty list to store the results

  for (i in 1:length(list_data)){ # Iterates over the assays included in list_data
    message(paste("- Analyzing:", names(list_data)[i]))

    if (verbose) message("- Filtering (if necessary) the `data` and `metadata` for common features.")
    filtered <- FilterByMetadata(data = list_data[[i]], metadata = metadata, verbose = verbose)

    ####
    if (verbose) message("- Running differential expression analysis.")
    if (!is.null(target_factor)) corr_factor <- target_factor
    DESeq2_run <- suppressWarnings(ComputeDESeq2(data = filtered[[1]], metadata = filtered[[2]], corr_factor = corr_factor, reduce = reduce))

    if (verbose) message("     Computing size correction to account for sequencing depth.")
    size_corrected_output_matrix <- suppressMessages(DESeq2::counts(DESeq2_run, normalized = TRUE))
    CheckDataFrame(data = size_corrected_output_matrix)

    if (verbose) message("     Computing variance stabilization.")
    vst <- DESeq2::varianceStabilizingTransformation(DESeq2_run)
    vst <- SummarizedExperiment::assay(vst)
    CheckDataFrame(data = vst)
    if (verbose) message("- Differential expression analysis completed.")
    ####

    if (is.null(target_factor)){ # EXPLORATORY APPROACH
      if (is.null(color_factor)){
        color_factor <- corr_factor
        message("No `color`_factor` given, the `corr_factor` will be used instead.")
      }
      results.list <- ExploratoryApproach(metadata = filtered[[2]], vst = vst, heatmap = heatmap, dim.reduct = dim.reduct, color_factor = color_factor,
                                          numPC = numPC, shape_factor = shape_factor, color_palette = color_palette, label_factor = label_factor, verbose = verbose)
    } else { # TARGETED APPROACH
      if (verbose) message("- Executing a targeted analysis.")
      if (length(unique(filtered[[2]][[target_factor]])) < 2) stop("A single class was identified.\n", paste("Current `target_factor` = ", target_factor, ".\n"),
                                                                   paste("Class distribution:", as.character(length(unique(filtered[[2]][[target_factor]])))))

      results.list <- TargetedApproach(DESeq2_run = DESeq2_run, target_factor = target_factor,verbose = verbose,
                                       fc_threshold = fc_threshold, pval_threshold = pval_threshold, sig_axis = sig_axis)
    }

    results.list[["size_corrected"]] <- size_corrected_output_matrix
    results.list[["vst"]] <- vst
    DESeq2_results_list[[names(list_data)[i]]] <- results.list
    message(paste("- COMPLETED analyzing:", names(list_data)[i]))
  }
  message("  2 . COMPLETED\n")
  return (DESeq2_results_list) # The output consists of a nested list, which elements depend on the parameters enabled.
}
