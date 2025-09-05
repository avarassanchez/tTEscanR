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

ComputeSizeCorrection <- function(data, metadata = NULL, corr_factor = NULL, reduce = 100, verbose = TRUE){

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
  DESeq2_run <- suppressWarnings(ComputeDESeq2(data = filtered[[1]], metadata = filtered[[2]], corr_factor = corr_factor, reduce = reduce, targets = NULL, verbose = verbose))

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
#' @param targets Optional; a \code{data.frame} with one column indicating the conditions to select, and another column with the labels for the comparisons.
#' @param target_factor Optional; a factor based on \code{metadata} columns to define the variable to classify the entries in \code{targets}. Required if \code{targets} provided.
#' @param fc_threshold Numeric; fold change threshold used for highlighting significant features in the volcano plot (if \code{targets} is specified). Defaults to 1.
#' @param pval_threshold Numeric; p-value threshold used for highlighting significant features in the volcano plot. Required if \code{targets} is specified. Defaults to 0.05.
#' @param reduce Numeric; a scaling factor used to normalize large expression values that exceed R's handling capacity. Defaults to 100.
#' @param heatmap Logical; if \code{TRUE}, generates a heatmap for exploratory analysis. Not applicable if \code{targets} is specified. Defaults to \code{TRUE}.
#' @param PCA Logical; if \code{TRUE}, generates a principal component analysis. Not applicable if \code{targets} is specified. Defaults to \code{TRUE}.
#' @param numPC Numeric; number of principal components to include in the PCA analysis. Required if \code{PCA} is \code{TRUE}. Defaults to 2.
#' @param corr_factor Optional; name of the categorical variable to correct the data. Required if \code{PCA} is \code{TRUE}.
#' @param color_factor Optional; name of the categorical variable to group by color the data points. Required if \code{PCA} is \code{TRUE}.
#' @param shape_factor Optional; name of the categorical variable to group by shape the data points. Required if \code{PCA} is \code{TRUE}.
#' @param color_palette Optional; a \code{vector} of color codes to customize PCA plot appearance.
#' @param labels Logical; if \code{TRUE} includes the data points labels in the PCA plot. Required if \code{PCA} is \code{TRUE}. Defaults to \code{FALSE}.
#' @param label_factor Optional; name of the categorical variable to label the data points. Required if \code{PCA} and \code{labels} are \code{TRUE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return A \code{list} of outputs per each matrix in \code{list_data}, based on the enabled parameters: (i) heatmap, (ii) PCA plot, and (iii) size corrected \code{data}.
#' @export
#'
#' @examples
#' data(subset_mRNA_data, metadata)
#' DE_analysis <- ExecuteDESeq2runner(list_data = list(mRNA = subset_mRNA_data),
#'                                    metadata = metadata, corr_factor = "tissue")

ExecuteDESeq2runner <- function(list_data, metadata, targets = NULL, target_factor = NULL, corr_factor = NULL, fc_threshold = 1, pval_threshold = 0.05, reduce = 100, heatmap = TRUE,
                                PCA = TRUE, numPC = 2, color_factor = NULL, shape_factor = NULL, color_palette = NULL, labels = FALSE, label_factor = NULL, verbose = TRUE){

  ###
  # CALL: User and Run_tTEscanR_pipeline()
  # DESCRIPTION: This function executes the differential expression analysis.
  ###

  if (verbose) message("1 . Checking the input parameters.")
  if (is.null(names(list_data)) || all(names(list_data) == "")) names(list_data) <- paste("list_data", seq_along(list_data)) # If list_data is an unnamed list the name is assigned based on appearance order
  CheckDataFrame(data = metadata)
  if (isTRUE(PCA) && (is.null(corr_factor) || !(corr_factor %in% colnames(metadata)))) stop("The correction factor was not found in the metadata.") # Checking the variable that will be used to correct the data when running DESeq2

  # Checking the variable that will be used to classify the entries in a targeted approach
  if (!(is.null(targets)) && !(target_factor %in% colnames(metadata))) stop("Please indicate a suitable `target_factor`.\n",
                                                                          paste("Current `target_factor` =", target_factor, "not in `metadata`.\n"),
                                                                          paste("Consider as potential `target_factor`:", paste(colnames(metadata), collapse = ",")))

  if (isTRUE(labels) &&  (is.null(label_factor) || !(label_factor %in% colnames(metadata)))){
    message("No label_factor was input or it was not found in the metadata.\n", "The PCA plot(s) will not display any label.")
    labels <- FALSE
  }

  message("  1 . COMPLETED\n", "2 . Running the differential expression analysis.")
  DESeq2_results_list <- list() # Create an empty list to store the results
  for (i in 1:length(list_data)){ # Iterates over the assays included in list_data
    message(paste("- Analyzing:", names(list_data)[i]))
    DESeq2_results_list[[names(list_data)[i]]] <- DESeq2runner(data = list_data[[i]], metadata = metadata, targets = targets, target_factor = target_factor, corr_factor = corr_factor,
                                                               fc_threshold = fc_threshold, pval_threshold = pval_threshold, reduce = reduce, heatmap = heatmap, PCA = PCA, numPC = numPC,
                                                               color_factor = color_factor, shape_factor = shape_factor, color_palette = color_palette, labels = labels, label_factor = label_factor, verbose = verbose)
    message(paste("- COMPLETED analyzing:", names(list_data)[i]))
  }
  message("  2 . COMPLETED\n")
  return (DESeq2_results_list) # The output consists of a nested list, which elements depend on the parameters enabled.
}
