#' Compute the Theoretical Translation Efficiency (tTE) score
#' @description
#' This function calculates the theoretical translation efficiency (tTE) score by integrating codon-anticodon usage and/or amino acid demand-supply across conditions.
#'
#' @param object A \code{tTEscanR_Object} containing a mRNA and a codon usage assay and/or a tRNA and an anticodon usage assay.
#' @param level Either \code{"codon"}, \code{"aa"} (default) or \code{"both"} to indicate which analysis to perform. Defaults to \code{"aa"}.
#' @param corr_method A correlation method accepted by \code{\link{cor}}. Defaults to \code{"spearman"}.
#' @param compute.significance Logical; if \code{TRUE}, computes the statistical significance (p-value) of the tTE scores. Defaults to \code{TRUE}.
#' @param overwrite Logical; if \code{TRUE}, overwrites any existing tTE table in the \code{tTEscanR_Object}. Defaults to \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return An updated \code{tTEscanR_Object} containing a new layer of information representing the translation efficiency table for the matching conditions in the mRNA and tRNA data.
#' @export
#'
#' @examples
#' data(subset_mRNA_data, subset_tRNA_data, metadata)
#' tTEscanR_obj <- Create_tTEscanR_Object(counts = list(mRNA = subset_mRNA_data,
#'                                                      tRNA = subset_tRNA_data),
#'                                        meta.data = list(metadata, "tissue"),
#'                                        meta.data.ids = list("ConditionsLabels", "CorrectionFactor"))
#' tTEscanR_obj <- ComputeCodonUsage(object = tTEscanR_obj, species = "hg38",
#'                                   additional.metrics = FALSE, reduce = 10000)
#' tTEscanR_obj <- ComputeAnticodonUsage(object = tTEscanR_obj)
#' tTEscanR_obj <- Compute_tTE(object = tTEscanR_obj, level = "codon", compute.significance = FALSE)

Compute_tTE <- function(object, level = "aa", corr_method = "spearman", compute.significance = TRUE, overwrite = FALSE, verbose = TRUE){

  ###
  # CALL: User or Run_tTEscanR_pipeline()
  # DESCRIPTION: This function computes the theoretical translation efficiency (tTE) score by correlating codon-anticodon usage and/or amino acid demand-supply across conditions.
  ###

  # This table links each `level` parameter to the different layers of data that will need to be retrieved
  assay_map_TE <- list(codon = list(mRNA_data = list("CodonUsage"), tRNA_data = list("AnticodonUsage"), level_list = list("codon")),
                       aa = list(mRNA_data = list("AADemand"), tRNA_data = list("AASupply"), level_list = list("aa")),
                       both = list(mRNA_data = list("CodonUsage", "AADemand"), tRNA_data = list("AnticodonUsage", "AASupply"), level_list = list("codon", "aa")))

  message("1 . Checking the format of the input data.")
  if (!(inherits(object, 'tTEscanR_Object'))) stop("`object` must be a tTEscanR object.")
  if (verbose) message("- The input contains a proper tTEscanR object.")
  IsIn_tTEscanR_Object(object = object, slot = "meta.data", section = "ConditionsLabels", verbose = FALSE)
  IsIn_tTEscanR_Object(object = object, slot = "meta.data", section = "CorrectionFactor")
  if (!(object@meta.data$CorrectionFactor %in% colnames(object@meta.data$ConditionsLabels))) stop("The correction factor was not found in the metadata.")
  if (verbose) message("- The ConditionsLabels and CorrectionFactor metadata have been properly loaded.")
  if (!(corr_method %in% c("spearman", "pearson", "kendall"))) stop("Please specify a suitable `corr_method`.\n", "For further details check the cor() documentation.")
  if (verbose) message("- The corr_method has been porperly specified.")
  if (!(level %in% names(assay_map_TE)) || is.null(level)) stop("Please specify a valid `level` input parameter: codon, aa, both.")
  if (verbose) message("- The level has been porperly specified.")

  level_list <- assay_map_TE[[level]]$level_list
  factor <- data_from_mRNA <- data_from_tRNA <- list() # Empty lists to be filled with the corresponding data sets

  for (i in 1:length(assay_map_TE[[level]]$level_list)){ # Iterates for the amount of assays that need to be retrieved
    # Check if the required assays are present in the tTEscanR object
    IsIn_tTEscanR_Object(object = object, slot = "assays", section = assay_map_TE[[level]]$mRNA_data[[i]], verbose = verbose)
    IsIn_tTEscanR_Object(object = object, slot = "assays", section = assay_map_TE[[level]]$tRNA_data[[i]],verbose = verbose)
    # Check that the data is in a suitable format
    CheckDataFrame(data = object@assays[[assay_map_TE[[level]]$mRNA_data[[i]]]])
    CheckDataFrame(data = object@assays[[assay_map_TE[[level]]$tRNA_data[[i]]]])
    # Store the data into the previously defined lists
    data_from_mRNA <- append(data_from_mRNA, list(object@assays[[assay_map_TE[[level]]$mRNA_data[[i]]]]))
    data_from_tRNA <- append(data_from_tRNA, list(object@assays[[assay_map_TE[[level]]$tRNA_data[[i]]]]))

    factor <- if (assay_map_TE[[level]]$level_list[[i]] == "codon") append(factor, sense_codons) else  append(factor, amino_acids)
  }

  message("  1 . COMPLETED\n", "2 . Computing the theoretical translation efficiency (tTE) analysis.")
  for (i in 1:length(data_from_mRNA)){ # Iterate over the dataset
    if (verbose) message(paste("Level:", level_list[[i]], "\n"), "\n------------------------------\n", "A . Assessing the data and metadata.")

    # Retain matching conditions in the mRNA and tRNA data
    filter_results <- FilterMatrix(data_mRNA = data_from_mRNA[[i]], data_tRNA = data_from_tRNA[[i]], level = level_list[[i]], verbose = verbose)

    # Filter the data and metadata based on the shared content - The generated matrices are already checked by the helper function
    if (verbose) message("- Filtering the metadata by the intersecting conditions.")
    filtered_mRNA <- FilterByMetadata(data = filter_results[[1]], metadata = object@meta.data$ConditionsLabels, verbose = FALSE) # Evaluating the mRNA data
    mRNA_filtered <- filtered_mRNA[[1]]
    metadata_filtered <- filtered_mRNA[[2]] # The metadata is just retrieved once as it has been previously checked the consistency between mRNA and tRNA datasets
    filtered_tRNA <- FilterByMetadata(data = filter_results[[2]], metadata = object@meta.data$ConditionsLabels, verbose = FALSE) # Evaluating the tRNA data
    tRNA_filtered <- filtered_tRNA[[1]]

    # Compute the size-correction of the data
    mRNA_filtered <- ComputeSizeCorrection(data = mRNA_filtered, metadata = metadata_filtered, corr_factor = object@meta.data$CorrectionFactor, verbose = verbose)
    tRNA_filtered <- ComputeSizeCorrection(data = tRNA_filtered, metadata = metadata_filtered, corr_factor = object@meta.data$CorrectionFactor, verbose = verbose)

    # Compute the correlation between the mRNA and the tRNA data
    if (verbose) message("B . Computing the tTE scores.")
    tTE_results_tibble <- ComputeCorrelation(data_mRNA = mRNA_filtered, data_tRNA = tRNA_filtered, corr_method = corr_method, verbose = verbose)

    if (isTRUE(compute.significance)){
      if (verbose) message("C . Computing the statistical significance of the tTE scores.")
      IsIn_tTEscanR_Object(object = object, slot = "assays", section = "tRNA", update.assay = FALSE, overwrite = FALSE, verbose = FALSE) # Check if the tRNA expression matrix is present in the tTEscanR object

      # Filter the tRNA expression data based on the shared conditions extracted above
      tRNA_exp_data <- object@assays$tRNA[, colnames(filter_results[[1]])] # Usage of tRNA expression matrix to increment the number of features considered in the shuffling analysis
      CheckDataFrame(data = tRNA_exp_data)

      # Compute the size-correction
      tRNA_exp_data <- suppressMessages(ComputeSizeCorrection(data = tRNA_exp_data, metadata = metadata_filtered, corr_factor = object@meta.data$CorrectionFactor, verbose = FALSE))

      tTE_results_tibble <- ComputeStatisticalSignificance(level = level_list[[i]], tTE_scores = tTE_results_tibble, data_mRNA = mRNA_filtered,
                                                           data_tRNA = tRNA_filtered, tRNA_exp = tRNA_exp_data, corr_method = corr_method, verbose = verbose)
    }

    ids_results <- if (level_list[[i]] == "aa") "tTEresults_AA" else "tTEresults_codon" # Define the labels to identify the outputs

    if (verbose) message("- Updating the tTEscanR object.")
    object <- suppressMessages(Update_tTEscanR_Object(object = object, meta.data = list(tTE_results_tibble), meta.data.ids = list(ids_results), overwrite.metadata = overwrite, verbose = FALSE))
    if (verbose) message(paste("  Level:", level_list[[i]], "\n"), "COMPLETED\n", "------------------------------\n")
  }
  message("  2 . COMPLETED")
  return(object) # The output tTEscanR object has been validated in Update_tTEscanR_Object()
}
