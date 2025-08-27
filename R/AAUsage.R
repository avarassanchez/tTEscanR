#' Compute Amino Acid (AA) Demand or Supply
#' @description
#' This function calculates the amino acid (AA) demand and/or supply from a codon and/or anticodon usage matrices of a \code{tTEscanR_Object}.
#' It aggregates the contribution of their features based on the standard genetic code (mapping codons/anticodons to AA),
#' The resulting values reflect total usage (demand) or availability (supply) of each amino acid, depending on the input type.
#'
#' @param object A \code{tTEscanR_Object} containing codon and/or anticodon usage assays to be analyzed.
#' @param level Either \code{"demand"}, \code{"supply"} or \code{"both"} to indicate which analysis to perform.
#' @param overwrite.assay Logical; if \code{TRUE}, overwrites any existing assay in the \code{object}. Defaults to \code{FALSE}
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}
#'
#' @return An updated \code{tTEscanR_Object} containing a new layer of information in the \code{assays} slot representing the AA demand and/or supply.
#' @export
#'
#' @examples
#' data(subset_tRNA_data)
#' tTEscanR_obj <- Create_tTEscanR_Object(counts = subset_tRNA_data, assay = "tRNA")
#' tTEscanR_obj <- ComputeAnticodonUsage(object = tTEscanR_obj)
#' tTEscanR_obj <- ComputeAAUsage(object = tTEscanR_obj, level = "supply")

ComputeAAUsage <- function(object, level, overwrite.assay = FALSE, verbose = TRUE){

  ###
  # CALL: User or Compute_tTE()
  # DESCRIPTION: This function computes separately or simultaneously the AA demand/supply of codon/anticodon usage matrices.
  ###

  # This table links each `level` parameter to the different layers of data that will need to be retrieved
  assay_map_AA <- list(demand = list(assay_id = list("AADemand"), check_assays = list("CodonUsage"), AAfunction = "sense"),
                       supply = list(assay_id = list("AASupply"), check_assays = list("AnticodonUsage"), AAfunction = "reverse"),
                       both = list(assay_id = list("AADemand", "AASupply"), check_assays = list("CodonUsage", "AnticodonUsage"), AAfunction = list("sense", "reverse")))

  message("1 . Checking the format of the input data.")
  if (!(inherits(object, 'tTEscanR_Object'))) stop("`object` must be a tTEscanR object.")
  if (verbose) message("- The input contains a proper tTEscanR object.")
  if (!level %in% names(assay_map_AA) || is.null(level)) stop("Please specify a valid `level` input parameter: demand, supply, both.")
  if (verbose) message("- The level has been porperly specified.")

  amino_acids_function <- data <- assay_id <- list() # Create empty lists to store the outputs

  for (i in 1:length(assay_map_AA[[level]]$check_assays)){ # Iterates for the amount of assays that need to be retrieved
    # Check if the required assays are present in the tTEscanR object
    retrieve_data <- RetrieveAAUsageData(object = object, data_section = assay_map_AA[[level]]$check_assays[[i]],
                                         assay_id = assay_map_AA[[level]]$assay_id[[i]], data_function = assay_map_AA[[level]]$AAfunction[[i]])
    data <- append(data, retrieve_data[[1]])
    assay_id <- append(assay_id, retrieve_data[[2]])
    amino_acids_function <- append(amino_acids_function, retrieve_data[[3]])
  }

  if (verbose) message("- The mRNA and/or tRNA data matrices have been properly loaded.")
  message("  1 . COMPLETED\n", "2 . Pooling counts from each amino acid.")

  for(i in 1:length(data)){ # Iterate over each assay stored in data
    if (verbose) message(paste("- Performing the", assay_id[[i]], "analysis.")) # Reports the name of the assay analyzed in each iteration
    AAmatrix <- GroupAA(data = data[[i]], amino_acids_function = amino_acids_function[[i]]) # Group the codons/anticodons into AA
    object <- suppressMessages(Update_tTEscanR_Object(object = object, counts = AAmatrix, assay = assay_id[[i]], overwrite.assay = overwrite.assay, verbose = FALSE))
  }
  message("  2 . COMPLETED\n", "3 . Retrieving the tTEscanR object.\n", "  3 . COMPLETED")
  return(object) # The output tTEscanR object has been validated in Update_tTEscanR_Object()
}
