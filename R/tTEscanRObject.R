#' @importFrom methods setClass new validObject
NULL

#' The tTEscanR Class
#' @description
#' The **tTEscanR** object is dynamically updated to store \code{assays} and \code{meta.data} at each analysis step.
#' Ensures efficient tracking and organization of inputs and outputs throughout the pipeline.
#' In order to ensure robustness throughout the pipeline, specific ids have been assigned and should be respected by the user.
#'
#' @slot assays A list of assays.
#' @slot meta.data A list of meta-information associated with the assays.
#'
#' @name tTEscanR_Object-class
#' @rdname tTEscanR_Object-class
#' @exportClass tTEscanR_Object

setClass(Class = 'tTEscanR_Object', slots = c(assays = 'list', meta.data = 'list'))

#' Create a tTEscanR Object
#' @description
#' This function initializes a \code{tTEscanR_Object}, a structured container designed to hold translational efficiency-related data.
#' The object consists of two main components; \code{assays}, which stores data matrices (e.g. expression, codon usage),
#' and \code{meta.data}, which stores associated information and additional intermediate calculations.
#' A \code{tTEscanR_Object} can be created using a single dataset, or initialized with a list of datasets provided at once.
#' Additional assays and metadata layers can be appended later using the \code{\link{Update_tTEscanR_Object}} function.
#'
#' In order to ensure robustness throughout the pipeline **specific ids** have been assigned and should be respected by the user.
#' **\code{assays} slot:**
#' - mRNA and tRNA count matrices as \code{"mRNA"} and \code{"tRNA"}
#' - Codon and anticodon usage count matrices as \code{"CodonUsage"} and \code{"AnticodonUsage"}
#' - Amino acid demand and supply count matrices as \code{"AADemand"} and \code{"AASupply"}
#' **\code{meta.data} slot:**
#' - Table with the conditions of the mRNA and tRNA data as \code{"ConditionsLabels"}
#' - Active correction factor to use when running differential expression analyses as \code{"CorrectionFactor"}
#'
#' @param counts A count \code{matrix} (or \code{list} of matrices) that will be stored in the \code{assays} slot.
#' @param assay Optional; a \code{character} string (or \code{list} of characters) to identify the \code{counts}. Required if \code{counts} is not a named list.
#' @param meta.data Optional; a variable (or \code{list} of variables) with additional information that will be stored in \code{meta.data} slot.
#' @param meta.data.ids Optional; a \code{character} string (or \code{list} with the labels) to identify the \code{meta.data}. Required if \code{meta.data} is not a named list.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return A \code{tTEscanR_Object}.
#' @export
#'
#' @examples
#' data(subset_mRNA_data, subset_tRNA_data, metadata)
#' tTEscanR_obj <- Create_tTEscanR_Object(counts = list(mRNA = subset_mRNA_data,
#'                                                      tRNA = subset_tRNA_data),
#'                                        meta.data = metadata, meta.data.ids = "ConditionsLabels")

Create_tTEscanR_Object <- function(counts, assay = NULL, meta.data = NULL, meta.data.ids = NULL, verbose = TRUE){

  ###
  # CALL: User
  # DESCRIPTION: This function creates a tTEscanR object.
  # Inputs are accepted in multiple formats.
  # More data can be added at the time of creation or afterwards using Update_tTEscanR_Object().
  # Using this function the metadata elements (if dealing with lists) are added independently.
  ###

  message("1 . Checking the format of the input data and adding it to the object.")
  if (verbose) message("- Adding the counts to a new tTEscanR object.") # Mandatory
  assay.list <- DefineNewData(data = counts, id = assay, mode = "fix", verbose = verbose)

  if (!is.null(meta.data)){
    if (verbose) message("- Adding the metadata to a new tTEscanR object.")
    metadata.list <- DefineNewData(data = meta.data, id = meta.data.ids, mode = "flexible", verbose = verbose)
    object <- suppressWarnings(expr = new(Class = 'tTEscanR_Object', assays = assay.list, meta.data = metadata.list)) # Generate the object
  } else {
    object <- suppressWarnings(expr = new(Class = 'tTEscanR_Object', assays = assay.list)) # Generate the object
  }

  message("  1 . COMPLETED\n", "2 . Validating and returning the tTEscanR object.")
  validObject(object = object)
  message("  2 . COMPLETED")
  return(object)
}
