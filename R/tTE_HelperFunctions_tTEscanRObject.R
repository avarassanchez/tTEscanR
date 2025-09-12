#' tTEscanR Object Update
#' @description
#' Updates an already created \code{tTEscanR_object} using \code{\link{Create_tTEscanR_Object}}. For more details, refer to \code{\link{Create_tTEscanR_Object}}.
#'
#' @param object An existing \code{tTEscanR_Object}.
#' @param counts Optional; a count matrix (or \code{list} of matrices) that will be stored in the \code{assays} slot of the input \code{object}. Supported formats: \code{matrix}, \code{data.frame} and \code{list}.
#' @param assay Optional; a character string (or vector of strings) specifying the name of the \code{counts}. Supported formats: \code{character} and \code{list}.
#' @param main_name Optional; a character string specifying the name of the \code{meta.data} if dealing with a \code{list} that needs to be added as a single element.
#' @param meta.data Optional; a \code{list} with additional information that will be stored in \code{meta.data} slot of the input \code{object}.
#' @param meta.data.ids Optional; a \code{list} with the labels to identify the \code{meta.data} (if \code{meta.data} is given).
#' @param overwrite.assay Logical; if \code{TRUE}, overwrites any existing assay in the \code{object} if the \code{assay} label coincides. Defaults to \code{FALSE}.
#' @param overwrite.metadata Logical; if \code{TRUE}, overwrites any existing metadata in the \code{object} if the \code{meta.data.ids} label coincides. Defaults to \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return An updated \code{tTEscanR_Object}.
#' @export
#'
#' @examples
#' data(subset_mRNA_data, subset_tRNA_data)
#' tTEscanR_obj <- Create_tTEscanR_Object(counts = subset_mRNA_data, assay = "mRNA")
#' tTEscanR_obj <- Update_tTEscanR_Object(object = tTEscanR_obj, counts = subset_tRNA_data,
#'                                        assay = "tRNA")

Update_tTEscanR_Object <- function(object, counts = NULL, assay = NULL, main_name = NULL, meta.data = NULL, meta.data.ids = NULL,
                                   overwrite.assay = FALSE, overwrite.metadata = FALSE, verbose = TRUE){

  ###
  # CALL: User or Create_tTEscanR_Object()
  # DESCRIPTION: This function adds new data and metadata to an existing tTEscanR object.
  # The code follows the same background structure as Create_tTEscanR_Object(), with the difference that there is already an object present.
  # In all cases it is required to check that the data is not already present, otherwise it is required to specify the overwrite parameters.
  ###

  message("1 . Checking the format of the input data.")
  if (!inherits(object, 'tTEscanR_Object')) stop("`object` must be a tTEscanR object.")
  if (verbose) message("- The input consists of a proper tTEscanR object.")

  # NOTHING TO ADD
  if (is.null(counts) && is.null(meta.data)) stop("Please specify either the parameter `counts` or `meta.data` to be added to the `object`.")
  message("  1 . COMPLETED")
  count <- 2

  # ADDING COUNTS
  if (!is.null(counts)){
    message(paste(as.character(count), ". Adding the `counts` to the tTEscanR object."))
    object <- DefineNewData(object = object, slot = "assays", data = counts, id = assay, mode = "fix", action.update = TRUE, overwrite = overwrite.assay, verbose = verbose)
    message(paste(" ", as.character(count), ". COMPLETED"))
    count <- count + 1
  }

  # ADDING METADATA
  if (!is.null(meta.data)){
    message(paste(as.character(count), ". Adding the `meta.data` to the tTEscanR object."))
    object <- DefineNewData(object = object, slot = "meta.data", data = meta.data, id = meta.data.ids, mode = "flexible",
                            main.name = main_name, action.update = TRUE, overwrite = overwrite.metadata, verbose = verbose)
    message(paste(" ", as.character(count), ". COMPLETED"))
    count <- count + 1
  }

  message(paste(as.character(count), ". Validating and returning the updated tTEscanR object."))
  validObject(object = object)
  message(paste(" ", as.character(count), ". COMPLETED"))
  return(object) # The output is an updated tTEscanR object
}

DefineNewData <- function(object = NULL, slot = NULL, data, id = NULL, action.update = FALSE, overwrite = FALSE, main.name = NULL, mode, verbose){

  ###
  # CALL: Create_tTEscanR_Object() and Update_tTEscanR_Object()
  # DESCRIPTION: This function controls the addition of data and metadata into a tTEscanR object.
  ###

  if (verbose) message(paste("- Evaluating the parameters", deparse(substitute(data)), "and", deparse(substitute(id)), "."))
  input_format <- IdentifyInputFormat(data = data, mode = mode)
  if (!input_format %in% c("list", "single")) stop("Internal error: IdentifyInputFormat returned an invalid mode.")

  # Check inputs and get final, validated IDs
  id <- CheckInputCombinations(data = data, labels = id, mode = input_format)

  if (input_format == "list"){
    if (isTRUE(action.update)){
      for (i in seq_along(data)){
        if (is.null(main.name)){
          IsIn_tTEscanR_Object(object = object, slot = slot, section = id[[i]], update.assay = TRUE, overwrite = overwrite, verbose = verbose)
          if (slot == "assays") object@assays[[id[[i]]]] <- data[[i]] else object@meta.data[[id[[i]]]] <- data[[i]]
        } else {
          if (is.null(object@meta.data[[main.name]])) object@meta.data[[main.name]] <- list()
          IsIn_tTEscanR_Object(object = object, slot = "meta.data", section = main.name, subset = id[[i]], update.assay = TRUE, overwrite = overwrite, verbose = verbose)
          object@meta.data[[main.name]][[id[[i]]]] <- data[[i]]
        }
      }
      return(object)
    } else {
      results.list <- data
      names(results.list) <- id
      return(results.list)
    }
  } else { # Dealing with a single input
    if (isTRUE(action.update)){
      IsIn_tTEscanR_Object(object = object, slot = slot, section = id, update.assay = TRUE, overwrite = overwrite, verbose = verbose)
      if (slot == "assays") object@assays[[id]] <- data else object@meta.data[[id]] <- data
      return(object)
    } else {
      results.list <- list(data)
      names(results.list) <- id
      return(results.list)
    }
  }
}

CheckInputCombinations <- function(data, labels, mode){

  ###
  # CALL: Create_tTEscanR_Object() and Update_tTEscanR_Object()
  # DESCRIPTION: This function checks all input combinations that are suitable to be added to a tTEscanR object to avoid future inconsistencies.
  # If dealing with lists we can have named or unnamed lists. Named lists already have an identifier for each piece of data.
  # However, unnamed lists will require an extra input parameter to specify the identifiers.
  ###

  if (mode == "list"){
    if (is.null(labels)){
      if (is.null(names(data)) || all(names(data) == "")) stop("The parameter is an unnamed list.\nPlease provide a named list or specify a suitable set of labels.")
      return(names(data))
    } else {
      if(!(is.null(names(data)))) stop("The parameter is a named list, but a labels parameter has also been input.\nPlease identify the data just in one way.")

      labels_vec <- unlist(labels)

      if (length(data) != length(labels_vec)) stop("The parameters require a 1:1 relation.")
      if (!is.character(labels_vec)) stop(paste("Element in `meta.data.ids` is not a string.\n Found class:", class(labels_vec)))
      if (anyDuplicated(labels_vec)) {
        duplicated_labels <- unique(labels_vec[duplicated(labels_vec)])
        stop("Please provide unique labels to identify the each piece of data.\n", "Duplicated labels: ", paste(duplicated_labels, collapse = ","))
      }
      return(labels_vec)
    }
  } else if (mode == "single"){
    if (is.null(labels) || !is.character(labels) || length(labels) != 1) stop("Please provide a suitable parameter to identify the `data`.\n", "Supported formats: string.")
    return(labels)
  }

  stop("Internal error: unhandled input mode.")
}

IsIn_tTEscanR_Object <- function(object, slot, section, subset = NULL, compute.assay = FALSE, update.assay = FALSE, overwrite = FALSE, verbose = TRUE){

  ###
  # CALL: Multiple
  # DESCRIPTION: This function has 2 modes: (i) compute mode - check if all the inputs required are available, (ii) update mode - check if the data will need to be overwritten.
  ###

  # Validate the slots
  if (is.null(slot) || !(slot %in% c("assays", "meta.data"))) stop("Wrong `slot` specified. Supported formats: assays and meta.data.")

  # Handle the case where the section does NOT exist
  if (!isTRUE(section %in% names(slot(object, slot)))) {
    if (isTRUE(update.assay)) {
      if (verbose) message(paste("- The", section, "section will be defined."))
      return(invisible(NULL))
    }

    if (isFALSE(compute.assay)) stop(paste("A", section, "section does not exist in the input object."))
    if (verbose) message(paste("- A", section, "section does not exist in the input object."))

    return(isFALSE(compute.assay))
  }

  # Handle the case where the section DOES exist
  if (verbose) message(paste("- A", section, "section exists in the input object."))

  # Check for nested subset if applicable
  if (!is.null(subset)) {
    subset_exists <- subset %in% names(object@meta.data[[section]])
    if (isTRUE(subset_exists)) {
      if (verbose) message(paste("- A", subset, "subset exists in the input object."))
      if (isTRUE(update.assay)) {
        if (isFALSE(overwrite)) stop(paste("Specify another name for", subset, "or change the input parameter `overwrite` to TRUE."))
        if (verbose) message(paste("- The", subset, "subset will be overwritten."))
      }
    }
  } else {
    # No subset specified, check if main section should be overwritten
    if (isTRUE(update.assay)) {
      if (isFALSE(overwrite)) stop(paste("Specify another name for", section, "or change the input parameter `overwrite` to TRUE."))
      if (verbose) message(paste("- The", section, "section will be overwritten."))
    }
  }

  return(invisible(NULL))
}
