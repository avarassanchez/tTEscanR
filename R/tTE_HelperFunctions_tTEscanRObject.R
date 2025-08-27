DefineNewData <- function(object = NULL, slot = NULL, data, id = NULL, action.update = FALSE, overwrite = FALSE, main.name = NULL, mode, verbose){

  ###
  # CALL: Create_tTEscanR_Object() and Update_tTEscanR_Object()
  # DESCRIPTION: This function controls the addition of data and metadata into a tTEscanR object.
  ###

  results.list <- list()
  if (verbose) message(paste("- Evaluating the parameters",  deparse(substitute(data)),  "and", deparse(substitute(id)), "."))
  input_format <- IdentifyInputFormat(data = data, mode = mode) # Only 2 formats accepted: lists or character (single string)

  if (input_format == "list"){ # Dealing with a list (single or multiple elements to add)
    if (isTRUE(CheckInputCombinations(data = data, labels = id, mode = input_format))) id <- names(data) # Check if the formats of the 2 parameters are compatible

    for (i in 1:length(data)){ # Goes through each data piece in the list
      if (is.null(main.name)){
        if (isTRUE(action.update)){
          IsIn_tTEscanR_Object(object = object, slot = slot, section = id[[i]], update.assay = TRUE, overwrite = overwrite, verbose = verbose)
          if (slot == "assays") object@assays[[id[[i]]]] <- data else if (slot == "meta.data") object@meta.data[[id[[i]]]] <- data[[i]]
        } else {
          results.list[i] <- list(data[[i]])
          names(results.list)[i] <- id[[i]]
        }
      } else {

        if (isTRUE(action.update)){
          if (is.null(object@meta.data[[main.name]])) object@meta.data[[main.name]] <- list()
          IsIn_tTEscanR_Object(object = object, slot = "meta.data", section = main.name, subset = id[[i]], update.assay = TRUE, overwrite = overwrite, verbose = verbose)
          object@meta.data[[main.name]][[id[[i]]]] <- data[[i]]
        }
      }
    }
  } else { # Dealing with a single input to add
    CheckInputCombinations(data = data, labels = id, mode = "single") # Check if the formats of the 2 parameters are compatible
    if (mode =="fix" && is.data.frame(data)) data <- as.matrix(data)

    if (isTRUE(action.update)){
      IsIn_tTEscanR_Object(object = object, slot = slot, section = id, update.assay = TRUE, overwrite = overwrite, verbose = verbose)
      if (slot == "assays") object@assays[[id]] <- data else if (slot == "meta.data") object@meta.data[[id]] <- data
    } else {
      results.list <- list(data)
      names(results.list) <- id
    }
  }
  if (isTRUE(action.update)) return(object) else return(results.list) # Returns a named list with the data and their corresponding ids
}

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
  if (verbose) message("- The input contains a proper tTEscanR object.")

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

CheckInputCombinations <- function(data, labels, mode){

  ###
  # CALL: Create_tTEscanR_Object() and Update_tTEscanR_Object()
  # DESCRIPTION: This function checks all input combinations that are suitable to be added to a tTEscanR object to avoid future inconsistencies.
  # If dealing with lists we can have named or unnamed lists. Named lists already have an identifier for each piece of data.
  # However, unnamed lists will require an extra input parameter to specify the identifiers.
  ###

  if (mode == "list"){ # Dealing with a list (single or multiple elements to add)

    if (is.null(labels)){ # No extra input parameter given - the list needs to be NAMED
      if (is.null(names(data)) || all(names(data) == "")) stop(paste("The parameter is an unnamed", mode, ".\n"), paste("Please provide a named", mode, "or specify a suitable set of labels."))
      return(TRUE) # Only returns TRUE if the list is named and the names need to be assigned to the ids variable

    } else { # Extra input parameter given - the list need to be UNNAMED
      # It could happen that we have a named list and extra input parameter with the labels
      if(!is.null(names(data))) stop(paste("The parameter is a named", mode, "but also a labels parameter has been input.\n"), "Please identify the data just in one way.")

      input_format_labels <- IdentifyInputFormat(data = labels, mode = "flexible")
      if (input_format_labels == "list") { # Dealing with a list (single or multiple elements to add)

        if (length(data) != length(labels)) stop("The parameters require a 1:1 relation.") # There has to be the same number of ids than data pieces

        for (i in 1:length(labels)){ # Iterates over the list of ids to check: (i) all of them are strings, (ii) there are no duplicated ids
          if (!is.character(labels[[i]])) stop(paste("Element in `meta.data.ids` is not a string.\n Found class:", class(labels[[i]]))) # Check format of the ids

          # Looks for duplicated ids
          if (input_format_labels == "list") labels <- unlist(labels)
          duplicated_labels <- labels[duplicated(labels)]
          if (length(duplicated_labels) > 0) stop("Please provide unique labels to identify the each piece of data.\n", "Duplicated labels: ", paste(duplicated_labels, collapse = ","))
        }
      } else { # Inconsistent format between the data and the ids (i.e. if data is in a list, the ids need to be in a list)
        stop("The parameters require the same format.\n", paste("Current formats: ", mode, " and ", input_format_labels))
      }
    }
  } else if (mode == "single"){ # Dealing with a single input to add

    # When dealing just with one piece of data it is mandatory to give the extra input parameter with the ids
    if (is.null(labels) || !is.character(labels)) stop("Please provide a suitable parameter to identify the `data`.\n", "Supported formats: string.")
  }
}

IsIn_tTEscanR_Object <- function(object, slot, section, subset = NULL, compute.assay = FALSE, update.assay = FALSE, overwrite = FALSE, verbose = TRUE){

  ###
  # CALL: Multiple
  # DESCRIPTION: This function has 2 modes: (i) compute mode - check if all the inputs required are available, (ii) update mode - check if the data will need to be overwritten.
  ###

  # Extracts the name of the section in the assays or meta.data slots of a tTEscanR object
  section_list <- list(assays = list(present = names(object@assays)), meta.data = list(present = names(object@meta.data)))

  if (!(slot %in% names(section_list)) || is.null(slot)) stop("Wrong `slot` specified. Supported formats: assays and meta.data.")
  section_exists <- section %in% section_list[[slot]]$present # Check if the input `slot` is already in the tTEscanR object

  if (isTRUE(section_exists)){ # IT ALREADY EXISTS
    if (verbose) message("- A ", section, " section exists in the input object.")

    if(!is.null(subset)){ # Check the elements of a section (dealing with nested lists)

      if(subset %in% names(object@meta.data[[section]])){
        if (verbose) message("- A ", subset, " subset exists in the input object.")

        if (isTRUE(update.assay)){
          if (isFALSE(overwrite)) stop(paste("Specify another name for", subset, "or change the input parameter `overwrite` to TRUE."))
          message(paste("- The", subset, "subset will be overwritten"))
        }
      }
    } else {
      if (isTRUE(update.assay)){
        if (isFALSE(overwrite)) stop(paste("Specify another name for", section, "or change the input parameter `overwrite` to TRUE."))
        message(paste("- The", section, "section will be overwritten"))
      }
    }
  } else { # IT DOES NOT EXIST
    if (isFALSE(update.assay)){
      if (isFALSE(compute.assay)) stop(paste("A", section, "section does not exist in the input object."))
      message(paste("- A", section, "section does not exist in the input object."))
      return(compute.assay)
    } else {
      message(paste("- The", section, "section will be defined."))
    }
  }
}
