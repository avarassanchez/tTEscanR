#' @importFrom methods setClass new validObject
NULL

#' The tTEscanR Class
#' @description
#' The \strong{tTEscanR} object is dynamically updated to store \code{assays}
#' and \code{meta.data} at each analysis step. Ensures efficient tracking and
#' organization of inputs and outputs throughout the pipeline. In order to
#' ensure robustness throughout the pipeline, specific ids have been assigned
#' and should be respected by the user.
#'
#' @slot assays A list of assays.
#' @slot meta.data A list of meta-information associated with the assays.
#'
#' @name tTEscanR_Object-class
#' @rdname tTEscanR_Object-class
#' @exportClass tTEscanR_Object

setClass(Class = "tTEscanR_Object", slots = c(
    assays = "list",
    meta.data = "list"
))

#' Get Assay Data from a tTEscanR Object
#' @description
#' This function safely retrieves the specified data from a
#' \code{tTEscanR_Object}.
#'
#' @param object A \code{tTEscanR_Object}.
#' @param name A character string specifying the name of the assay to retrieve
#'     (e.g. "mRNA", "tRNA").
#'
#' @returns The requested assay data (typically a matrix or data.frame).
#' @export
#'
#' @examples
#' data(default_tTEscanR_mRNA_data)
#' tTEscanR_obj <- createObject(
#'     counts = default_tTEscanR_mRNA_data,
#'     assay = "mRNA"
#' )
#' mRNA_data <- getAssay(tTEscanR_obj, "mRNA")
setGeneric("getAssay", function(object, name) standardGeneric("getAssay"))

#' @rdname getAssay
#' @export

setMethod("getAssay", "tTEscanR_Object", function(object, name) {
    if (!name %in% names(object@assays)) {
        stop(gettextf("Assay '%s' not found in tTEscanR_Object", name))
    }

    return(object@assays[[name]])
})

#' Get Metadata from a tTEscanR Object
#' @description
#' This function safely retrieves the specified metadata from a
#' \code{tTEscanR_Object}.
#'
#'
#' @param object A \code{tTEscanR_Object}.
#' @param name Optional; A character string specifying the name of the metadata
#'     to retrieve (e.g. "ConditionsLabels", "CorrectionFactor").
#'
#' @returns A \code{data.frame} or a vector depending of the \code{name}
#'     parameter.
#' @export
#'
#' @examples
#' data(default_tTEscanR_mRNA_data, default_tTEscanR_metadata)
#' tTEscanR_obj <- createObject(
#'     counts = default_tTEscanR_mRNA_data,
#'     assay = "mRNA",
#'     meta.data = default_tTEscanR_metadata,
#'     meta.data.ids = "ConditionsLabels"
#' )
#' conditions <- getMetadata(tTEscanR_obj, "ConditionsLabels")
setGeneric("getMetadata", function(object, name) standardGeneric("getMetadata"))

#' @rdname getMetadata
#' @export
setMethod("getMetadata", "tTEscanR_Object", function(object, name) {
    if (missing(name)) {
        return(object@meta.data)
    }
    return(object@meta.data[[name]])
})

setValidity("tTEscanR_Object", function(object) {
    ## setValidity is called by the new() used to generate a tTEscanR_Object
    errors <- character() # Empty vector to store any validation errors

    ## Set the names of the assays that can be included in a tTEscanR_Object
    valid_assay_ids <- c(
        "mRNA", "tRNA", "CodonUsage", "AnticodonUsage",
        "AADemand", "AASupply", "SizeCorrected_mRNA", "SizeCorrectedCodonUsage",
        "SizeCorrectedAADemand", "SizeCorrected_tRNA",
        "SizeCorrectedAnticodonUsage", "SizeCorrectedAASupply"
    )

    if (length(object@assays) > 0) { # Check if the assay list is not empty

        assay_names <- names(object@assays)

        if (is.null(assay_names)) { # Check all elements of the list are named
            errors <- c(errors, "The 'assays' list must have named elements.")
        } else {
            invalid_names <- setdiff(assay_names, valid_assay_ids)
            if (length(invalid_names) > 0) { # Check the validity of the names
                error_msg <- paste(
                    "Invalid 'assay' names detected:",
                    paste(invalid_names, collapse = ", "),
                    "\n Valid names are:",
                    paste(valid_assay_ids, collapse = ", ")
                )
                errors <- c(errors, error_msg)
            }

            duplicated_names <- duplicated(assay_names)
            if (any(duplicated_names)) { # Check for duplicated names
                error_msg <- paste(
                    "Duplicated 'assay' names detected:",
                    paste(unique(assay_names[duplicated_names]), collapse = ",")
                )
                errors <- c(errors, error_msg)
            }
        }
    }

    if (length(errors) == 0) {
        return(TRUE)
    } else {
        return(errors)
    }
})


#' Create a tTEscanR Object
#' @description
#' This function initializes a \code{tTEscanR_Object}, a structured container
#' designed to hold translational efficiency-related data. The object consists
#' of two main components; \code{assays}, which stores data matrices
#' (e.g. expression, codon usage), and \code{meta.data}, which stores associated
#' information and additional intermediate calculations. A
#' \code{tTEscanR_Object} can be created using a single dataset, or initialized
#' with a list of datasets provided at once. Additional assays and metadata
#' layers can be appended later using the \code{\link{updateObject}}
#' function.
#'
#' @details
#' In order to ensure robustness throughout the pipeline \strong{specific ids}
#' have been assigned and should be respected by the user.
#' \strong{\code{assays} slot:}
#' - mRNA and tRNA count matrices as \code{"mRNA"} and \code{"tRNA"}
#' - Codon and anticodon usage count matrices as \code{"CodonUsage"} and
#' \code{"AnticodonUsage"}
#' - Amino acid demand and supply count matrices as \code{"AADemand"} and
#' \code{"AASupply"}
#' - Size corrected count matrices contain the prefix \code{"SizeCorrected"}
#' added to the raw count matrices names (e.g. \code{"SizeCorrected_mRNA} or
#' \code{"SizeCorrectedCodonUsage})
#' \strong{\code{meta.data} slot:}
#' - Table with the conditions of the mRNA and tRNA data as
#' \code{"ConditionsLabels"}
#' - Active correction factor to use when running differential expression
#' analyses as \code{"CorrectionFactor"}
#' - Optional; Identifier \code{DataMetadataIndex} to indicate the column in
#' \code{"ConditionsLabels"} that contains the labels of the conditions of the
#' \code{assay}.
#'
#' @param counts A count \code{matrix} (or \code{list} of matrices) that will
#'     be stored in the \code{assays} slot.
#' @param assay Optional; a \code{character} string (or \code{list} of
#'     characters) to identify the \code{counts}. Required if \code{counts}
#'     is not a named \code{list}.
#' @param meta.data Optional; a variable (or \code{list} of variables) with
#'     additional information that will be stored in \code{meta.data} slot.
#' @param meta.data.ids Optional; a \code{character} string (or \code{list}
#'     with the labels) to identify the \code{meta.data}. Required if
#'     \code{meta.data} is not a named \code{list}.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#'     Defaults to \code{TRUE}.
#'
#' @return A \code{tTEscanR_Object}.
#' @export
#'
#' @examples
#' data(
#'     default_tTEscanR_mRNA_data, default_tTEscanR_tRNA_data,
#'     default_tTEscanR_metadata
#' )
#' tTEscanR_obj <- createObject(
#'     counts = list(
#'         mRNA = default_tTEscanR_mRNA_data,
#'         tRNA = default_tTEscanR_tRNA_data
#'     ),
#'     meta.data = default_tTEscanR_metadata,
#'     meta.data.ids = "ConditionsLabels"
#' )
createObject <- function(counts, assay = NULL, meta.data = NULL,
    meta.data.ids = NULL, verbose = TRUE) {
    if (verbose) {
        message(
            "--- Creation of a tTEscanR Object ---",
            "\n1 . Initial assessment and addition of the input data."
        )
    }
    if (verbose) message("- Adding the 'counts' to the tTEscanR object.")
    assay_list <- defineNewData(
        data = counts, id = assay, mode = "fix",
        verbose = verbose, action_update = FALSE
    )

    metadata_list <- list()
    if (!is.null(meta.data)) {
        if (verbose) message("- Adding the 'metadata' to the tTEscanR object.")
        metadata_list <- defineNewData(
            data = meta.data, id = meta.data.ids,
            mode = "flexible", verbose = verbose
        )
    }

    ## Generate the object
    object <- new(
        Class = "tTEscanR_Object", assays = assay_list,
        meta.data = metadata_list
    )
    if (verbose) {
        message(
            "1 . COMPLETED\n",
            "--- The tTEscanR Object has been successfully created ---"
        )
    }
    return(object)
}

#' tTEscanR Object Update
#' @description
#' Updates an existing \code{tTEscanR_object} using
#' \code{\link{createObject}}. For more details, refer to
#' \code{\link{createObject}}.
#'
#' @param object An existing \code{tTEscanR_Object}.
#' @param counts Optional; a count matrix (or \code{list} of matrices) that
#'     will be stored in the \code{assays} slot of the input \code{object}.
#'     Supported formats: \code{matrix}, \code{data.frame} and \code{list}.
#' @param assay Optional; a \code{character} string (or \code{list} of strings)
#'     specifying the name of the \code{counts}. Supported
#'     formats: \code{character} and \code{list}.
#' @param main_name Optional; a \code{character} string specifying the name of
#'     the \code{meta.data} if dealing with a \code{list} that needs to be
#'     added as a single element.
#' @param meta.data Optional; a \code{list} with additional information that
#'     will be stored in \code{meta.data} slot of the input \code{object}.
#' @param meta.data.ids Optional; a \code{list} with the labels to identify the
#'     \code{meta.data} (if \code{meta.data} is given).
#' @param overwrite Logical; if \code{TRUE}, overwrites any existing assay or
#'     metadata in the \code{object} if the \code{assay} or
#'     \code{meta.data.ids} label coincides. Defaults to \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#'     Defaults to \code{TRUE}.
#'
#' @return An updated \code{tTEscanR_Object}.
#' @export
#'
#' @examples
#' data(default_tTEscanR_mRNA_data, default_tTEscanR_tRNA_data)
#' tTEscanR_obj <- createObject(
#'     counts = default_tTEscanR_mRNA_data,
#'     assay = "mRNA"
#' )
#' tTEscanR_obj <- updateObject(
#'     object = tTEscanR_obj,
#'     counts = default_tTEscanR_tRNA_data,
#'     assay = "tRNA"
#' )
updateObject <- function(object, counts = NULL, assay = NULL, main_name = NULL,
    meta.data = NULL, meta.data.ids = NULL, overwrite = FALSE, verbose = TRUE) {
    if (verbose) {
        message("1 . Initial assessment and addition of the input data.")
    }
    if (!inherits(object, "tTEscanR_Object")) {
        stop("'object' must be a tTEscanR object.")
    }
    if (verbose) message("- The input consists of a proper tTEscanR object.")

    if (is.null(counts) && is.null(meta.data)) { # NOTHING TO ADD
        stop(
            "Please specify either the parameter 'counts' or 'meta.data' to ",
            "be added to the 'object'."
        )
    }

    if (!is.null(counts)) { # ADDING COUNTS
        if (verbose) message("- Adding the 'counts' to the tTEscanR object.")
        object <- defineNewData(
            object = object, slot = "assays", data = counts, id = assay,
            mode = "fix", action_update = TRUE, overwrite = overwrite,
            verbose = verbose
        )
    }

    if (!is.null(meta.data)) { # ADDING METADATA
        if (verbose) message("- Adding the 'meta.data' to the tTEscanR object.")
        object <- defineNewData(
            object = object, slot = "meta.data", data = meta.data,
            id = meta.data.ids, mode = "flexible", main_name = main_name,
            action_update = TRUE, overwrite = overwrite, verbose = verbose
        )
    }

    if (verbose) {
        message(
            "1 . COMPLETED", "\n2 . Validating and returning the ",
            "updated tTEscanR object."
        )
    }
    validObject(object = object)
    if (verbose) message("2 . COMPLETED\n")
    return(object)
}

defineNewData <- function(object = NULL, slot = NULL, data, id = NULL,
    action_update = FALSE, overwrite = FALSE, main_name = NULL, mode, verbose) {
    if (verbose) message("- Evaluating the parameters...")
    input_format <- identifyInputFormat(data = data, mode = mode)
    if (!input_format %in% c("list", "single")) {
        stop(
            "Internal problem in identifyInputFormat, returned an ",
            "invalid 'mode'."
        )
    }

    ## Check inputs and get final, validated IDs
    id <- checkInputCombinations(data = data, labels = id, mode = input_format)
    if (input_format == "single") data <- list(data)
    if (!action_update) {
        names(data) <- id
        return(data)
    }

    for (i in seq_along(data)) {
        curr_id <- id[i]
        curr_data <- data[[i]]
        if (is.null(main_name)) {
            isInObject(
                object = object, slot = slot, section = curr_id,
                update_assay = TRUE, overwrite = overwrite, verbose = verbose
            )
            if (slot == "assays") {
                object@assays[[curr_id]] <- curr_data
            } else {
                object@meta.data[[curr_id]] <- curr_data
            }
        } else {
            if (is.null(object@meta.data[[main_name]])) {
                object@meta.data[[main_name]] <- list()
            }
            isInObject(
                object = object, slot = "meta.data", section = main_name,
                subset = curr_id, update_assay = TRUE, overwrite = overwrite,
                verbose = verbose
            )
            object@meta.data[[main_name]][[curr_id]] <- curr_data
        }
    }
    return(object)
}

checkInputCombinations <- function(data, labels, mode = c("single", "list")) {
    mode <- match.arg(mode)
    if (mode == "single") data <- list(data)

    if (is.null(labels)) {
        data_names <- names(data)
        if (is.null(data_names) || all(data_names == "") ||
            any(is.na(data_names))) {
            stop(
                "Please provide a named list or specify a suitable set of ",
                "labels to identify the 'data'."
            )
        }
        labels <- data_names
    } else {
        labels <- unlist(labels)
    }

    if (length(data) != length(labels)) {
        stop(sprintf(
            "1:1 Relation required: %d items vs %d labels.",
            length(data), length(labels)
        ))
    }
    if (!is.character(labels)) {
        stop(
            "The provided 'labels' must be character strings.\n Found class: ",
            class(labels)
        )
    }

    if (anyDuplicated(labels)) {
        duplicated_labels <- unique(labels[duplicated(labels)])
        stop(
            "Please provide unique 'labels' to identify the each piece ",
            "of data.\n", "Duplicated labels: ",
            paste(duplicated_labels, collapse = ",")
        )
    }
    return(labels)
}

isInObject <- function(object, slot = c("assays", "meta.data"), section,
    subset = NULL, compute_assay = FALSE, update_assay = FALSE,
    overwrite = FALSE, verbose = TRUE) {
    slot_name <- match.arg(slot) # Validate the slots
    targets <- slot(object, slot_name) # Extracts sections inside the slots
    if (!section %in% names(targets)) { # The section does NOT exists
        if (update_assay) {
            if (verbose) message("- The ", section, " section will be defined.")
            return(invisible(TRUE))
        }
        if (compute_assay) {
            if (verbose) message("- Section ", section, " is missing.")
            return(FALSE)
        }
        stop("Section '", section, "' not found in slot '", slot_name, "'.")
    }
    if (verbose) {
        message("- Section '", section, "' exists in '", slot_name, "'.")
    }
    if (update_assay) { # Cases where the section ALREADY exists
        already_exists <- FALSE
        if (!is.null(subset)) {
            if (subset %in% names(targets[[section]])) already_exists <- TRUE
        } else {
            already_exists <- TRUE
        }
        if (already_exists && !overwrite) {
            stop(
                "Section '", section,
                "' exists. Set 'overwrite = TRUE' to replace it."
            )
        }
        if (verbose) {
            if (!is.null(subset)) {
                message(
                    "- Subset '", subset, "' in section '",
                    section, "' will be overwritten."
                )
            } else {
                message("- Section '", section, "' will be overwritten.")
            }
        }
    }
    if (compute_assay) {
        return(TRUE)
    } else {
        return(invisible(TRUE))
    }
}
