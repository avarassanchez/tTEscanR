#' Create a tTEscanR object as a MultiAssayExperiment
#' @description
#' Initializes a \code{\link[MultiAssayExperiment]{MultiAssayExperiment}} object
#' to store translation efficiency data layers.
#'
#' @param counts A matrix, data.frame or list of count matrices.
#' @param assay Optional; a character string or vector specifying the assay
#'      name(s) when \code{counts} is a matrix, data.frame or unnamed list.
#'      Defaults to \code{"mRNA"}.
#' @param meta_data Optional; a data.frame of sample metadata.
#' @param sample_id Optional; a character string naming the column in
#'      \code{meta.data} to use as sample row names.
#' @param params Optional; a list of global run parameters.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#'     Defaults to \code{TRUE}.
#'
#' @return A \code{\link[MultiAssayExperiment]{MultiAssayExperiment}} object.
#' @export
#'
#' @examples
#' data("default_tTEscanR_mRNA_data", package = "tTEscanR")
#' data("default_tTEscanR_tRNA_data", package = "tTEscanR")
#' data("default_tTEscanR_metadata", package = "tTEscanR")
#'
#' tTEscanR_obj <- createObject(
#'     counts = list(
#'         mRNA = default_tTEscanR_mRNA_data,
#'         tRNA = default_tTEscanR_tRNA_data
#'     ),
#'     meta_data = default_tTEscanR_metadata,
#'     sample_id = "conditions"
#' )
#'
createObject <- function(counts, assay = "mRNA", meta_data = NULL,
    sample_id = NULL, params = list(), verbose = TRUE) {
    if (verbose) message("--- Creation of a tTEscanR Object ---")
    counts_list <- checkCountsInput(counts = counts, assay = assay)

    if (verbose) message("- Converting counts to SummarizedExperiment objects.")
    se_list <- lapply(counts_list, function(mat) {
        SummarizedExperiment::SummarizedExperiment(
            assays = list(counts = as.matrix(mat))
        )
    })

    if (verbose) message("- Adding sample metadata to colData")
    if (!is.null(meta_data)) {
        col_data <- addColData(meta_data = meta_data, sample_id = sample_id)
    } else {
        all_samples <- unique(unlist(lapply(counts_list, colnames)))
        if (is.null(all_samples)) {
            stop(
                "'counts' must have column names when 'meta_data' is NULL."
            )
        }
        col_data <- S4Vectors::DataFrame(row.names = all_samples)
    }


    mae <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = MultiAssayExperiment::ExperimentList(se_list),
        colData = col_data, metadata = params
    )

    if (verbose) {
        message("--- The tTEscanR object has been successfully created ---")
    }
    return(mae)
}

#' Update Assay in a tTEscanR object
#'
#' @description
#' Updates an existing \code{\link[MultiAssayExperiment]{MultiAssayExperiment}}
#' object by adding or replacing one or more count matrices, sample metadata, or
#' global parameters.
#'
#' @param object An existing \code{tTEscanR_Object}.
#' @param counts A count matrix or data.frame to append or update.
#' @param assay Character string specifying the target assay name.
#' @param meta_data Optional; a data.frame of sample metadata
#' @param sample_id Optional; a character string naming the column in
#'      \code{meta.data} to use as sample row names.
#' @param params Optional; a list of global run parameters.
#' @param overwrite Logical; if \code{TRUE}, overwrites existing assays or
#'      metadata if the label coincides. Defaults to \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#'     Defaults to \code{TRUE}.
#'
#' @return An updated \code{\link[MultiAssayExperiment]{MultiAssayExperiment}}
#'     object.
#' @export
#'
#' @examples
#' data("default_tTEscanR_mRNA_data", package = "tTEscanR")
#' data("default_tTEscanR_tRNA_data", package = "tTEscanR")
#'
#' tTEscanR_obj <- createObject(
#'     counts = default_tTEscanR_mRNA_data,
#'     assay = "mRNA",
#'     sample_id = "conditions"
#' )
#' tTEscanR_obj <- updateObject(
#'     object = tTEscanR_obj,
#'     counts = default_tTEscanR_tRNA_data,
#'     assay = "tRNA"
#' )
updateObject <- function(object, counts = NULL, assay = NULL, meta_data = NULL,
    sample_id = NULL, params = NULL, overwrite = FALSE, verbose = TRUE) {
    if (verbose) message("--- Updating a tTEscanR Object ---")
    checkObject(object = object, verbose = verbose)
    if (is.null(counts) && is.null(meta_data) && is.null(params)) {
        stop("Specify 'counts', 'meta_data', or 'params' to update.")
    }
    ## Extract existing components
    exp_list <- MultiAssayExperiment::experiments(object)
    current_col_data <- MultiAssayExperiment::colData(object)
    current_metadata <- S4Vectors::metadata(object)

    if (!is.null(counts)) {
        if (verbose) message("- Updating object counts.")
        counts_list <- checkCountsInput(counts = counts, assay = assay)
        exp_list <- checkPresentData(
            action = "counts", data = exp_list, new_data = counts_list,
            overwrite = overwrite
        )
    }
    if (!is.null(meta_data)) {
        if (verbose) message("- Updating sample metadata (colData).")
        new_col_data <- addColData(meta_data = meta_data, sample_id = sample_id)
        current_col_data <- checkPresentData(
            action = "meta_data", data = current_col_data,
            new_data = new_col_data, overwrite = overwrite
        )
    }
    all_exp_samples <- unique(unlist(lapply(exp_list, colnames)))
    missing_samples <- setdiff(all_exp_samples, rownames(current_col_data))
    if (length(missing_samples) > 0) {
        new_rows <- S4Vectors::DataFrame(row.names = missing_samples)
        current_col_data <- rbind(current_col_data, new_rows)
    }
    if (!is.null(params)) {
        if (verbose) message("- Updating global metadata/parameters.")
        if (!is.list(params) || is.null(names(params))) {
            stop("'params' must be a named list.")
        }
        for (nm in names(params)) current_metadata[[nm]] <- params[[nm]]
    }
    object <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = exp_list, colData = current_col_data,
        metadata = current_metadata
    ) ## Re-construct MultiAssayExperiment
    if (verbose) {
        message("--- The tTEscanR object has been successfully updated ---")
    }
    return(object)
}

checkPresentData <- function(action, data, new_data, overwrite){

    if (action == "counts") {
        labels <- names(data)
        new_labels <- names(new_data)
    } else { # action == "meta_data"
        labels <- colnames(data)
        new_labels <- colnames(new_data)
    }
    exist <- intersect(labels, new_labels)

    if (length(exist) > 0 && !overwrite) {
        stop(
            sprintf(
                "%s already in 'object'.\nSet 'overwrite = TRUE' to replace.",
                paste(exist, collapse = ", ")
            )
        )
    }

    for (i in new_labels) {
        if (action == "counts") {
            mat <- as.matrix(new_data[[i]])
            new <- SummarizedExperiment::SummarizedExperiment(
                assays = stats::setNames(list(mat), "counts")
            )
        } else {
            new <- new_data[[i]]
        }

        data[[i]] <- new
    }

    return(data)

}

addColData <- function(meta_data, sample_id) {
    ## If the metadata does not contain rownames uses the sample_id column
    has_default_rownames <- is.null(rownames(meta_data)) ||
        identical(rownames(meta_data), as.character(seq_len(nrow(meta_data))))

    if (has_default_rownames) {
        if (is.null(sample_id)) {
            stop(
                "The 'meta_data' requires matching row names to the ",
                "columns in 'counts'. Use 'sample_id' to specify the ",
                "column to consider."
            )
        } else {
            if (sample_id %in% colnames(meta_data)){
                rownames(meta_data) <- meta_data[[sample_id]]
            } else {
                stop("The 'sample_id' is not present in 'meta_data'.")
            }
        }
    }
    col_data <- S4Vectors::DataFrame(meta_data)

    return(col_data)
}

checkCountsInput <- function(counts, assay) {
    if (is.matrix(counts) || is.data.frame(counts)) {
        if (is.null(assay)) {
            stop("As 'counts' is not a named list 'assay' is required.")
        }

        counts_list <- stats::setNames(list(counts), assay)

    } else if (is.list(counts)) {
        if (is.null(names(counts))) {
            if (is.null(assay) || length(assay) != length(counts)) {
                stop(
                    "'assay' vector length must match 'counts' list length ",
                    "when elements are unnamed."
                )
            }
            counts_list <- stats::setNames(counts, assay)
        } else {
            counts_list <- counts
        }

    } else {
        stop("'counts' must be a matrix, data.frame, or list of matrices.")
    }

    checkValidAssayNames(counts_list) # Check if the assay name is allowed

    return(counts_list)
}

checkValidAssayNames <- function(assay_names) {
    if (is.list(assay_names) && !is.data.frame(assay_names)) {
        assay_names <- names(assay_names)
    }

    if (is.null(assay_names) || length(assay_names) == 0 ||
        any(is.na(assay_names)) || any(nchar(assay_names) == 0)) {
        stop("Assay name(s) cannot be NULL, NA, empty, or unnamed.")
    }

    valid_assay_ids <- c(
        "mRNA", "tRNA", "CodonUsage", "AnticodonUsage",
        "AADemand", "AASupply", "SizeCorrected_mRNA", "SizeCorrectedCodonUsage",
        "SizeCorrectedAADemand", "SizeCorrected_tRNA",
        "SizeCorrectedAnticodonUsage", "SizeCorrectedAASupply"
    )

    invalid_names <- setdiff(assay_names, valid_assay_ids)
    if (length(invalid_names) > 0) {
        stop(
            sprintf(
                "Invalid 'assay' name(s) detected: %s\nValid names are: %s",
                paste(invalid_names, collapse = ", "),
                paste(valid_assay_ids, collapse = ", ")
            )
        )
    }

    if (anyDuplicated(assay_names)) {
        duplicated_names <- unique(assay_names[duplicated(assay_names)])
        stop(
            sprintf(
                "Duplicated 'assay' name(s) detected: %s",
                paste(duplicated_names, collapse = ", ")
            )
        )
    }
}

checkAssayPresent <- function(object, assay_name, compute = FALSE) {
    is_present <- assay_name %in% names(object)

    if (!is_present) {
        if (isFALSE(compute)) {
            stop(
                sprintf(
                    "Required assay '%s' not in 'object'. Available assays: %s",
                    assay_name,
                    if (length(names(object)) > 0) {
                        paste(sQuote(names(object)), collapse = ", ")
                    } else {
                        "none"
                    }
                )
            )
        } else {
            return(FALSE) # Assay does not exist, but will be computed
        }
    }

    return(TRUE) # Assay exists
}

checkParamPresent <- function(object, param_name, compute = FALSE) {
    meta_names <- names(S4Vectors::metadata(object))
    is_present <- param_name %in% meta_names

    if (!is_present) {
        if (isFALSE(compute)) {
            stop(
                sprintf(
                    "Required param '%s' not in 'object'. Available params: %s",
                    param_name,
                    if (length(meta_names) > 0) {
                        paste(sQuote(meta_names), collapse = ", ")
                    } else {
                        "none"
                    }
                )
            )
        } else {
            return(FALSE) # Param does not exist, but will be computed
        }
    }

    return(TRUE) # Param exists
}

checkObject <- function(object, verbose) {
    if (!inherits(object, "MultiAssayExperiment")) {
        stop("'object' must be a MultiAssayExperiment object.")
    }
    if (verbose) {
        message("- The input consists of a proper MultiAssayExperiment object.")
    }

    return(invisible(NULL))
}
