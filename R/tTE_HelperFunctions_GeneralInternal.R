###
# General functions used internally by the tTEscanR main functions
# 3 parts: (i) checking functions
#          (ii) input data control
#          (iii) translation functions
###

# (i) Checking functions

CheckCodonFreqTable <- function(data, species, verbose = TRUE) {
    if (is.null(data)) { # No user-defined codon frequency table has been given
        if (is.null(species)) {
            stop(
                "No 'codon_freq' provided and no 'species' specified.\n",
                "Specify 'hg38' for human or 'mm39' for mouse default data."
            )
        }

        # Extract from the tTEscanR memory the default codon frequency table
        return(SelectDefaultData(species = species, verbose = verbose))
    }

    # A user-defined codon-frequency table has been given
    # Check that the format of the table is suitable
    CheckDataFrame(data = data, required_names = TRUE)
    actual_codons <- rownames(data)

    return(data) # Returns the codon frequency per gene table
}

CheckDataFrame <- function(data, required_names = TRUE) {
    # Control - in case the data object is empty
    if (is.null(data)) stop("The dataset is empty or could not be found.")
    dims <- dim(data)

    # We use at least 2 columns and 2 rows - evaluate it is not empty
    if (is.null(dims) || dims[1] < 1 || dims[2] < 1) {
        row_val <- ifelse(is.null(dims), 0, dims[1])
        col_val <- ifelse(is.null(dims), 0, dims[2])
        stop(sprintf(
            "Incorrect dimensions. Expected at least 2x2. Found: %d x %d",
            row_val, col_val
        ))
    }

    # required_names is used to specify if the data needs row and column names
    if (required_names) {
        row_names <- rownames(data)
        col_names <- colnames(data)
        if (is.null(row_names) || is.null(col_names)) {
            stop(
                "Dataset labels missing. ",
                "Both rownames and colnames are required."
            )
        }

        if (any(row_names == "", na.rm = TRUE) ||
            any(col_names == "", na.rm = TRUE)) {
            stop(
                "Dataset labels are incomplete. ",
                "Empty strings detected in rownames or colnames."
            )
        }
    }
    return(invisible(NULL))
}

CheckGeneAnnotation <- function(vector1, vector2, verbose = TRUE) {
    # Check in which format each vector of genes is
    if (verbose) message("- Vector 1: codon frequency per gene table.")
    format1 <- IsEnsemblID(gene_vector = vector1, verbose = verbose)
    if (verbose) message("- Vector 1 was properly loaded\n")

    if (verbose) message("- Vector 2: mRNA gene expression data.")
    format2 <- IsEnsemblID(gene_vector = vector2, verbose = verbose)
    if (verbose) message("- Vector 2 was properly loaded")

    if (format1 != format2) {
        stop(
            "Annotation mismatch detected!\n",
            "- Codon frequency table: ", format1,
            "\n", "- mRNA data: ", format2, ".\n",
            "Please ensure both datasets use the same gene naming convention."
        )
    }

    if (verbose) message("Both vectors are in the same format: ", format1, ".")
}

CheckIntegerLength <- function(data, reduce, verbose) {
    if (any(data > .Machine$integer.max)) { # R limit in the integer length

        data <- round(data / reduce)
        CheckDataFrame(data = data)

        if (verbose) {
            message(
                "- Some values exceed R's integer limit. Matrix divided by ",
                reduce, "."
            )
        }
    }

    return(data)
}

# (ii) Input data control

SelectDefaultData <- function(species, action = "codon_freq", verbose = TRUE) {
    valid_species <- c("hg38", "mm39")
    if (!(species %in% valid_species)) {
        stop(
            "Incorrect 'species'.\n",
            "Supported formats: hg38 for human and mm39 for mouse."
        )
    }

    # Construct the name of the object to retrieve
    if (action == "codon_freq") {
        if (verbose) {
            message("- The default ", species, " 'codon_freq' will be used.")
        }
        data_name <- paste("codon_freq_table_canonical", species, sep = "_")
    }

    if (action == "confidence_set") {
        if (verbose) {
            message("The default ", species, " 'confidence_set' will be used.")
        }
        data_name <- paste(species, "set", sep = "_")
    }

    if (action == "tRNA_map") {
        if (verbose) {
            message("The default ", species, " 'tRNA_name_map' will be used.")
        }
        data_name <- paste(species, "map", sep = "_")
    }

    data <- get(data_name, envir = asNamespace("tTEscanR"))
    return(data) # Return the retrieved data
}

IdentifyInputFormat <- function(data, mode = c("fix", "flexible")) {
    mode <- match.arg(mode) # Validate the string arguments

    is_table_like <- function(x) {
        is.data.frame(x) || is.matrix(x) || inherits(x, "dgCMatrix")
    }

    if (is.list(data) && !is.data.frame(data)) {
        if (mode == "fix") {
            valid_types <- vapply(data, is_table_like, logical(1))
            if (!all(valid_types)) {
                stop(
                    "Incorrect 'data' format.\nSupported formats: ",
                    "list of dataframes or matrices."
                )
            }
        }
        return("list")
    }

    if (is_table_like(data)) {
        return("single")
    }
    if (is.character(data)) {
        return("single")
    }
    if (is.numeric(data) && mode == "flexible") {
        return("single")
    }

    if (mode == "fix") {
        stop(
            "Incorrect 'data' format.\n", "Supported formats: data.frame, ",
            "matrix, sparse matrix, or a list of these."
        )
    } else {
        stop(
            "Incorrect 'data' format.\n", "Supported formats: table-like ",
            "objects, numeric/chracter vectors, or lists."
        )
    }
}

FilterByMetadata <- function(data, metadata, id_col = NULL, verbose = TRUE) {
    conditions <- colnames(data)

    if (is.null(id_col)) {
        matching_counts <- vapply(
            metadata, function(col) sum(conditions %in% col), numeric(1)
        )
        matching_col <- names(which.max(matching_counts))
        if (matching_counts[matching_col] == 0) {
            stop("No matching Sample IDs found between 'data' and 'metadata.'")
        }
        if (verbose) {
            message(sprintf("- Detected metadata column: '%s'", matching_col))
        }
    } else {
        matching_col <- id_col
    }

    target_col <- metadata[[matching_col]]
    data_dups <- unique(conditions[duplicated(conditions)])
    meta_dups <- unique(target_col[duplicated(target_col)])
    all_forbidden <- unique(c(data_dups, meta_dups))
    if (verbose && length(all_forbidden) > 0) {
        message(
            "- Found ", length(all_forbidden),
            " duplicated IDs. Removing all instances from analysis."
        )
    }

    shared_unique <- setdiff(intersect(conditions, target_col), all_forbidden)
    if (length(shared_unique) == 0) {
        stop(
            "No unique, non-duplicated matching IDs found between 'data' ",
            "and 'metadata.'"
        )
    }

    filtered_data <- data[, shared_unique, drop = FALSE]
    reorder_idx <- match(shared_unique, target_col)
    filtered_metadata <- metadata[reorder_idx, , drop = FALSE]
    rownames(filtered_metadata) <- shared_unique

    if (verbose) {
        message(
            "- Final dataset contains ", ncol(filtered_data), " unique samples."
        )
    }

    return(list(data = filtered_data, metadata = filtered_metadata))
}

# (iii) Translation functions

PerformTranslation <- function(
    data, notation_from, notation_to, genetic_code, verbose
) {
    data <- gsub("U", "T", data)

    col_from <- switch(
        notation_from, "codon" = "Codon", "anticodon" = "Anticodon"
    )
    col_to <- switch(
        notation_to, "codon" = "Codon", "anticodon" = "Anticodon",
        "aa" = genetic_code
    )
    if (notation_to == "aa" &&
        !(col_to %in% colnames(final_matrix_genetic_code))) {
        stop(
            "The genetic code '", genetic_code,
            "' is not a valid column in final_matrix_genetic_code."
        )
    }

    dict_from <- final_matrix_genetic_code[[col_from]]
    dict_to <- final_matrix_genetic_code[[col_to]]

    index <- match(data, dict_from)
    valid_index <- !is.na(index)
    matches <- sum(valid_index)

    if (matches == 0) {
        stop("No overlapping features identified for ", notation_from)
    }
    if (verbose) {
        message(sprintf(
            "- %d out of %d will be translated.",
            matches, length(data)
        ))
    }

    translated_features <- data.frame(
        from = data, to = dict_to[index], stringsAsFactors = FALSE
    )
    colnames(translated_features) <- c(notation_from, notation_to)

    # Return a table with the input and the translated features
    return(translated_features)
}

RetrieveTranslation <- function(
    format_input, position, translated_features, data_to_translate,
    notation_to, notation_from
) {
    new_values <- translated_features[[notation_to]]
    has_translation <- !is.na(new_values) # Update names with valid translation

    if (format_input) {
        if (position == "row") {
            original_names <- rownames(data_to_translate)
            original_names[has_translation] <- new_values[has_translation]
            rownames(data_to_translate) <- original_names
        } else if (position == "column") {
            original_names <- colnames(data_to_translate)
            original_names[has_translation] <- new_values[has_translation]
            colnames(data_to_translate) <- original_names
        } else { # Update specific column in-place
            original_names <- data_to_translate[[position]]
            data_to_translate[[position]] <- ifelse(
                is.na(new_values), original_names, new_values
            )
        }
    } else { # Working with a vector
        data_to_translate <- ifelse(is.na(new_values), data_to_translate,
                                    new_values
        )
    }

    return(data_to_translate)
}

#' Relates Codons, Anticodons and their corresponding Amino Acids
#' @description
#' This function converts: codons and anticodons into anticodons, codons or
#' amino acids based on the genetic code.
#'
#' @param data A \code{character} vector or \code{data.frame}
#' containing gene names or features to be translated.
#' @param position Optional; either \code{"row"}, \code{"column"} or
#' \code{<column_name>} to specify the location of the genes or features in
#' \code{data}. Required if \code{data_to_transalte} is a
#' \code{data.frame}.
#' @param genetic_code A \code{character} string to specify the genetic code to
#' be used. Defaults to \code{"Standard"}.
#' @param notation_to Either \code{"codon"}, \code{"anticodon"} or \code{"aa"}
#' to select the output format of the features in \code{data}.
#' @param notation_from Either \code{"codon"} or \code{"anticodon"} to select
#' the input format of the features in \code{data_to_translate}.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#' Defaults to \code{TRUE}.
#'
#' @return Translated features (codons, anticodons or amino acids) from
#' \code{data_to_translate}.
#' @export
#'
#' @examples
#' data(default_tTEscanR_mRNA_data, default_tTEscanR_metadata)
#' tTEscanR_obj <- CreateObject(
#'     counts = default_tTEscanR_mRNA_data,
#'     assay = "mRNA",
#'     meta.data = list(default_tTEscanR_metadata, "tissue"),
#'     meta.data.ids = list("ConditionsLabels", "CorrectionFactor")
#' )
#' tTEscanR_obj <- ComputeCodonUsage(
#'     object = tTEscanR_obj, codon_freq = NULL,
#'     species = "hg38", additional_metrics = FALSE
#' )
#' codons <- rownames(getAssay(tTEscanR_obj, "CodonUsage"))
#' codons_to_AA <- FeaturesToAA(
#'     data = codons, notation_from = "codon", notation_to = "aa"
#' )
#' codons_to_anticodons <- FeaturesToAA(
#'     data = codons,
#'     notation_from = "codon", notation_to = "anticodon"
#' )
FeaturesToAA <- function(
    data, position = NULL, genetic_code = "Standard", verbose = TRUE,
    notation_from = c("codon", "anticodon"),
    notation_to = c("aa", "anticodon", "codon")
) {
    if (verbose) {
        message(
            "--- Translating the input features ---",
            "\n1 . Checking the format of the input data."
        )
    }
    notation_from <- match.arg(notation_from)
    notation_to <- match.arg(notation_to)
    if (notation_from == notation_to) {
        stop("Parameters 'notation_from' and 'notation_to' can not be equal")
    }
    is_table <- is.data.frame(data) || is.matrix(data)
    if (verbose) {
        message("1 . COMPLETED\n", "2 . Extracting the features to translate.")
    }
    if (is_table) {
        if (is.null(position)) stop("Parameter 'position' required.")
        if (!position %in% c("row", "column", colnames(data))) {
            stop(
                "Incorrect value for 'position' input.\n",
                "Supported formats: 'row', 'column', or a valid column name."
            )
        }
        features <- switch(position, "row" = rownames(data),
            "column" = colnames(data), data[[position]]
        )
    } else features <- data
    if (verbose) message("2 . COMPLETED\n", "3 . Translating the features.")
    translated_features <- PerformTranslation(
        data = features, notation_from = notation_from, verbose = verbose,
        notation_to = notation_to, genetic_code = genetic_code
    )
    data <- RetrieveTranslation(
        format_input = is_table, position = position, notation_to = notation_to,
        data_to_translate = data, notation_from = notation_from,
        translated_features = translated_features
    )
    if (verbose) {
        message(
            "3 . COMPLETED\n",
            "--- The features have been successfully translated ---"
        )
    }
    return(data)
}

IsEnsemblID <- function(gene_vector, verbose) {
    # Checking input parameters
    # if (!is.character(gene_vector) || !is.vector(gene_vector)) {
    #     stop("Input 'gene_vector' must be a character vector.")
    #}

    matches <- grepl(
        "^ENS[A-Z]*[GPT][0-9]+(\\.[0-9]+)?$",
        gene_vector, perl = TRUE
    )

    # Evaluates if the Ensembl format contains the gene version
    if (any(matches)) {
        if (any(grepl("\\.", gene_vector[matches]))) {
            if (verbose) message("- Gene versions detected.")
        }
    }

    all_ens <- all(matches)

    if (any(matches) && !all_ens) {
        stop("There are gene annotation inconsistencies in the gene vector.")
    }

    # Same annotation in all genes
    return(if (all_ens) "Ensembl ID" else "Gene Symbol")
}
