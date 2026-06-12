#' Compute Amino Acid (AA) Demand or Supply
#' @description
#' This function calculates the amino acid (AA) demand and/or supply from a
#' codon and/or anticodon usage matrices of a \code{tTEscanR_Object}. It
#' aggregates the contribution of their features based on the standard genetic
#' code (mapping codons/anticodons to AA), The resulting values reflect total
#' usage (demand) or availability (supply) of each amino acid, depending on the
#' input type.
#'
#' @details
#' In order to generalize the analysis to any organism available in Ensembl,
#' tTEscanR can be used with the 33 different genetic codes described by NCBI
#' Taxonomy:
#'
#' “Standard”, "Vertebrate Mitochondrial”, "Yeast Mitochondrial”,
#' ”Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial;
#' Mycoplasma; Spiroplasma”, "Invertebrate Mitochondrial”, "Ciliate Nuclear;
#' Dasycladacean Nuclear; Hexamita Nuclear”, "Echinoderm Mitochondrial;
#' Flatworm Mitochondrial”, ”Euplotid Nuclear”, "Bacterial, Archaeal and Plant
#' Plastid”, "Alternative Yeast Nuclear”, "Ascidian Mitochondrial”, "Alternative
#' Flatworm Mitochondrial”, "Blepharisma Macronuclear”, "Chlorophycean
#' Mitochondrial”, "Trematode Mitochondrial”, "Scenedesmus obliquus
#' mitochondrial”, "Thraustochytrium mitochondrial code”, "Rhabdopleuridae
#' Mitochondrial”, "Candidate Division SR1 and Gracilibacteria”, "Pachysolen
#' tannophilus Nuclear”, "Karyorelict Nuclear”, "Condylostoma Nuclear”,
#' "Mesodinium Nuclear”, "Peritrich Nuclear”, "Blastocrithidia Nuclear”,
#' "Balanophoraceae Plastid”, "Cephalodiscidae Mitochondrial"
#'
#' @param object A \code{tTEscanR_Object} containing codon and/or anticodon
#' usage assays to be analyzed.
#' @param level Either \code{"demand"}, \code{"supply"} or \code{"both"} to
#' indicate which analysis to perform.
#' @param genetic_code A \code{character} string to specify the genetic code to
#' be used. Defaults to \code{"Standard"}.
#' @param overwrite Logical; if \code{TRUE}, overwrites any existing assay in
#' the \code{object}. Defaults to \code{FALSE}
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#' Defaults to \code{TRUE}
#'
#' @return An updated \code{tTEscanR_Object} containing a new layer of
#' information in the \code{assays} slot representing the AA demand and/or
#' supply.
#' @export
#'
#' @examples
#' data(default_tTEscanR_tRNA_data)
#' tTEscanR_obj <- CreateObject(
#'     counts = default_tTEscanR_tRNA_data,
#'     assay = "tRNA"
#' )
#' tTEscanR_obj <- ComputeAnticodonUsage(object = tTEscanR_obj)
#' tTEscanR_obj <- ComputeAAUsage(object = tTEscanR_obj, level = "supply")
ComputeAAUsage <- function(
    object, level, genetic_code = "Standard", overwrite = FALSE, verbose = TRUE
) {
    level_name <- GeneralChecksAAUsage(
        object = object, level = level, verbose = verbose
    )
    map <- assay_map_AA[[level]]
    n <- length(map$check_assays)
    res <- vector("list", n) # Create empty lists to store the outputs

    for (i in seq_len(n)) { # Iterates by the assays that need to be retrieved
        res[[i]] <- RetrieveAAUsageData(
            object = object, genetic_code = genetic_code,
            data_section = map$check_assays[[i]],
            data_function = map$AAfunction[[i]]
        )
    }
    if (verbose) {
        message(
            "- The mRNA and/or tRNA data matrices have been properly loaded.\n",
            "1 . COMPLETED\n", "2 . Pooling counts from each amino acid."
        )
    }
    for (i in seq_along(res)) { # Iterate over each assay stored
        current_id <- map$assay_id[[i]]
        if (verbose) message("- Performing the ", current_id, " analysis.")
        AAmatrix <- GroupAA(
            data = res[[i]]$data, codons = res[[i]]$codons,
            genetic_code = genetic_code
        ) # Group features into AA
        object <- UpdateObject(
            object = object, counts = AAmatrix, assay = current_id,
            overwrite = overwrite, verbose = FALSE
        )
    }
    if (verbose) {
        message(
            "2 . COMPLETED\n", "--- The amino acid ", level_name,
            " has been successfully computed ---\n"
        )
    }
    return(object) # The object has been validated in UpdateObject()
}

GeneralChecksAAUsage <- function(object, level, verbose) {
    if (level == "both") level_name <- "usage" else level_name <- level
    if (verbose) {
        message(
            "\n--- Computation of the amino acid ", level_name,
            " ---", "\n1 . Checking the format of the input data."
        )
    }
    if (!(inherits(object, "tTEscanR_Object"))) {
        stop("'object' must be a tTEscanR object.")
    }
    if (is.null(level) || !level %in% names(assay_map_AA)) {
        stop("Specify a valid 'level' input parameter: demand, supply, both.")
    }
    if (verbose) {
        message(
            "- The input consists of a proper tTEscanR object.\n",
            "- The 'level' has been porperly specified."
        )
    }

    return(level_name)
}

GroupAA <- function(data, codons, genetic_code) {
    aa_map <- stats::setNames(
        final_matrix_genetic_code[[genetic_code]],
        final_matrix_genetic_code$Codon
    )
    AAextracted <- unname(aa_map[codons])

    if (length(AAextracted) == 0 || all(is.na(AAextracted))) {
        stop("No amino acids could be retrieved.")
    }
    if (nrow(data) != length(AAextracted)) {
        stop(sprintf(
            "Dimension mismatch: 'data' has %d rows, but %d amino ",
            "acids were extracted.", nrow(data), length(AAextracted)
        ))
    }

    # Vectorized aggregation
    # Compute sum of rows for each amino acid
    AAmatrix <- rowsum(as.matrix(data), group = AAextracted)
    # Sort the amino acids alphabetically
    AAmatrix <- AAmatrix[order(rownames(AAmatrix)), , drop = FALSE]

    CheckDataFrame(data = AAmatrix) # Ensure that the generated matrix is valid
    return(AAmatrix) # Return amino acid by condition matrix
}

RetrieveAAUsageData <- function(
    object, data_section, data_function, genetic_code
) {
    # Check if the required assays are present in the tTEscanR object
    IsInObject(
        object = object, slot = "assays",
        section = data_section, verbose = FALSE
    )
    raw_data <- getAssay(object, data_section)
    CheckDataFrame(data = raw_data) # Check that data is in a suitable format

    if (data_function == "sense") { # Extract AA from the codons - AADemand

        if (!genetic_code %in% colnames(final_matrix_genetic_code)) {
            valid_codes <- paste(colnames(final_matrix_genetic_code)[-1],
                collapse = ", "
            )
            stop(sprintf(
                "Genetic code '%s' not found. Available codes are: %s",
                genetic_code, valid_codes
            ))
        }

        # Extract the amino acid translations for the selected code
        aa_translations <- final_matrix_genetic_code[[genetic_code]]
        # Filter out the stop codons ("*") to get only the sense codons
        sense_codons <- final_matrix_genetic_code$Codon[aa_translations != "*"]
        raw_data <- raw_data[
            rownames(raw_data) %in% sense_codons, ,
            drop = FALSE
        ]

        target_codons <- rownames(raw_data)
    } else { # Extract AA from the anticodons
        target_codons <-
            as.character(Biostrings::
            reverseComplement(Biostrings::
            DNAStringSet(rownames(raw_data))))
    }

    # The output consists of list with the data, and the codon/anticodons
    return(list(data = raw_data, codons = target_codons))
}
