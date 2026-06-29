#' Compute the Theoretical Translation Efficiency (tTE) score
#' @description
#' This function calculates the theoretical translation efficiency (tTE) score
#' by integrating codon-anticodon usage and/or amino acid demand-supply across
#' conditions.
#'
#' @param object A \code{tTEscanR_Object} containing a mRNA and a codon usage
#'     assay and/or a tRNA and an anticodon usage assay.
#' @param level Either \code{"codon"}, \code{"aa"} or \code{"both"} to indicate
#'     which analysis to perform.
#' @param genetic_code A \code{character} string to specify the genetic code to
#'     be used. Defaults to \code{"Standard"}.
#' @param corr_method A correlation method accepted by \code{\link{cor}}.
#'     Defaults to \code{"spearman"}.
#' @param compute_significance Logical; if \code{TRUE}, computes the
#'     statistical significance (p-value) of the tTE scores. Defaults
#'     to \code{TRUE}.
#' @param overwrite Logical; if \code{TRUE}, overwrites any existing tTE table
#'     in the \code{tTEscanR_Object}. Defaults to \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#'     Defaults to \code{TRUE}.
#'
#' @return An updated \code{tTEscanR_Object} containing a new layer of
#'     information representing the translation efficiency table for the
#'     matching conditions in the mRNA and tRNA data.
#' @export
#'
#' @examples
#' data(
#'     default_tTEscanR_mRNA_data, default_tTEscanR_tRNA_data,
#'     default_tTEscanR_metadata
#' )
#' tTEscanR_obj <- createObject(
#'     counts = list(
#'         mRNA = default_tTEscanR_mRNA_data, tRNA = default_tTEscanR_tRNA_data
#'     ),
#'     meta.data = list(default_tTEscanR_metadata, "tissue"),
#'     meta.data.ids = list("ConditionsLabels", "CorrectionFactor")
#' )
#' tTEscanR_obj <- computeCodonUsage(
#'     object = tTEscanR_obj, species = "hg38",
#'     additional_metrics = FALSE, reduce = 10000
#' )
#' tTEscanR_obj <- computeAnticodonUsage(object = tTEscanR_obj)
#' tTEscanR_obj <- computeTheoreticalTE(
#'     object = tTEscanR_obj, level = "codon", compute_significance = FALSE
#' )
computeTheoreticalTE <- function(object, level = c("codon", "aa", "both"),
    genetic_code = "Standard",
    corr_method = c("spearman", "pearson", "kendall"),
    compute_significance = TRUE, overwrite = FALSE, verbose = TRUE) {
    if (verbose) {
        message(
            "\n--- Computation of the theoretical translation efficiency ",
            "(tTE) ---"
        )
    }
    sections <- generalChecksComputetTE(object = object, verbose = verbose)
    meta <- sections$meta
    batch <- sections$batch
    corr_method <- match.arg(corr_method)
    level <- match.arg(level)
    target_levels <- if (level == "both") names(assay_map_TE) else level
    if (is.null(target_levels) ||
        !all(target_levels %in% names(assay_map_TE))) {
        stop("Please specify a valid 'level' input parameter: codon, aa, both.")
    }
    results <- list() # Empty list to store the results
    if (verbose) message("1 . COMPLETED\n", "2 . Computing the tTE analysis.")
    for (i in target_levels) {
        results <- helperComputetTE(
            level = i, results = results, verbose = verbose, meta = meta,
            object = object, batch = batch, genetic_code = genetic_code,
            sig = compute_significance, corr = corr_method
        )
    }
    if (verbose) {
        message("\n2 . COMPLETED\n", "--- The tTE was properly computed ---\n")
        message("- Updating the tTEscanR object.")
    }
    object <- updateObject(
        object = object, meta.data = results, verbose = FALSE,
        meta.data.ids = names(results), overwrite = overwrite
    )
    return(object) # tTEscanR object was validated in updateObject()
}

generalChecksComputetTE <- function(object, verbose) {
    if (verbose) message("1 . Checking the format of the input data.")
    if (!(inherits(object, "tTEscanR_Object"))) {
        stop("'object' must be a tTEscanR object.")
    }
    if (verbose) message("- The input contains a proper tTEscanR object.")

    ## Evaluate that the required data and metadata are present in the object
    isInObject(
        object = object, slot = "assays", section = "tRNA",
        update_assay = FALSE, overwrite = FALSE, verbose = FALSE
    )
    isInObject(
        object = object, slot = "meta.data", section = "ConditionsLabels"
    )
    isInObject(
        object = object, slot = "meta.data", section = "CorrectionFactor"
    )

    meta <- getMetadata(object, "ConditionsLabels")
    batch <- getMetadata(object, "CorrectionFactor")
    if (!(batch %in% colnames(meta))) {
        stop("The correction factor was not found in the metadata.")
    }

    return(list(batch = batch, meta = meta))
}

helperComputetTE <- function(results, level, verbose, object, meta, batch,
    corr, sig, genetic_code) {
    if (verbose) message("\nProcessing level: ", level)
    config <- assay_map_TE[[level]]
    mRNA <- getAssay(object, config$mRNA)
    tRNA <- getAssay(object, config$tRNA)
    checkDataFrame(data = mRNA)
    checkDataFrame(data = tRNA)
    filt1 <- filterMatrix( # Retain matching conditions
        data_mRNA = mRNA, data_tRNA = tRNA, level = level, verbose = verbose
    )
    if (verbose) message("- Filtering the metadata by intersecting conditions.")
    filt2 <- filterByMetadata( # The metadata was previously matched
        data = filt1$mRNA, metadata = meta, verbose = FALSE
    )
    if (verbose) message("- Size-correcting the data.")
    mRNA_filt <- computeSizeCorrection(
        data = filt2$data, metadata = filt2$metadata,
        batch = batch, verbose = FALSE
    )
    tRNA_filt <- computeSizeCorrection(
        data = filt1$tRNA, metadata = filt2$metadata,
        batch = batch, verbose = FALSE
    )
    tTE_res <- computeCorrelation( # Correlation mRNA and the tRNA data
        data_mRNA = mRNA_filt, data_tRNA = tRNA_filt, corr_method = corr,
        verbose = verbose
    )
    if (isTRUE(sig)) {
        if (verbose) {
            message("- Computing statistical significance (may take time)...")
        }
        ## Filter tRNA by shared cond.- tRNA used to increase num. ft shuffling
        tRNA <- getAssay(object, "tRNA")[, colnames(filt1$mRNA), drop = FALSE]
        checkDataFrame(data = tRNA)
        tRNA <- computeSizeCorrection(
            data = tRNA, metadata = filt2$metadata, batch = batch,
            verbose = FALSE
        )
        tTE_res <- computeStatisticalSignificance(
            level = level, tTE_scores = tTE_res, data_mRNA = mRNA_filt,
            genetic_code = genetic_code, data_tRNA = tRNA_filt,
            tRNA_exp = tRNA, corr_method = corr, verbose = verbose
        )
    }
    results[[config$id]] <- tTE_res
    if (verbose) message("\nCompleted the processing level: ", level)
    return(results)
}

assignMapLevel <- function(level, tRNA_exp, data_mRNA, genetic_code) {
    tRNA_names <- rownames(tRNA_exp)

    ## Extract the anticodons from the tRNA data
    tRNAsAnticodons <- vapply(
        strsplit(tRNA_names, "-"), "[[", 3,
        FUN.VALUE = character(1)
    )

    if (level == "aa") {
        if (!(genetic_code %in% colnames(final_matrix_genetic_code))) {
            stop(
                "The genetic code '", genetic_code, "' is not a valid ",
                "identifier."
            )
        }

        aa_dictionary <- stats::setNames(
            final_matrix_genetic_code[[genetic_code]],
            final_matrix_genetic_code$Anticodon
        )
        tRNA_map_level <- unname(aa_dictionary[tRNAsAnticodons])

        if (any(is.na(tRNA_map_level))) {
            warning(
                "Some anticodons mapped to codons missing from the genetic",
                " code table. They will be ignored."
            )
        }
    } else {
        codon_to_anti_dict <- stats::setNames(
            final_matrix_genetic_code$Anticodon,
            final_matrix_genetic_code$Codon
        )
        tRNA_map_level <- tRNAsAnticodons
        rownames(data_mRNA) <- unname(codon_to_anti_dict[rownames(data_mRNA)])
    }

    return(list(data_mRNA = data_mRNA, tRNA_map_level = tRNA_map_level))
}

computePvalues <- function(tTE_conditions, data_mRNA, tRNA_exp_filtered,
    N_perm, map_factor, corr_method, tTE_scores) {
    n_cond <- length(tTE_conditions)
    tTE_p_values <- numeric(n_cond) # Empty vector, store significance results
    M <- Matrix::sparse.model.matrix(~ 0 + map_factor)
    M <- Matrix::t(M) # Features as rows and columns as tRNA genes

    pb <- utils::txtProgressBar(min = 0, max = n_cond, style = 3) # Time counter
    for (i in seq_len(n_cond)) { # Iterate over conditions present in the data
        curr_cond <- tTE_conditions[i] # Get condition interrogated each time
        mRNA_vec <- data_mRNA[, curr_cond] # Observed (true) usage
        ## Observed tRNA usage: shuffled N_perm times determine stat. signif.
        tRNA_vec_raw <- tRNA_exp_filtered[, curr_cond]
        n_tRNA_genes <- length(tRNA_vec_raw)
        null_correlations <- numeric(N_perm) # Vector to store permutation val.

        for (j in seq_len(N_perm)) { # Iterates over number of iterations
            # Randomize order
            shuffled_tRNA <- tRNA_vec_raw[sample.int(n_tRNA_genes)]
            # Aggregated usage - sum of tRNA genes per feature
            agg_usage <- as.numeric(M %*% shuffled_tRNA)
            # Correlation: shuffled tRNA data and unshuffled mRNA data
            null_correlations[j] <- stats::cor( # Fills matrix with tTE score
                agg_usage, mRNA_vec,
                method = corr_method
            )
        }
        ## Fit parameters (mu and sigma for null distribution)
        mean_random_correlations <- mean(null_correlations)
        sd_random_correlations <- stats::sd(null_correlations)

        ## Calculate p-value for the actual tTE
        tTE_p_values[i] <- stats::pnorm(
            tTE_scores[i],
            mean = mean_random_correlations,
            sd = sd_random_correlations, lower.tail = FALSE
        )
        utils::setTxtProgressBar(pb, i) # Increase the progress bar
    }
    close(pb)
    tTEresult_table <- tibble::tibble(
        condition = tTE_conditions, tTE = as.numeric(tTE_scores),
        p_value = tTE_p_values, neg_log10_tTE_p_value = -log10(tTE_p_values)
    )
    return(tTEresult_table)
}

computeStatisticalSignificance <- function(level, tTE_scores, data_mRNA,
    data_tRNA, tRNA_exp, corr_method, genetic_code, N_perm = 1000, verbose) {
    map_results <- assignMapLevel(
        tRNA_exp = tRNA_exp, data_mRNA = data_mRNA,
        genetic_code = genetic_code, level = level
    )
    data_mRNA <- map_results$data_mRNA
    tRNA_map_level <- map_results$tRNA_map_level

    common_features <- intersect(rownames(data_mRNA), unique(tRNA_map_level))
    data_mRNA <- as.matrix(data_mRNA[common_features, , drop = FALSE])

    keep_tRNA_indices <- which(tRNA_map_level %in% common_features)
    tRNA_exp_filtered <- as.matrix(tRNA_exp[keep_tRNA_indices, , drop = FALSE])
    tRNA_map_level_filtered <- tRNA_map_level[keep_tRNA_indices]

    ## Initialize the output matrix
    map_factor <- factor(tRNA_map_level_filtered, levels = common_features)

    ## Generates summary table with the input tTE scores and their p-values
    tTE_conditions <- colnames(data_tRNA)
    tTE_result <- computePvalues(
        data_mRNA = data_mRNA, map_factor = map_factor, tTE_scores = tTE_scores,
        tTE_conditions = tTE_conditions, N_perm = N_perm,
        corr_method = corr_method, tRNA_exp_filtered = tRNA_exp_filtered
    )
    return(tTE_result)
}

computeCorrelation <- function(data_mRNA, data_tRNA, corr_method, verbose) {
    if (verbose) message("- Computing the ", corr_method, " correlation.")

    ## Extract the features and the conditions to consider
    ## Both datasets have been filtered for common conditions
    tTE_conditions <- colnames(data_mRNA)

    ## Selects condition in each iteration & filters the data accordingly
    tTEs_at_level <- vapply(seq_along(tTE_conditions), function(i) {
        stats::cor(data_mRNA[, i], data_tRNA[, i], method = corr_method)
    }, FUN.VALUE = numeric(1))
    names(tTEs_at_level) <- tTE_conditions

    if (any(is.na(tTEs_at_level))) {
        warning(
            "Some conditions produced NA correlation scores.",
            " Check for zero variance in those columns."
        )
    }
    if (is.null(tTEs_at_level)) stop("The tTE scores could not be computed.")
    return(tTEs_at_level) # Vector: features and correlation value (tTE score)
}

filterMatrix <- function(data_mRNA, data_tRNA, level = c("codon", "aa"),
    verbose) {
    level <- match.arg(level)
    ## Extracting the matching conditions
    tTE_conditions <- intersect(colnames(data_mRNA), colnames(data_tRNA))
    if (length(tTE_conditions) == 0) {
        stop(
            "No intersecting conditions were found between the matrices.\n",
            "If needed check that the annotation (columns) between the ",
            "matrices follows the same format."
        )
    }
    row_mRNA <- rownames(data_mRNA)
    row_tRNA <- rownames(data_tRNA)
    if (level == "codon") { # Find the proper codon-anticodon pairs
        codon_to_anti_dict <- stats::setNames(
            final_matrix_genetic_code$Anticodon, final_matrix_genetic_code$Codon
        )
        translated_codons <- unname(codon_to_anti_dict[row_mRNA])
        tTE_factors <- intersect(translated_codons, row_tRNA)
        if (length(tTE_factors) == 0) {
            stop("No intersecting features found for codons.")
        }
        valid_index <- translated_codons %in% tTE_factors
        mRNA_filtered_rows <- row_mRNA[valid_index]
        tRNA_filtered_rows <- translated_codons[valid_index]
    } else { # Dealing with AA, perform a direct intersection
        tTE_factors <- intersect(row_mRNA, row_tRNA)
        if (length(tTE_factors) == 0) {
            stop("No intersecting features found for amino acids")
        }
        mRNA_filtered_rows <- tTE_factors
        tRNA_filtered_rows <- tTE_factors
    }
    if (verbose) {
        message(
            "- Found ", length(tTE_conditions), " conditions and ",
            length(tRNA_filtered_rows), " features."
        )
    }
    if (verbose) message("- Filtering the data by the intersecting conditions.")
    filtered_mRNA <- data_mRNA[mRNA_filtered_rows, tTE_conditions, drop = FALSE]
    filtered_tRNA <- data_tRNA[tRNA_filtered_rows, tTE_conditions, drop = FALSE]
    checkDataFrame(data = filtered_mRNA)
    checkDataFrame(data = filtered_tRNA)
    return(list(mRNA = filtered_mRNA, tRNA = filtered_tRNA))
}
