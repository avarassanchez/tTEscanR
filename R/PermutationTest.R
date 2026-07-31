#' Run Codon Usage Permutation Test (Background)
#' @description
#' This function performs a \strong{permutation test} to compare a mRNA dataset
#' to the reference codon frequency-per-gene matrix.
#'
#' @param n_permut Numeric; number of permutations to perform. Defaults to 1000.
#' @param n_features Numeric; number of features to select in each permutation.
#' Defaults to 100. If \code{target_data} is given this parameter will take its
#' length.
#' @param target_data Optional; a mRNA expression count matrix with features as
#'      rows and conditions as columns.
#' @param codon_freq Optional; a user-provided codon frequency per gene table.
#'      If necessary, it can be computed using \code{\link{getCodonFreq}}.
#' @param species Optional, a \code{character} string specifying the species
#'      reference genome version (used if \code{codon_freq} is not provided or
#'      \code{translate} is \code{TRUE}). Supported values include \code{"hg38"}
#'      (human) and \code{"mm39"} (mouse).
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#'      Defaults to \code{TRUE}.
#'
#' @return A table with the codons and their frequencies after computing all
#'      the permutations.
#' @export
#'
#' @examples
#' data(default_tTEscanR_mRNA_data)
#' genes <- default_tTEscanR_mRNA_data[1:20, ]
#' permut <- getPermutationDist(
#'     n_permut = 100, target_data = genes,
#'     species = "hg38"
#' )
getPermutationDist <- function(n_permut = 1000, n_features = 100,
    target_data = NULL, codon_freq = NULL, species = NULL, verbose = TRUE) {
    if (n_permut < 2) stop("The parameter 'n_permut' must be at least 2.")
    message(
        "--- Computing the permutation analysis ---",
        "\n1 . Loading the codon frequency per gene table."
    )
    if (!is.null(target_data)) { # TARGETED APPROACH
        codon_freq <- consistencyWithCodonFreq(
            data = target_data, codon_freq = codon_freq, species = species,
            verbose = verbose
        )
        target_names <- rownames(target_data)
        is_target <- colnames(codon_freq) %in% target_names
        if (!any(is_target)) {
            stop(
                "None of the genes in 'target_data' are ",
                "present in the codon frequency table."
            )
        }
        n_features <- sum(is_target) # Isolate background, remove targets
        codon_freq <- codon_freq[, !is_target, drop = FALSE]
    } else { # CONTROL APPROACH - ONLY TAKES THE codon_freq
        codon_freq <- checkCodonFreqTable(
            data = codon_freq, species = species, verbose = verbose
        )
    }
    codon_freq <- as.matrix(codon_freq)
    n_total_genes <- ncol(codon_freq)
    message("1 . COMPLETED\n", "2 . Starting permutation test...")
    perm_results <- replicate(n_permut, {
        index <- sample.int(n_total_genes, n_features) # Sampling of indices
        sub_sums <- rowSums(codon_freq[, index, drop = FALSE]) # Subset & sum
        return(sub_sums / sum(sub_sums)) # Normalized relative contribution
    })
    codons <- rownames(codon_freq)
    permut <- data.frame(
        codon = rep(codons, times = n_permut),
        freq = as.numeric(perm_results), stringsAsFactors = FALSE
    )
    permut <- permut[order(permut$codon, permut$freq), ]
    rownames(permut) <- NULL
    message(
        "2 . COMPLETED\n",
        "--- The permutation analysis has been successfully executed ---"
    )
    return(permut) # Returns permutation matrix
}

#' Assess Significance (P-value) & Corrects for Multiple Hypothesis Testing
#'
#' @param dist A table with the codons and their frequencies after completing a
#'      permutation test. Output from \code{\link{getPermutationDist}}.
#' @param value A \code{list} of \code{data.frame} of the codon exonic
#'      background of a mRNA gene expression matrix. The codon exonic
#'      background can be computed in \code{\link{computeCodonUsage}} using the
#'      additional_metrics parameter or directly running
#'      \code{\link{computeExonicBackground}}.
#' @param padj_threshold Numeric; p-value threshold used for highlighting
#'      significant features in the volcano plot. Defaults to 0.05.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#'      Defaults to \code{TRUE}.
#'
#' @return A table with the codon exonic background and their significance
#'      level before (p-value) and after the correction (p-adjusted value).
#' @export
#'
#' @examples
#' data(default_tTEscanR_mRNA_data, default_tTEscanR_metadata)
#' selected_genes <- default_tTEscanR_mRNA_data[1:20, ]
#' permutation_test <- getPermutationDist(
#'     n_permut = 100, target_data = selected_genes, species = "hg38"
#' )
#' tTEscanR_obj <- createObject(
#'     counts = default_tTEscanR_mRNA_data,
#'     assay = "mRNA",
#'     meta.data = list(default_tTEscanR_metadata, "tissue"),
#'     meta.data.ids = list("ConditionsLabels", "CorrectionFactor")
#' )
#' tTEscanR_obj <- computeCodonUsage(
#'     object = tTEscanR_obj, species = "hg38",
#'     additional_metrics = FALSE, reduce = 1000
#' )
#'
#' codon_usage <- getAssay(tTEscanR_obj, "CodonUsage")
#' codon_background <- rowSums(codon_usage) / sum(rowSums(codon_usage))
#' codons_to_AA <- featuresToAA(
#'     data = names(codon_background),
#'     notation_from = "codon", notation_to = "aa"
#' )
#' codon_background <- data.frame(
#'     group = codons_to_AA, codon = names(codon_background),
#'     freq = as.numeric(codon_background), row.names = NULL
#' )
#' significance <- obtainSignificance(
#'     dist = permutation_test, value = codon_background
#' )
obtainSignificance <- function(dist, value, padj_threshold = 0.05,
    verbose = TRUE) {
    if (verbose) {
        message("--- Computing the statistical significance ---")
        message("1 . Checking the input data.")
    }
    value <- as.data.frame(value)
    dist <- as.data.frame(dist)
    if (ncol(value) != 3 || ncol(dist) != 2) {
        stop(
            "Incorrect number of columns.\n", "'value' must have: group, ",
            "codon, frequency.\n", "'dist' must have: codon, frequency."
        )
    }
    colnames(value) <- c("group", "codon", "freq")
    colnames(dist) <- c("codon", "freq")
    if (!is.numeric(value$freq) || !is.numeric(dist$freq)) {
        stop("The frequency values in 'value' and 'dist' must be numeric.")
    }
    if (verbose) message("1 . COMPLETED\n2 . Calculating empirical p-values.")
    dist_list <- split(dist$freq, dist$codon) # Median of all codons as ref
    dist_medians <- vapply(dist_list, stats::median, numeric(1))
    results_list <- lapply(seq_len(nrow(value)), function(i) {
        current_codon <- value$codon[i]
        obs_freq <- as.numeric(value$freq[i])
        d_freq <- dist_list[[current_codon]]
        if (is.null(d_freq)) {
            return(data.frame(p_val = NA, tail = NA))
        }
        if (obs_freq < dist_medians[current_codon]) { # Get tail based on median
            p <- mean(d_freq <= obs_freq)
            t <- "left"
        } else {
            p <- mean(d_freq >= obs_freq)
            t <- "right"
        }
        return(data.frame(p_val = p, tail = t))
    })
    res_value <- cbind(value, do.call(rbind, results_list)) # res_df
    if (verbose) message("2 . COMPLETED\n3 . Performing FDR correction.")
    res_value$p_val_adj <- stats::p.adjust(res_value$p_val, method = "BH")
    res_value$p_val_adj <- round(res_value$p_val_adj, 4)
    res_value$sig_adj <- res_value$p_val_adj < padj_threshold
    if (verbose) {
        message("3 . COMPLETED")
        message("--- The significance has been successfully computed ---")
    }
    return(res_value)
}
