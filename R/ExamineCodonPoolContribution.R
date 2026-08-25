#' Examine the Codon Pool Contribution
#' @description
#' This function analyzes the contribution of the most highly expressed genes
#' to the overall codon pool across conditions. It is particularly useful for
#' evaluating codon bias in highly expressed genes and how it varies across
#' conditions. If needed, gene annotations can be translated for consistency,
#' and internal species-specific for human (\code{"hg38"}) and mouse
#' (\code{"mm39"}) are supported.
#'
#' @param object A \code{tTEscanR_Object} containing a mRNA and CodonUsage
#'      assay.
#' @param N Numeric; number of top genes to consider in the codon pool
#'      contribution. Defaults to 10.
#' @param corr_method A correlation method accepted by \code{\link{cor}}.
#'      Defaults to \code{"spearman"}.
#' @param codon_freq Optional; a user-provided codon frequency per gene table.
#'      If necessary, it can be computed using \code{\link{getCodonFreq}}.
#' @param species A character string specifying the species reference genome
#'      version (used if \code{codon_freq} is not provided or \code{translate}
#'      is \code{TRUE}). Supported values include \code{"hg38"} (human) and
#'      \code{"mm39"} (mouse).
#' @param overwrite Logical; if \code{TRUE}, overwrites any existing anticodon
#'      usage assay and metadata in the \code{tTEscanR_Object}. Defaults to
#'      \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#'      Defaults to \code{TRUE}.
#'
#' @return An updated \code{tTEscanR_Object} containing new layers of
#'      information in the \code{meta.data} slot, representing the codon pool
#'      contribution.
#' @export
#'
#' @examples
#' data(default_tTEscanR_mRNA_data, default_tTEscanR_metadata)
#' tTEscanR_obj <- createObject(
#'     counts = list(mRNA = default_tTEscanR_mRNA_data),
#'     meta.data = list(default_tTEscanR_metadata, "tissue"),
#'     meta.data.ids = list("ConditionsLabels", "CorrectionFactor")
#' )
#' tTEscanR_obj <- computeCodonUsage(
#'     object = tTEscanR_obj, species = "hg38", additional_metrics = FALSE
#' )
#' tTEscanR_obj <- showPoolContribution(
#'     object = tTEscanR_obj, species = "hg38"
#' )
showPoolContribution <- function(object, codon_freq = NULL, species = NULL,
    N = 10, corr_method = "spearman", overwrite = FALSE, verbose = TRUE) {
    if (verbose) message("--- Computation of the codon pool contribution ---")
    get <- extractSection(object = object, verbose = verbose)
    mRNA <- get$raw_mRNA
    if (verbose) message("2 . Checking the codon frequency per gene table.")
    codon_freq <- consistencyWithCodonFreq(
        data = mRNA, codon_freq = codon_freq, species = species
    )
    if (verbose) message("- Retrieving the mRNAs in common.")
    mRNAs_in_common <- intersect(colnames(codon_freq), rownames(mRNA))
    if (length(mRNAs_in_common) == 0) stop("No mRNAs in common found.")
    if (verbose) message("2 . COMPLETED")
    get2 <- getData(
        object = object, codon = get$cu, mRNA = mRNA, meta = get$meta,
        batch = get$batch, verbose = verbose
    )
    scc <- get2$codon # Size-corrected codon usage
    scmRNA <- get2$mRNA[mRNAs_in_common, , drop = FALSE] # Size-corrected mRNA
    codon_freq <- codon_freq[, mRNAs_in_common, drop = FALSE]
    if (verbose) message("4 . Computing each gene's codon pool contribution.")
    pool_cont <- getOrCompute(
        object, "meta.data", "CodonPoolContribution_Results",
        "CodonPoolContribution", verbose, function() {
            contrib <- as.matrix(scmRNA) * colSums(codon_freq)
            t(t(contrib) / colSums(contrib))
        }
    )
    new_meta <- helperPoolContribution(
        corr = corr_method, codon = scc, mean = get$mean_cu, verbose = verbose,
        mRNA = scmRNA, codon_freq = codon_freq, pool_cont = pool_cont, N = N
    )
    if (verbose) message("6 . Updating the object")
    res <- list()
    name <- names(object@assays)
    if (!"SizeCorrectedCodonUsage" %in% name) res$SizeCorrectedCodonUsage <- scc
    if (!"SizeCorrected_mRNA" %in% name) res$SizeCorrected_mRNA <- scmRNA
    object <- updateObject(
        overwrite = overwrite, object = object,
        counts = if (length(res) > 0) unname(res) else NULL,
        assay = if (length(res) > 0) names(res) else NULL,
        main_name = "CodonPoolContribution_Results", meta.data = new_meta
    )
    if (verbose) {
        message("6 . COMPLETED")
        message("--- The pool contribution was computed successfully ---")
    }
    return(object)
}

getData <- function(object, codon, mRNA, meta, batch, verbose) {
    if (verbose) {
        message("3 . Calculating each condition's correlation to the mean.")
    }
    codon_sc <- getOrCompute(
        object, "assays", "SizeCorrectedCodonUsage", NULL, verbose, function() {
            computeSizeCorrection(
                data = codon, metadata = meta, batch = batch, verbose = FALSE
            )
        }
    )

    data_mRNA <- getOrCompute(
        object, "assays", "SizeCorrected_mRNA", NULL, verbose, function() {
            computeSizeCorrection(
                data = mRNA, metadata = meta, batch = batch, verbose = FALSE
            )
        }
    )

    if (verbose) message("3 . COMPLETED")
    return(list(codon = codon_sc, mRNA = data_mRNA))
}

helperPoolContribution <- function(corr, codon, mean, mRNA, codon_freq,
    pool_cont, N, verbose) {
    if (verbose) {
        message(
            "4 . COMPLETED\n", "5 . Computing the codon pool diversity ",
            "with/without top ", N, " genes."
        )
    }
    correction <- computeIndividualGeneCorrelation(
        corr_method = corr, codon_usage = codon,
        mean_codon_usage = mean
    )

    res <- analizeTopGeneImpact(
        data = mRNA, codon_freq = codon_freq, corr_method = corr,
        pool_contribution = pool_cont, mean_codon_usage = mean, N = N
    )
    impact <- res$summary
    corr_topN <- round(stats::cor(
        impact$codon_diversity, correction,
        method = corr
    ), 3)
    corr_no_topN <- round(stats::cor(
        impact$sum_top_contribution,
        impact$correlation,
        method = corr
    ), 3)

    if (verbose) message("5 . COMPLETED")
    new_meta <- list(
        CodonPoolContribution = pool_cont,
        PoolContributor_NO_TopGenes = res$removed_contr,
        PoolContributorTopGenes = res$top_contributors,
        CorrelationTopGenes = corr_topN, Correlation_NO_TopGenes = corr_no_topN
    )

    new_meta <- addMeta(
        meta = new_meta, impact = impact, corr = correction, N = N
    )
    return(new_meta)
}

addMeta <- function(meta, impact, corr, N) {
    names(meta)[2:5] <- paste0(c(
        "PoolContributor_NO_Top", "PoolContributorTop",
        "CorrelationTop", "Correlation_NO_Top"
    ), N, "Genes")

    if (all(names(corr) == impact$condition)) {
        meta[[paste0("top", N, "GenesCodonPoolDiversity")]] <- data.frame(
            condition = impact$condition,
            codon_diversity = impact$codon_diversity,
            condition_correlations_to_mean_codon_usage = corr
        )

        meta[[paste0(
            "NO_top", N, "GenesCodonPoolDiversity"
        )]] <- data.frame(
            condition = impact$condition,
            codon_diversity = impact$sum_top_contribution,
            correlation = impact$correlation,
            condition_correlations_to_mean_codon_usage = corr
        )
    }
    return(metadata_list = meta)
}

extractSection <- function(object, verbose) {
    if (verbose) message("1 . Evaluating the input tTEscanR object.")
    if (!inherits(object, "tTEscanR_Object")) {
        stop("'object' must be a tTEscanR object.")
    }

    for (i in c("mRNA", "CodonUsage")) {
        isInObject(
            object = object, slot = "assays", section = i, verbose = verbose
        )
    }

    codon_usage <- getAssay(object, "CodonUsage")
    raw_mRNA <- getAssay(object, "mRNA")

    for (i in c("ConditionsLabels", "CorrectionFactor")) {
        isInObject(
            object = object, slot = "meta.data", section = i, verbose = verbose
        )
    }

    batch <- getMetadata(object, "CorrectionFactor")
    metadata <- getMetadata(object, "ConditionsLabels")
    if (!(batch %in% colnames(metadata))) {
        stop("The correction factor was not found in the metadata.")
    }

    if (verbose) message("- Extracting/Computing the mean codon usage.")
    mean_codon_usage <- getOrCompute(
        object, "meta.data", "CodonUsage_AdditionalMetrics", "MeanCodonUsage",
        verbose, function() {
            computeMeanUsage(
                data = codon_usage, assay = "CodonUsage", mode = "raw",
                metadata = metadata, batch = batch, verbose = FALSE
            )
        }
    )
    if (verbose) message("1 . COMPLETED")
    return(list(
        cu = codon_usage, raw_mRNA = raw_mRNA,
        batch = batch, meta = metadata, mean_cu = mean_codon_usage
    ))
}

getOrCompute <- function(object, slot, section, subset = NULL, verbose,
    compute_fun) {
    exists <- isInObject(
        object, slot, section, subset,
        compute_assay = TRUE, verbose = FALSE
    )
    if (isTRUE(exists)) {
        if (is.null(subset)) {
            return(slot(object, slot)[[section]])
        }
        return(slot(object, slot)[[section]][[subset]])
    }
    if (verbose) message("- Computing missing component: ", section, subset)
    return(compute_fun())
}

computeIndividualGeneCorrelation <- function(codon_usage, mean_codon_usage,
    corr_method) {
    correlation_vec <- stats::cor(
        mean_codon_usage$mean_usage_across_conditions,
        y = as.matrix(codon_usage), method = corr_method
    ) # Vectorized correlation
    correlation_to_mean_codon_usage <- as.numeric(correlation_vec)
    names(correlation_to_mean_codon_usage) <- colnames(codon_usage)

    ## Returns a matrix that stores how similar the codon usage in a particular
    ## condition is to the mean codon usage across all conditions
    return(correlation_to_mean_codon_usage)
}

analizeTopGeneImpact <- function(data, codon_freq, pool_contribution,
    mean_codon_usage, N, corr_method) {
    # Identifying top N contributors based on the pool contribution
    top_contributors <- apply(pool_contribution, 2, function(x) {
        top_values <- sort(x, partial = N, decreasing = TRUE)[seq_len(N)] # [1:N]
        return(list(sum = sum(top_values), names = names(top_values)))
    })

    sums <- vapply(top_contributors, `[[`, "sum", FUN.VALUE = numeric(1))
    top_N_names <- do.call(cbind, lapply(top_contributors, `[[`, "names"))
    colnames(top_N_names) <- colnames(data)
    rownames(top_N_names) <- paste("top", seq_len(N), "gene", sep = "")

    ## Calculating correlation without top N genes
    data_copy <- as.matrix(data)
    for (i in seq_len(ncol(data_copy))) {
        data_copy[top_N_names[, i], i] <- 0
    }

    removed_top_N_genes_usage <- as.matrix(codon_freq) %*% data_copy
    removed_top_N_genes_usage <- t(
        t(removed_top_N_genes_usage) / colSums(removed_top_N_genes_usage)
    )

    corr_no_top <- computeIndividualGeneCorrelation(
        codon_usage = removed_top_N_genes_usage,
        mean_codon_usage = mean_codon_usage, corr_method = corr_method
    )
    impact_summary <- tibble::tibble(
        condition = colnames(data), sum_top_contribution = sums,
        codon_diversity = 1 - sums, correlation = corr_no_top
    )

    return(list(
        summary = impact_summary, top_contributors = top_N_names,
        removed_contr = removed_top_N_genes_usage
    ))
}
