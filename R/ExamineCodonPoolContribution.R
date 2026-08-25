#' Examine the Codon Pool Contribution
#' @description
#' This function analyzes the contribution of the most highly expressed genes
#' to the overall codon pool across conditions. It is particularly useful for
#' evaluating codon bias in highly expressed genes and how it varies across
#' conditions. If needed, gene annotations can be translated for consistency,
#' and internal species-specific for human (\code{"hg38"}) and mouse
#' (\code{"mm39"}) are supported.
#'
#' @param object A \code{\link[MultiAssayExperiment]{MultiAssayExperiment}}
#'      containing a \code{"mRNA"} and \code{"CodonUsage"} assay.
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
#'      usage assay and metadata in the \code{object}. Defaults to \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#'      Defaults to \code{TRUE}.
#'
#' @return An updated \code{\link[MultiAssayExperiment]{MultiAssayExperiment}}
#'      containing new layers of information in the \code{params} slot,
#'      representing the codon pool contribution,
#'      \code{"CodonPoolContribution_Results"}.
#' @export
#'
#' @examples
#' data(default_tTEscanR_mRNA_data, default_tTEscanR_metadata)
#' tTEscanR_obj <- createObject(
#'     counts = list(mRNA = default_tTEscanR_mRNA_data),
#'     meta_data = default_tTEscanR_metadata, sample_id = "conditions",
#'     params = list("CorrectionFactor" = "tissue")
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
        object, "meta.data", "CodonPoolContribution_Results", verbose,
        function() {
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
    name <- names(SummarizedExperiment::assay(object))
    if (!"SizeCorrectedCodonUsage" %in% name) res$SizeCorrectedCodonUsage <- scc
    if (!"SizeCorrected_mRNA" %in% name) res$SizeCorrected_mRNA <- scmRNA
    object <- updateObject(
        overwrite = overwrite, object = object,
        counts = if (length(res) > 0) unname(res) else NULL,
        assay = if (length(res) > 0) names(res) else NULL,
        params = list("CodonPoolContribution_Results" = new_meta)
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
        object, "assays", "SizeCorrectedCodonUsage", verbose, function() {
            computeSizeCorrection(
                data = codon, metadata = meta, batch = batch, verbose = FALSE
            )
        }
    )

    data_mRNA <- getOrCompute(
        object, "assays", "SizeCorrected_mRNA", verbose, function() {
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

    res <- analyzeTopGeneImpact(
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
    checkObject(object = object, verbose = verbose)

    checkAssayPresent(object = object, assay_name = "mRNA")
    checkAssayPresent(object = object, assay_name = "CodonUsage")

    codon_usage <- SummarizedExperiment::assay(object, "CodonUsage")
    raw_mRNA <- SummarizedExperiment::assay(object, "mRNA")

    checkParamPresent(object = object, param_name = "CorrectionFactor")
    batch <- S4Vectors::metadata(object)[["CorrectionFactor"]]

    metadata <- MultiAssayExperiment::colData(object)
    if (S4Vectors::isEmpty(metadata)) {
        stop("The 'object' does not contain metadata.")
    }

    if (!(batch %in% colnames(metadata))) {
        stop("The correction factor was not found in the metadata.")
    }

    if (verbose) message("- Extracting/Computing the mean codon usage.")
    mean_codon_usage <- getOrCompute(
        object, "meta_data", "CodonUsage_AdditionalMetrics",
        verbose, function() {
            computeMeanUsage(
                data = codon_usage, assay = "CodonUsage", mode = "raw",
                metadata = metadata, batch = batch, verbose = FALSE
            )
        }
    )
    if (inherits(mean_codon_usage, "list")) {
        mean_codon_usage <- mean_codon_usage[["MeanCodonUsage"]]
    }

    if (verbose) message("1 . COMPLETED")
    return(list(
        cu = codon_usage, raw_mRNA = raw_mRNA,
        batch = batch, meta = metadata, mean_cu = mean_codon_usage
    ))
}

getOrCompute <- function(object, slot, section, verbose, compute_fun) {
    if (slot == "meta_data") {
        exists <- checkParamPresent(
            object = object, param_name = section, compute = TRUE
        )
        if (isTRUE(exists)) {
            data <- S4Vectors::metadata(object)[[section]]
            return(data)
        }
    }

    if (slot == "assays") {
        exists <- checkAssayPresent(
            object = object, assay_name = section, compute = TRUE
        )
        if (isTRUE(exists)) {
            data <- SummarizedExperiment::assay(object, section)
            return(data)
        }
    }

    if (verbose) message("- Computing missing component: ", section)
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

analyzeTopGeneImpact <- function(data, codon_freq, pool_contribution,
    mean_codon_usage, N, corr_method) {
    pool_mat <- as.matrix(pool_contribution)
    n_rows   <- nrow(pool_mat)
    n_cols   <- ncol(pool_mat)
    ## Identifying top N contributors based on the pool contribution
    ranks <- matrixStats::colRanks(
            pool_mat, ties.method = "first", preserveShape = TRUE
        )
    row_m <- row(ranks)
    target_ranks <- seq(n_rows, n_rows - N + 1)
    top_idx <- matrix(
        vapply(
            target_ranks, function(r) row_m[ranks == r],
            FUN.VALUE = integer(n_cols)
        ), nrow = N, ncol = n_cols, byrow = TRUE
    )
    gene_names  <- rownames(pool_mat)
    top_N_names <- matrix(
        gene_names[top_idx], nrow = N, ncol = n_cols,
        dimnames = list(paste0("top", seq_len(N), "gene"), colnames(data))
    )
    row_col_idx <- cbind(as.vector(top_idx), rep(seq_len(n_cols), each = N))
    top_vals_mat <- matrix(pool_mat[row_col_idx], nrow = N, ncol = n_cols)
    sums <- matrixStats::colSums2(top_vals_mat)

    ## Calculating correlation without top N genes
    data_copy <- as.matrix(data)
    data_copy[row_col_idx] <- 0
    removed_top_N_genes_usage <- as.matrix(codon_freq) %*% data_copy
    col_sums <- matrixStats::colSums2(removed_top_N_genes_usage)
    removed_top_N_genes_usage <- sweep(
        removed_top_N_genes_usage, 2, col_sums, FUN = "/"
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
