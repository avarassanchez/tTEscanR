#' Runs the Theoretical Translation Efficiency (tTE) Pipeline
#' @description
#' This function wraps up all the independent functions of theoretical
#' translation efficiency pipeline. Requires an mRNA and tRNA count matrices to
#' compute the codon and anticodon usage and further derive the amino acid
#' demand and supply. With matching condition in the former matrices the
#' theoretical translation efficiency would be computed.
#'
#' @param mRNA_data A count matrix of mRNA genes (rows) per conditions
#'     (columns).
#' @param tRNA_data A count matrix of tRNA genes (rows) per conditions
#'     (columns).
#' @param meta_data A \code{data.frame} with the meta-information related with
#'     the conditions in \code{mRNA_data} or \code{tRNA_data}. Just the
#'     intersecting conditions will be considered.
#' @param sample_id Optional; a character string naming the column in
#'      \code{meta.data} to use as sample row names.
#' @param genetic_code A \code{character} string to specify the genetic code
#'     to be used. Defaults to \code{"Standard"}.
#' @param batch A factor based on \code{meta_data} columns to define the
#'     variable to correct for when size correcting the count matrices.
#' @param corr_method A correlation method accepted by \code{\link{cor}}.
#'     Defaults to \code{"spearman"}.
#' @param additional_metrics Logical; if \code{TRUE}, computes: (i) codon
#'     exonic background, (ii) mean codon usage, and (iii) correlation
#'     between the previous metrics. Defaults to \code{TRUE}.
#' @param compute_significance Logical; if \code{TRUE}, computes the statistical
#'     significance (p-value) of the tTE scores. Defaults to \code{TRUE}.
#' @param codon_freq Optional; a user-provided codon frequency-per-gene table.
#'     If necessary, it can be computed using \code{\link{getCodonFreq}}.
#' @param species Either \code{"hg38"} (human) or \code{"mm39"} (mouse) to
#'     specify which default codon frequency-per-gene table to use. Required if
#'     \code{codon_freq} not provided or if the gene annotation is inconsistent
#'     between both inputs.
#' @param runDESeq Logical; if \code{TRUE}, performs differential expression
#'     analysis to each assay in the \code{tTEscanR_Object}. Defaults to
#'     \code{TRUE}
#' @param reduce Numeric; a scaling factor used to normalize large expression
#'     values that exceed R's handling capacity. Defaults to 100.
#' @param ... Additional arguments passed directly to
#'     \code{\link{runDEAnalysis}} (e.g., \code{dim_reduct},
#'     \code{compute_pairwise}, \code{color_factor}).
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#'     Defaults to \code{TRUE}.
#'
#' @return A \code{\link[MultiAssayExperiment]{MultiAssayExperiment}} with the
#'     assays and metadata computed through the tTE pipeline.
#' @export
#'
#' @examples
#' data("default_tTEscanR_mRNA_data", package = "tTEscanR")
#' data("default_tTEscanR_tRNA_data", package = "tTEscanR")
#' data("default_tTEscanR_metadata", package = "tTEscanR")
#'
#' tTEscanR_obj <- runPipeline(
#'     mRNA_data = default_tTEscanR_mRNA_data,
#'     tRNA_data = default_tTEscanR_tRNA_data,
#'     meta_data = default_tTEscanR_metadata,
#'     species = "hg38", batch = "tissue", additional_metrics = FALSE,
#'     sample_id = "conditions", compute_significance = FALSE, runDESeq = FALSE
#' )
runPipeline <- function(mRNA_data, tRNA_data, meta_data, batch,
    sample_id = NULL, corr_method = "spearman", additional_metrics = TRUE,
    runDESeq = TRUE, compute_significance = TRUE, codon_freq = NULL,
    species = NULL, genetic_code = "Standard", reduce = 100, verbose = TRUE,
    ...) {
    if (verbose) {
        message(
            "\n--- Execution of the tTEscanR pipeline ---",
            "\n The following steps will be executed:\n",
            "(A) Codon usage\n(B) Anticodon usage\n(C) Amino acid demand and ",
            "supply\n(D) Theoretical translation efficiency at the codon and ",
            "amino acid level"
        )
    }
    if (runDESeq) if (verbose) message("(E) Differential expression analysis\n")

    tTEscanR_obj <- createObject(
        counts = list(mRNA = mRNA_data, tRNA = tRNA_data),
        meta_data = meta_data, sample_id = sample_id, verbose = verbose,
        params = list("CorrectionFactor" = batch)
    )
    tTEscanR_obj <- computeCodonUsage(
        object = tTEscanR_obj, codon_freq = codon_freq, species = species,
        additional_metrics = additional_metrics, verbose = verbose
    )
    tTEscanR_obj <- computeAnticodonUsage(
        object = tTEscanR_obj, verbose = verbose
    )
    tTEscanR_obj <- computeAAUsage(
        object = tTEscanR_obj, level = "both",
        genetic_code = genetic_code, verbose = verbose
    )
    tTEscanR_obj <- computeTheoreticalTE(
        object = tTEscanR_obj, level = "both", corr_method = corr_method,
        compute_significance = compute_significance, verbose = verbose
    )
    if (isTRUE(runDESeq)) {
        tTEscanR_obj <- runDEpipeline(
            object = tTEscanR_obj, verbose = verbose, ...)
    }
    if (verbose) message("--- The pipeline has been properly executed ---\n")
    return(tTEscanR_obj)
}

runDEpipeline <- function(object, verbose, ...) {
    targets <- c(
        "mRNA", "CodonUsage", "AADemand", "tRNA", "AnticodonUsage", "AASupply"
    )

    assays <- intersect(names(object), targets)
    selected_assays_list <- MultiAssayExperiment::assays(object)[assays]

    meta <- MultiAssayExperiment::colData(object)
    batch <- S4Vectors::metadata(object)[["CorrectionFactor"]]

    DESeq2_results <- runDEAnalysis(
        list_data = selected_assays_list, metadata = meta, batch = batch,
        ..., verbose = verbose
    )

    object <- updateObject(
        object = object, params = list("Results_runDESeq" = DESeq2_results),
        verbose = verbose
    )

    return(object)
}
