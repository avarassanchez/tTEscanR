#' Examine the Codon Pool Contribution
#' @description
#' This function analyzes the contribution of the most highly expressed genes to the overall codon pool across conditions.
#' It is particularly useful for evaluating codon bias in highly expressed genes and how it varies across conditions.
#' If needed, gene annotations can be translated for consistency, and internal species-specific for human (\code{"hg38"}) and mouse (\code{"mm39"}) are supported.
#'
#' @param object A \code{tTEscanR_Object} containing a mRNA and CodonUsage assay.
#' @param N Numeric; number of top genes to consider in the codon pool contribution. Defaults to 10.
#' @param corr_method A correlation method accepted by \code{\link{cor}}. Defaults to \code{"spearman"}.
#' @param codon_freq Optional; a user-provided codon frequency per gene table. If necessary, it can be computed using \code{\link{GetCodonFreq}}.
#' @param species A character string specifying the species reference genome version (used if \code{codon_freq} is not provided or \code{translate} is \code{TRUE}). Supported values include \code{"hg38"} (human) and \code{"mm39"} (mouse).
#' @param overwrite Logical; if \code{TRUE}, overwrites any existing anticodon usage assay and metadata in the \code{tTEscanR_Object}. Defaults to \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return An updated \code{tTEscanR_Object} containing new layers of information in the \code{meta.data} slot, representing the codon pool contribution.
#' @export
#'
#' @examples
#' data(default_tTEscanR_mRNA_data, default_tTEscanR_metadata)
#' tTEscanR_obj <- Create_tTEscanR_Object(counts = list(mRNA = default_tTEscanR_mRNA_data),
#'                                        meta.data = list(default_tTEscanR_metadata, "tissue"),
#'                                        meta.data.ids = list("ConditionsLabels", "CorrectionFactor"))
#' tTEscanR_obj <- ComputeCodonUsage(object = tTEscanR_obj, species = "hg38",
#'                                   additional_metrics = FALSE)
#' tTEscanR_obj <-  ShowPoolContribution(object = tTEscanR_obj, species = "hg38")

ShowPoolContribution <- function(object, codon_freq = NULL, species = NULL, N = 10, corr_method = "spearman",
                                 overwrite = FALSE, verbose = TRUE){

  message("--- Computation of the codon pool contribution ---", "\n1 . Evaluating the input tTEscanR object.")
  if (!inherits(object, 'tTEscanR_Object')) stop("'object' must be a tTEscanR object.")

  for (i in c("mRNA", "CodonUsage")) IsIn_tTEscanR_Object(object = object, slot = "assays", section = i, verbose = verbose)
  codon_usage <- object@assays$CodonUsage
  raw_mRNA <- object@assays$mRNA

  for (i in c("ConditionsLabels", "CorrectionFactor")) IsIn_tTEscanR_Object(object = object, slot = "meta.data", section = i, verbose = verbose)
  batch <- object@meta.data$CorrectionFactor
  metadata <- object@meta.data$ConditionsLabels
  if (!(batch %in% colnames(metadata))) stop("The correction factor was not found in the metadata.")

  if (verbose) message("- Extracting/Computing the mean codon usage.")
  mean_codon_usage <- GetOrCompute(object, "meta.data", "CodonUsage_AdditionalMetrics", "MeanCodonUsage", verbose, function() {
    ComputeMeanUsage(data = codon_usage, assay = "CodonUsage", mode = "raw", metadata = metadata, batch = batch, verbose = FALSE)
  })

  message("1 . COMPLETED\n", "2 . Checking the codon frequency per gene table.")
  codon_frequency_per_gene_table <- ConsistencyWithCodonFreq(data = raw_mRNA, codon_freq = codon_freq, species = species, verbose = FALSE)

  if (verbose) message("- Retrieving the mRNAs in common between the codon frequency and the mRNA data matrix.")
  mRNAs_in_common <- intersect(colnames(codon_frequency_per_gene_table), rownames(raw_mRNA))
  if (length(mRNAs_in_common) == 0) stop("No mRNAs in common found.")

  message("2 . COMPLETED\n", "3 . Calculating each condition's correlation to the mean codon usage.")
  codon_usage_size_corrected <- GetOrCompute(object, "assays", "SizeCorrectedCodonUsage", NULL, verbose, function() {
    ComputeSizeCorrection(data = codon_usage, metadata = metadata, batch = batch, verbose = FALSE)
  })

  data_mRNA <- GetOrCompute(object, "assays", "SizeCorrected_mRNA", NULL, verbose, function() {
    ComputeSizeCorrection(data = raw_mRNA, metadata = metadata, batch = batch, verbose = FALSE)
  })

  message("3 . COMPLETED\n", "4 . Computing each gene's codon pool contribution.")
  data_mRNA <- data_mRNA[mRNAs_in_common, , drop = FALSE]
  codon_frequency_per_gene_table <- codon_frequency_per_gene_table[, mRNAs_in_common, drop = FALSE]

  codon_pool_contribution <- GetOrCompute(object, "meta.data", "CodonPoolContribution_Results", "CodonPoolContribution", verbose, function() {
    contrib <- as.matrix(data_mRNA) * colSums(codon_frequency_per_gene_table)
    t(t(contrib) / colSums(contrib))
  })

  message("4 . COMPLETED\n", paste("5 . Computing the codon pool diversity with/without top", N, "genes."))
  indiv_corr <- ComputeIndividualGeneCorrelation(codon_usage = codon_usage_size_corrected, mean_codon_usage = mean_codon_usage, corr_method = corr_method)

  impact_results <- AnalizeTopGeneImpact(data = data_mRNA, codon_freq = codon_frequency_per_gene_table, pool_contribution = codon_pool_contribution,
                                         mean_codon_usage = mean_codon_usage, N = N, corr_method = corr_method)

  impact_summary <- impact_results$summary
  correl_topN <- round(stats::cor(impact_summary$codon_diversity, indiv_corr, method = corr_method), 3)
  correl_no_topN <- round(stats::cor(impact_summary$sum_top_contribution, impact_summary$correlation, method = corr_method), 3)

  message("5 . COMPLETED\n", "6 . Updating the object")
  metadata_list <- list(CodonPoolContribution = codon_pool_contribution,
                        PoolContributor_NO_TopGenes = impact_results$removed_contr,
                        PoolContributorTopGenes = impact_results$top_contributors,
                        CorrelationTopGenes = correl_topN,
                        Correlation_NO_TopGenes = correl_no_topN)

  names(metadata_list)[2:5] <- paste0(c("PoolContributor_NO_Top", "PoolContributorTop", "CorrelationTop", "Correlation_NO_Top"), N, "Genes")

  if (all(names(indiv_corr) == impact_summary$condition)) {
    metadata_list[[paste0("top", N, "GenesCodonPoolDiversity")]] <- data.frame(
      condition = impact_summary$condition,
      codon_diversity = impact_summary$codon_diversity,
      condition_correlations_to_mean_codon_usage = indiv_corr)

    metadata_list[[paste0("NO_top", N, "GenesCodonPoolDiversity")]] <- data.frame(
      condition = impact_summary$condition,
      codon_diversity = impact_summary$sum_top_contribution,
      correlation = impact_summary$correlation,
      condition_correlations_to_mean_codon_usage = indiv_corr)
  }

  assays_to_update <- list()
  if (!"SizeCorrectedCodonUsage" %in% names(object@assays)) {
    assays_to_update$SizeCorrectedCodonUsage <- codon_usage_size_corrected
  }
  if (!"SizeCorrected_mRNA" %in% names(object@assays)) {
    assays_to_update$SizeCorrected_mRNA <- data_mRNA
  }

  object <- suppressMessages(Update_tTEscanR_Object(object = object, counts = if (length(assays_to_update) > 0) unname(assays_to_update) else NULL,
                                                    assay = if (length(assays_to_update) > 0) names(assays_to_update) else NULL,
                                                    main_name = "CodonPoolContribution_Results", meta.data = metadata_list,
                                                    overwrite = overwrite))

  message("6 . COMPLETED\n", "--- The computation of the codon pool contribution was performed successfully ---")
  return(object)
}

GetOrCompute <- function(object, slot, section, subset = NULL, verbose, compute_fun) {

  ###
  # CALL: ShowPoolContribution()
  # DESCRIPTION: This function checks is the requested data already exists in the input object or if it needs to be computed.
  ###

  exists <- IsIn_tTEscanR_Object(object, slot, section, subset, compute_assay = TRUE, verbose = FALSE)
  if (isTRUE(exists)) {
    if (is.null(subset)) return(slot(object, slot)[[section]])
    return(slot(object, slot)[[section]][[subset]])
  }
  if (verbose) message(paste("- Computing missing component:", section, subset))
  return(compute_fun())
}

ComputeIndividualGeneCorrelation <- function(codon_usage, mean_codon_usage, corr_method){

  ###
  # CALL: ShowPoolContribution()
  # DESCRIPTION: Identification of outliers: calculate correlation of each condition's codon usage to the mean codon usage across conditions.
  # This function takes as input the size corrected codon usage matrix and the mean codon usage.
  # It correlates each individual codon usage value with the mean codon usage (reference).
  ###

  correlation_vec <- stats::cor(mean_codon_usage$mean_usage_across_conditions, y = as.matrix(codon_usage), method = corr_method) # Vectorized correlation
  correlation_to_mean_codon_usage <- as.numeric(correlation_vec)
  names(correlation_to_mean_codon_usage) <- colnames(codon_usage)

  return(correlation_to_mean_codon_usage) # Returns a matrix that stores how similar the codon usage in a particular condition is to the mean codon usage across all conditions
}

AnalizeTopGeneImpact <- function(data, codon_freq, pool_contribution, mean_codon_usage, N, corr_method){

  # Identifying top N contributors based on the pool contribution
  top_contributors <- apply(pool_contribution, 2, function(x) {
    top_values <- sort(x, decreasing = TRUE)[1:N]
    return(list(sum = sum(top_values), names = names(top_values)))
  })

  sums <- sapply(top_contributors, `[[`, "sum")
  top_N_names <- do.call(cbind, lapply(top_contributors, `[[`, "names"))
  colnames(top_N_names) <- colnames(data)
  rownames(top_N_names) <- paste("top", seq_len(N), "gene", sep = "") # Restored old naming

  # Calculating correlation without top N genes
  data_copy <- as.matrix(data)
  for (i in seq_len(ncol(data_copy))) {
    data_copy[top_N_names[, i], i] <- 0
  }

  removed_top_N_genes_usage <- as.matrix(codon_freq) %*% data_copy
  removed_top_N_genes_usage <- t(t(removed_top_N_genes_usage) / colSums(removed_top_N_genes_usage))

  corr_no_top <- ComputeIndividualGeneCorrelation(codon_usage = removed_top_N_genes_usage, mean_codon_usage = mean_codon_usage, corr_method = corr_method)
  impact_summary <- tibble::tibble(condition = colnames(data), sum_top_contribution = sums, codon_diversity = 1 - sums, correlation = corr_no_top)

  return(list(summary = impact_summary, top_contributors = top_N_names, removed_contr = removed_top_N_genes_usage))
}
