#' Examine the Codon Pool Contribution
#' @description
#' This function analyzes the contribution of the most highly expressed genes to the overall codon pool across conditions.
#' It is particularly useful for evaluating codon bias in highly expressed genes and how it varies across conditions.
#' If needed, gene annotations can be translated for consistency, and internal species-specific for human (\code{"hg38"}) and mouse (\code{"mm39"}) are supported.
#'
#' @param object A \code{tTEscanR_Object} containing a mRNA and CodonUsage assay.
#' @param N Numeric; number of top genes to consider in the codon pool contribution. Defaults to 10.
#' @param corr_method A correlation method accepted by \code{\link{cor}}. Defaults to \code{"spearman"}.
#' @param codon_freq Optional; a user-provided codon frequency per gene table. If necessary, it can be computed using \code{\link{ObtainCodonFreqPerGene}}.
#' @param species A character string specifying the species reference genome version (used if \code{codon_freq} is not provided or \code{translate} is \code{TRUE}). Supported values include \code{"hg38"} (human) and \code{"mm39"} (mouse).
#' @param filter A character string specifying how to choose among multiple transcripts per gene either \code{"canonical"} (default) or \code{"length"} (longest transcript).
#' @param overwrite.assay Logical; if \code{TRUE}, overwrites any existing anticodon usage assay in the \code{tTEscanR_Object}. Defaults to \code{FALSE}.
#' @param overwrite.metadata Logical; if \code{TRUE}, overwrites any existing related metadata in the \code{tTEscanR_Object}. Defaults to \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return An updated \code{tTEscanR_Object} containing new layers of information in the \code{meta.data} slot, representing the codon pool contribution.
#' @export
#'
#' @examples
#' data(subset_mRNA_data, metadata)
#' tTEscanR_obj <- Create_tTEscanR_Object(counts = list(mRNA = subset_mRNA_data),
#'                                        meta.data = list(metadata, "tissue"),
#'                                        meta.data.ids = list("ConditionsLabels", "CorrectionFactor"))
#' tTEscanR_obj <- ComputeCodonUsage(object = tTEscanR_obj, species = "hg38",
#'                                   additional.metrics = FALSE)
#' tTEscanR_obj <-  ExamineCodonPoolContribution(object = tTEscanR_obj, species = "hg38")

ExamineCodonPoolContribution <- function(object, codon_freq = NULL, species = NULL, filter = "canonical", N = 10, corr_method = "spearman",
                                         overwrite.assay = FALSE, overwrite.metadata = FALSE, verbose = TRUE){

  ###
  # CALL: User
  # DESCRIPTION: This function takes the mRNA expression data and the codon usage computed to evaluate the contribution of the most highly expressed genes in the codon usage pool.
  # The codon pool contribution refers to how much each gene contributes to the total codon pool in a condition.
  # All the metrics computed in this function will be stored as meta.data in the tTEscanR object.
  ###

  message("1 . Evaluating the input tTEscanR object.")
  if (!inherits(object, 'tTEscanR_Object')) stop("`object` must be a tTEscanR object.")
  if (verbose) message("- The input contains a proper tTEscanR object.")

  # Checking that the tTEscanR object contains the suitable assays
  IsIn_tTEscanR_Object(object = object, slot = "assays", section = "mRNA", verbose = verbose)
  CheckDataFrame(data = object@assays$mRNA)
  IsIn_tTEscanR_Object(object = object, slot = "assays", section = "CodonUsage", verbose = verbose)
  CheckDataFrame(data = object@assays$CodonUsage)
  IsIn_tTEscanR_Object(object = object, slot = "meta.data", section = "ConditionsLabels", verbose = verbose)
  CheckDataFrame(data = object@meta.data$ConditionsLabels)
  IsIn_tTEscanR_Object(object = object, slot = "meta.data", section = "CorrectionFactor", verbose = verbose)
  if (!(object@meta.data$CorrectionFactor %in% colnames(object@meta.data$ConditionsLabels))) stop("The correction factor was not found in the metadata.")

  # Check if the mean codon usage is available, otherwise computes it
  check_mean_codon <- IsIn_tTEscanR_Object(object = object, slot = "meta.data", section = "CodonUsage_AdditionalMetrics", subset = "MeanCodonUsage", compute.assay = TRUE, verbose = verbose)

  if (isFALSE(check_mean_codon)){
    if (verbose) message("- Computing the mean codon usage.")
    mean_codon_usage <- suppressMessages(ComputeMeanUsage(data = object@assays$CodonUsage, assay = "CodonUsage", metadata = object@meta.data$ConditionsLabels, corr_factor = object@meta.data$CorrectionFactor, verbose = FALSE))
  } else { # Retrieves the mean codon usage from the tTEscanR object
    mean_codon_usage <- object@meta.data$CodonUsage_AdditionalMetrics$MeanCodonUsage
  }

  message("  1 . COMPLETED\n", "2 . Checking the codon frequency per gene table.")
  # Loading the codon frequency and assessing consistency in gene annotation between mRNA_data and codon_freq
  codon_frequency_per_gene_table <- ConsistencyWithCodonFreq(data = object@assays$mRNA, codon_freq = codon_freq, species = species, filter = filter, verbose = FALSE)

  if (verbose) message("- Retrieving the mRNAs in common between the codon frequency and the mRNA data matrix.")
  mRNAs_in_common <- intersect(colnames(codon_frequency_per_gene_table), rownames(object@assays$mRNA))
  if (is.null(mRNAs_in_common)) stop("No mRNAs in common found between the mRNA data (object@assays$mRNA) and the codon frequency table.")

  message("  2 . COMPLETED\n", "3 . Calculate each condition's correlation to the mean codon usage across conditions.")
  if (verbose) message("- Size-correcting the codon frequency matrix.")
  codon_usage_size_corrected <- suppressMessages(ComputeSizeCorrection(data = object@assays$CodonUsage, metadata = object@meta.data$ConditionsLabels, corr_factor = object@meta.data$CorrectionFactor, verbose = FALSE))

  # Identify outliers: calculate correlation of each condition's codon usage to the mean codon usage across conditions
  condition_correlations_to_mean_codon_usage <- ComputeIndividualGeneCorrelation(codon_usage = codon_usage_size_corrected, mean_codon_usage = mean_codon_usage, corr_method = corr_method)

  message("  3 . COMPLETED\n", "4 . Computing each gene's codon pool contribution for each condition.")
  if (verbose) message("- Size-correcting the mRNA matrix.")
  data_mRNA <- suppressMessages(ComputeSizeCorrection(data = object@assays$mRNA, metadata = object@meta.data$ConditionsLabels, corr_factor = object@meta.data$CorrectionFactor, verbose = FALSE))

  if (verbose) message("- Filtering the mRNA data and the codon frequency table based on common mRNAs.")
  data_mRNA <- data_mRNA[mRNAs_in_common, ]
  codon_frequency_per_gene_table <- codon_frequency_per_gene_table[, mRNAs_in_common]
  CheckDataFrame(data = data_mRNA)
  CheckDataFrame(data = codon_frequency_per_gene_table)

  # Multiply the mRNA data with the sum of all codons in each gene
  codon_pool_contribution <- as.matrix(data_mRNA) * colSums(codon_frequency_per_gene_table)
  codon_pool_contribution <- sweep(codon_pool_contribution, 2, colSums(codon_pool_contribution), FUN = "/") # Data normalization
  CheckDataFrame(data = codon_pool_contribution)

  message("  4 . COMPLETED\n", paste("5 . Computing the codon pool diversity with top", N, "genes."))
  extract_topN_genes <- ComputeTopNGenes(data = codon_pool_contribution, N = N) # How much of the codon pool is explained by the top N genes per condition

  if (verbose) message("- Examining correlation between codon pool diversity and correlation to mean codon usage.")
  correl_codon_pool_mean <- round(stats::cor(extract_topN_genes[[1]]$codon_diversity, condition_correlations_to_mean_codon_usage, method = corr_method), 3)

  message("  5 . COMPLETED\n", paste("6 . Repeat codon pool diversity computations without top", N, "genes."))
  removed_top_correlation_to_mean_codon_usage <- ComputeWithoutTopNGenes(data = data_mRNA, codon_freq = codon_frequency_per_gene_table, mean_codon_usage = mean_codon_usage,
                                                                         extract_topN_genes = extract_topN_genes, corr_method = corr_method)

  correl_codon_pool_mean_no_topN <- round(stats::cor(removed_top_correlation_to_mean_codon_usage$codon_diversity, removed_top_correlation_to_mean_codon_usage$correlation, method = corr_method), 3)

  message("  6 . COMPLETED\n", "7 . Updating the object.")
  metadata_list <- list() # Store all the computations into a named list

  if (length(table(names(condition_correlations_to_mean_codon_usage) == extract_topN_genes[[1]]$condition)) == 1){
    condition_correlations_to_mean_codon_usage <- as.data.frame(condition_correlations_to_mean_codon_usage)
    topN_codon_pool_diversity <- cbind(extract_topN_genes[[1]], condition_correlations_to_mean_codon_usage)
    metadata_list[[paste("top", N, "GenesCodonPoolDiversity", sep = "")]] <- topN_codon_pool_diversity
  }

  # metadata_list[["ConditionsCorrelationToMeanCodonUsage"]] <- condition_correlations_to_mean_codon_usage
  # metadata_list[[paste("CodonPoolDiversityVSCorrelation_Top", N, "Genes", sep = "")]] <- extract_topN_genes[[1]]
  # metadata_list[[paste("CodonPoolContributionTop", N, "Genes", sep = "")]] <- extract_topN_genes[[1]]$codon_diversity # included in topN_codon_pool_diversity
  metadata_list[[paste("PoolContributorTop", N, "Genes", sep = "")]] <- extract_topN_genes[[2]]

  metadata_list[["CodonPoolContribution"]] <- codon_pool_contribution
  metadata_list[[paste("CorrelationTop", N, "Genes", sep = "")]] <- correl_codon_pool_mean
  metadata_list[[paste("CodonPoolDiversityVSCorrelation_NO_Top", N, "Genes", sep = "")]] <- removed_top_correlation_to_mean_codon_usage[[1]]
  metadata_list[[paste("Correlation_NO_Top", N, "Genes", sep = "")]] <- correl_codon_pool_mean_no_topN

  object <- suppressMessages(Update_tTEscanR_Object(object = object, counts = codon_usage_size_corrected, assay = "SizeCorrectedCodonUsage", main_name = "CodonPoolContribution_Results",
                                                    meta.data = metadata_list, overwrite.assay = overwrite.assay, overwrite.metadata = overwrite.metadata))
  message("  7 . COMPLETED")
  return(object) # The output tTEscanR object has been validated in Update_tTEscanR_Object()
}
