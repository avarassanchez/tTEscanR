#' Runs the Theoretical Translation Efficiency (tTE) Pipeline
#' @description
#' This function wraps up all the independent functions of theoretical translation efficiency pipeline.
#' Requires an mRNA and tRNA count matrices to compute the codon and anticodon usage and further derive the amino acid demand and supply.
#' With matching condition in the former matrices the theoretical translation efficiency would be computed.
#'
#' @param mRNA_data A count matrix of mRNA genes (rows) per conditions (columns).
#' @param tRNA_data A count matrix of tRNA genes (rows) per conditions (columns).
#' @param metadata A \code{data.frame} with the meta-information related with the conditions in \code{mRNA_data} or \code{tRNA_data}. Just the intersecting conditions will be considered.
#' @param genetic_code A \code{character} string to specify the genetic code to be used. Defaults to \code{"Standard"}.
#' @param batch A factor based on \code{metadata} columns to define the variable to correct for when size correcting the count matrices.
#' @param corr_method A correlation method accepted by \code{\link{cor}}. Defaults to \code{"spearman"}.
#' @param additional_metrics Logical; if \code{TRUE}, computes: (i) codon exonic background, (ii) mean codon usage, and (iii) correlation between the previous metrics. Defaults to \code{TRUE}.
#' @param compute_significance Logical; if \code{TRUE}, computes the statistical significance (p-value) of the tTE scores. Defaults to \code{TRUE}.
#' @param codon_freq Optional; a user-provided codon frequency-per-gene table. If necessary, it can be computed using \code{\link{GetCodonFreq}}.
#' @param species Either \code{"hg38"} (human) or \code{"mm39"} (mouse) to specify which default codon frequency-per-gene table to use. Required if \code{codon_freq} not provided or if the gene annotation is inconsistent between both inputs.
#' @param runDESeq Logical; if \code{TRUE}, performs differential expression analysis to each assay in the \code{tTEscanR_Object}. Defaults to \code{TRUE}
#' @param dim_reduct Either \code{"PCA"}, \code{"UMAP"} or \code{"tSNE"} to specify the dimensionality reduction approach to be executed if \code{runDESeq} is \code{TRUE}. Defaults to \code{"PCA"}.
#' @param color_factor A factor based on \code{metadata} columns to define the colors in the PCA plots. Required if \code{runDESeq} is \code{TRUE}.
#' @param reduce Numeric; a scaling factor used to normalize large expression values that exceed R's handling capacity. Defaults to 100.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return A \code{tTEscanR_Object} with the assays and metadata computed through the tTE pipeline.
#' @export
#'
#' @examples
#' data(default_tTEscanR_mRNA_data, default_tTEscanR_tRNA_data, default_tTEscanR_metadata)
#' tTEscanR_obj <- Run_tTEscanR_pipeline(mRNA_data = default_tTEscanR_mRNA_data,
#'                                       tRNA_data = default_tTEscanR_tRNA_data,
#'                                       metadata = default_tTEscanR_metadata,
#'                                       species = "hg38", batch = "tissue",
#'                                       additional_metrics = FALSE,
#'                                       compute_significance = FALSE, runDESeq = FALSE)

Run_tTEscanR_pipeline <- function(mRNA_data, tRNA_data, metadata, batch = NULL, corr_method = "spearman", additional_metrics = TRUE,
                                  compute_significance = TRUE, codon_freq = NULL, species = NULL, genetic_code = "Standard", dim_reduct = NULL,
                                  reduce = 100, runDESeq = TRUE, color_factor = NULL, verbose = TRUE){

  ###
  # CALL: User
  # DESCRIPTION: This function runs the whole pipeline from the mRNA and tRNA expression data, to the computation of the tTE scores and their significance.
  # Additionally, if enabled, computes the differential expression analysis (DEA) over each of the count matrices.
  ###

  message("\n--- Execution of the tTEscanR pipeline ---", "\n The following steps will be executed:\n",
          "(A) Codon usage\n(B) Anticodon usage\n(C) Amino acid demand and supply\n(D) Theoretical translation efficiency at the codon and amino acid level")
  if (runDESeq) message("(E) Differential expression analysis\n")

  tTEscanR_obj <- Create_tTEscanR_Object(counts = list(mRNA = mRNA_data, tRNA = tRNA_data), meta.data = list(metadata, batch), meta.data.ids = list("ConditionsLabels", "CorrectionFactor"), verbose = verbose)
  tTEscanR_obj <- ComputeCodonUsage(object = tTEscanR_obj, codon_freq = codon_freq, species = species, additional_metrics = additional_metrics, verbose = verbose)
  tTEscanR_obj <- ComputeAnticodonUsage(object = tTEscanR_obj, verbose = verbose)
  tTEscanR_obj <- ComputeAAUsage(object = tTEscanR_obj, level = "both", genetic_code = genetic_code, verbose = verbose)
  tTEscanR_obj <- Compute_tTE(object = tTEscanR_obj, level = "both", corr_method = corr_method, compute_significance = compute_significance, verbose = verbose)

  if(isTRUE(runDESeq)){
    target_assays <- c("mRNA", "CodonUsage", "AADemand", "tRNA", "AnticodonUsage", "AASupply")
    available_assays <- tTEscanR_obj@assays[names(tTEscanR_obj@assays) %in% target_assays]

    DESeq2_results <- suppressMessages(RunDEAnalysis(list_data = available_assays, metadata = tTEscanR_obj@meta.data$ConditionsLabels, batch = tTEscanR_obj@meta.data$CorrectionFactor, dim_reduct = dim_reduct))
    tTEscanR_obj <- Update_tTEscanR_Object(object = tTEscanR_obj, main_name = "Results_runDESeq", meta.data = DESeq2_results, verbose = verbose)
  }

  message("--- The pipeline has been properly executed ---\n")
  return(tTEscanR_obj) # The tTEscanR object is validated every time it is updated
}
