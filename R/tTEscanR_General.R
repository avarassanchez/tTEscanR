#' Runs the Theoretical Translation Efficiency (tTE) Pipeline
#' @description
#' This function wraps up all the independent functions of theoretical translation efficiency pipeline.
#' Requires an mRNA and tRNA count matrices to compute the codon and anticodon usage and further derive the amino acid demand and supply.
#' With matching condition in the former matrices the theoretical translation efficiency would be computed.
#'
#' @param mRNA_data A count matrix of mRNA genes (rows) per conditions (columns).
#' @param tRNA_data A count matrix of tRNA genes (rows) per conditions (columns).
#' @param metadata A \code{data.frame} with the meta-information related with the conditions in \code{mRNA_data} or \code{tRNA_data}. Just the intersecting conditions will be considered.
#' @param corr_factor A factor based on \code{metadata} columns to define the variable to correct for when size correcting the count matrices.
#' @param corr_method A correlation method accepted by \code{\link{cor}}. Defaults to \code{"spearman"}.
#' @param additional.metrics Logical; if \code{TRUE}, computes: (i) codon exonic background, (ii) mean codon usage, and (iii) correlation between the previous metrics. Defaults to \code{TRUE}.
#' @param compute.significance Logical; if \code{TRUE}, computes the statistical significance (p-value) of the tTE scores. Defaults to \code{TRUE}.
#' @param codon_freq Optional; a user-provided codon frequency-per-gene table. If necessary, it can be computed using \code{\link{ObtainCodonFreqPerGene}}.
#' @param species Either \code{"hg38"} (human) or \code{"mm39"} (mouse) to specify which default codon frequency-per-gene table to use. Required if \code{codon_freq} not provided or if the gene annotation is inconsistent between both inputs.
#' @param filter Either \code{"canonical"} (default) or \code{"length"} (longest transcript) to specify which transcript to choose if several are available for the same gene.
#' @param runDESeq Logical; if \code{TRUE}, performs differential expression analysis to each assay in the \code{tTEscanR_Object}. Defaults to \code{TRUE}
#' @param color_factor A factor based on \code{metadata} columns to define the colors in the PCA plots. Required if \code{runDESeq} is \code{TRUE}.
#' @param reduce Numeric; a scaling factor used to normalize large expression values that exceed R's handling capacity. Defaults to 100.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return A \code{tTEscanR_Object} with the assays and metadata computed through the tTE pipeline.
#' @export
#'
#' @examples
#' data(subset_mRNA_data, subset_tRNA_data, metadata)
#' tTEscanR_obj <- Run_tTEscanR_pipeline(mRNA_data = subset_mRNA_data, tRNA_data = subset_tRNA_data,
#'                                       metadata = metadata, species = "hg38",
#'                                       corr_factor = "tissue", additional.metrics = FALSE,
#'                                       compute.significance = FALSE, runDESeq = FALSE)

Run_tTEscanR_pipeline <- function(mRNA_data, tRNA_data, metadata,
                                  corr_factor = NULL, corr_method = "spearman",
                                  additional.metrics = TRUE, compute.significance = TRUE,
                                  codon_freq = NULL, species = NULL, filter = "canonical",
                                  reduce = 100, runDESeq = TRUE, color_factor = NULL,
                                  verbose = TRUE){
  ###
  # CALL: User
  # DESCRIPTION: This function runs the whole pipeline from the mRNA and tRNA expression data, to the computation of the tTE scores and their significance.
  # Additionally, if enabled, computes the differential expression analysis (DEA) over each of the count matrices.
  ###

  message("1 . Defining the tTEscanR object.")
  if (verbose) message("\n------------------------------")
  tTEscanR_obj <- Create_tTEscanR_Object(counts = list(mRNA = mRNA_data, tRNA = tRNA_data), meta.data = list(metadata, corr_factor), meta.data.ids = list("ConditionsLabels", "CorrectionFactor"), verbose = verbose)
  if (verbose) message("------------------------------\n")

  message("  1 . COMPLETED\n", "2 . Computing the codon usage.")
  if (verbose) message("\n------------------------------")
  tTEscanR_obj <- ComputeCodonUsage(object = tTEscanR_obj, codon_freq = codon_freq, species = species, filter = filter, additional.metrics = additional.metrics, verbose = verbose)
  if (verbose) message("------------------------------\n")

  message("  2 . COMPLETED\n", "3 . Computing the anticodon usage.")
  if (verbose) message("\n------------------------------")
  tTEscanR_obj <- ComputeAnticodonUsage(object = tTEscanR_obj, verbose = verbose)
  if (verbose) message("------------------------------\n")

  message("  3 . COMPLETED\n", "4 . Computing the amino acid demand and supply.")
  if (verbose) message("\n------------------------------")
  tTEscanR_obj <- ComputeAAUsage(object = tTEscanR_obj, level = "both", verbose = verbose)
  if (verbose) message("------------------------------\n")

  message("  4 . COMPLETED\n", "5 . Computing the tTE scores.")
  if (verbose) message("\n------------------------------")
  tTEscanR_obj <- Compute_tTE(object = tTEscanR_obj, level = "both", corr_method = corr_method, compute.significance = compute.significance, verbose = verbose)
  if (verbose) message("------------------------------\n")

  message("  5 . COMPLETED")
  count <- 6

  if(isTRUE(runDESeq)){ # Perform differential expression analysis
    message(paste(as.character(count), "Running the differential expression analysis."))
    all_DESeq2 <- list(mRNA = tTEscanR_obj@assays$mRNA, CodonUsage = tTEscanR_obj@assays$CodonUsage, AADemand = tTEscanR_obj@assays$AADemand,
                       tRNA = tTEscanR_obj@assays$tRNA, AnticodonUsage = tTEscanR_obj@assays$AnticodonUsage, AASupply = tTEscanR_obj@assays$AASupply)

    DESeq2_results <- suppressMessages(RunDEAnalysis(list_data = all_DESeq2, metadata = tTEscanR_obj@meta.data$ConditionsLabels, corr_factor = tTEscanR_obj@meta.data$CorrectionFactor))
    tTEscanR_obj <- suppressMessages(Update_tTEscanR_Object(object = tTEscanR_obj, main_name = "Results_runDESeq", meta.data = DESeq2_results))

    count <- count + 1
    message(paste("  ", as.character(count), ". COMPLETED"))
  }

  message(paste(as.character(count), "Exporting the tTEscanR object."), paste("  ", as.character(count), ". COMPLETED"))
  return(tTEscanR_obj) # The tTEscanR object is validated every time it is updated
}
