#' Compute Codon Usage from mRNA Gene Expression Data
#' @description
#' This function estimates **codon usage profiles** based on gene-level mRNA expression data stored in a \code{tTEscanR_Object}.
#' It optionally accepts pre-computed codon frequency tables or uses internally generated default tables when not provided.
#' When enabled, it can evaluate the correlation between background codon composition and observed mean codon usage.
#' If the additional metrics are to be computed the input \code{tTEscanR_Object} needs to have a **"CorrectionFactor"** stored in the \code{"meta.data"} slot.
#'
#' @param object A \code{tTEscanR_Object} containing a mRNA assay.
#' @param codon_freq Optional; a user-provided codon frequency-per-gene table. If necessary, it can be computed using \code{\link{ObtainCodonFreqPerGene}}.
#' @param species Optional; either \code{"hg38"} (human) or \code{"mm39"} (mouse) to load the default settings. Required if \code{codon_freq} is not provided.
#' @param filter Optional; either \code{"canonical"} (default) or \code{"length"} (longest transcript) to specify which transcript to choose if several are available for the same gene.
#' @param reduce Numeric; a scaling factor used to normalize large expression values that exceed R's handling capacity. Defaults to 100.
#' @param additional.metrics Logical; if \code{TRUE}, computes: (i) codon exonic background, (ii) mean codon usage, and (iii) correlation between the previous metrics. Defaults to \code{TRUE}.
#' @param corr_method A correlation method accepted by \code{\link{cor}}. Required if \code{additional.metrics} is \code{TRUE}. Defaults to \code{"spearman"}.
#' @param overwrite.assay Logical; if \code{TRUE}, overwrites any existing assay in the \code{object}. Defaults to \code{FALSE}.
#' @param overwrite.metadata Logical; if \code{TRUE}, overwrites any existing metadata in the \code{object}. Defaults to \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return An updated \code{tTEscanR_Object} containing a new layer of information \code{"CodonUsage"} in the \code{assays} slot representing the codon usage.
#' Additional computations will be stored in the \code{meta.data} slot as \code{"CodonUsage_AdditionalMetrics"}.
#' @export
#'
#' @examples
#' data(subset_mRNA_data, metadata)
#' tTEscanR_obj <- Create_tTEscanR_Object(counts = subset_mRNA_data, assay = "mRNA",
#'                                        meta.data = list(metadata, "tissue"),
#'                                        meta.data.ids = list("ConditionsLabels", "CorrectionFactor"))
#' tTEscanR_obj <- ComputeCodonUsage(object = tTEscanR_obj, species = "hg38",
#'                                   additional.metrics = FALSE, reduce = 1000)

ComputeCodonUsage <- function(object, codon_freq = NULL, species = NULL, filter = "canonical", reduce = 100, additional.metrics = TRUE,
                              corr_method = "spearman", overwrite.assay = FALSE, overwrite.metadata = FALSE, verbose = TRUE){

  ###
  # CALL: User or Compute_tTE()
  # DESCRIPTION: This function takes a mRNA gene expression matrix and a codon frequency per gene table to compute the codon usage.
  # The codon usage represents how much of each codon (row) is needed for translation in a given condition (column).
  ###

  message("1 . Checking the format of the input data.")
  if (!inherits(object, 'tTEscanR_Object')) stop("`object` must be a tTEscanR object.")
  if (verbose) message("- The input contains a proper tTEscanR object.")
  IsIn_tTEscanR_Object(object = object, slot = "assays", section = "mRNA", verbose = FALSE)
  CheckDataFrame(data = object@assays$mRNA)
  if (verbose) message("- The mRNA assay has been properly loaded.")

  if (additional.metrics){
    IsIn_tTEscanR_Object(object = object, slot = "meta.data", section = "ConditionsLabels", verbose = FALSE)
    CheckDataFrame(data = object@meta.data$ConditionsLabels)
    IsIn_tTEscanR_Object(object = object, slot = "meta.data", section = "CorrectionFactor", verbose = FALSE)
    if (!(object@meta.data$CorrectionFactor %in% colnames(object@meta.data$ConditionsLabels))) stop("The correction factor was not found in the metadata.")
    if (verbose) message("- The ConditionsLabels and CorrectionFactor metadata have been properly loaded.")
  }

  # Loading the codon frequency and assess consistency in gene annotation mRNA_data and codon_freq
  if (!(is.null(codon_freq)) && !(class(codon_freq) %in% c("dgCMatrix", "data.frame", "matrix"))) stop("Wrong codon_freq format.\n Supported formats: dgCMatrix, data.frame, matrix")
  codon_frequency_per_gene_table <- suppressMessages(ConsistencyWithCodonFreq(data = object@assays$mRNA, codon_freq = codon_freq, species = species, filter = filter, verbose = FALSE))
  if (verbose) message("- The codon frequency per gene table has been properly loaded.")

  message("  1 . COMPLETED\n", "2 . Retrieving common mRNAs (gene expression and codon frequency table).")
  mRNAs_in_common <- intersect(colnames(codon_frequency_per_gene_table), rownames(object@assays$mRNA))
  if (is.null(mRNAs_in_common)) stop("No mRNAs in common found between the mRNA data (object@assays$mRNA) and the codon frequency table.")

  if (verbose) message("- Adding the mRNAsInCommon to the meta.data of the tTEscanR object.")
  object <- suppressMessages(Update_tTEscanR_Object(object = object, meta.data = list(mRNAs_in_common), meta.data.ids = list("mRNAsInCommon"), overwrite.assay = overwrite.assay, overwrite.metadata = overwrite.metadata, verbose = FALSE))

  if (verbose) message("- Filtering the datasets based on the common mRNAs.") # Used for computation but not saved
  filtered_mRNA_gene_expression <- as.matrix(object@assays$mRNA[object@meta.data$mRNAsInCommon, ])
  filtered_codon_frequency_per_gene_table <- as.matrix(codon_frequency_per_gene_table[ , object@meta.data$mRNAsInCommon])
  CheckDataFrame(data = filtered_mRNA_gene_expression)
  CheckDataFrame(data = filtered_codon_frequency_per_gene_table)

  message("  2 . COMPLETED\n", "3 . Computing the codon usage matrix.")
  codon_usage_per_condition <- filtered_codon_frequency_per_gene_table %*% filtered_mRNA_gene_expression
  CheckDataFrame(data = as.matrix(codon_usage_per_condition))

  if (verbose) message("- Evaluating the values in the codon usage matrix.")
  if (any(codon_usage_per_condition > .Machine$integer.max)){ # It can happen that the computed numbers are larger than those accepted by R
    codon_usage_per_condition <- codon_usage_per_condition / reduce # The reduce parameter will be used to divide all the values in the matrix.
    codon_usage_per_condition <- round(codon_usage_per_condition)
    CheckDataFrame(data = codon_usage_per_condition)

    message("- Some values of the computed matrix exceed the maximum value accepted by R.\n", paste("- The matrix has been divided by the factor, `reduce` = ", reduce,"."))
  }
  message("  3 . COMPLETED")
  count <- 4

  if (isTRUE(additional.metrics)){ # There are 3 additional metrics: (i) codon exonic background, (ii) mean codon usage, and (iii) correlation between (i) and (ii).
    message("4 . Computing the additional metrics.")
    additional.metrics <- ComputeAdditionalMetrics_CodonUsage(codon_usage = codon_usage_per_condition, codon_freq = filtered_codon_frequency_per_gene_table,
                                                              metadata = object@meta.data$ConditionsLabels, corr_method = corr_method,
                                                              corr_factor = object@meta.data$CorrectionFactor, verbose = verbose)

    if (verbose) message("- Adding the CodonUsage_AdditionalMetrics to the meta.data of the tTEscanR object.") # Nested list
    object <- suppressMessages(Update_tTEscanR_Object(object = object, main_name = "CodonUsage_AdditionalMetrics", meta.data = additional.metrics, overwrite.metadata = overwrite.metadata, verbose = FALSE))
    message("  4 . COMPLETED")
    count <- 5
  }

  message(paste(as.character(count), ". Updating the tTEscanR object."))
  object <- suppressMessages(Update_tTEscanR_Object(object = object, counts = codon_usage_per_condition, assay = "CodonUsage",
                                                    overwrite.assay = overwrite.assay, verbose = FALSE))
  message(paste(" ", as.character(count), ". COMPLETED"))
  return(object) # The output tTEscanR object has been validated in Update_tTEscanR_Object()
}

#' Compute Anticodon Usage from tRNA Gene Expression Data
#' @description
#' This function calculates **anticodon usage profiles** from tRNA gene expression data stored in a \code{tTEscanR_Object}.
#' It summarizes the expression of tRNAs by their anticodon identity, which can be used to estimate the tRNA supply landscape.
#' The tRNA gene names need to be properly annotated for proper recognition. Expected format: tRNA-Asn-GTT-5-1.
#'
#' @param object A \code{tTEscanR_Object} containing a tRNA assay.
#' @param overwrite.assay Logical; if \code{TRUE}, overwrites any existing assay in the \code{object}. Defaults to \code{FALSE}.
#' @param overwrite.metadata Logical; if \code{TRUE}, overwrites any existing metadata in the \code{object}. Defaults to \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return An updated \code{tTEscanR_Object} containing a new layer of information \code{"AnticodonUsage"} in the \code{assays} slot representing the anticodon usage.
#' @export
#'
#' @examples
#' data(subset_tRNA_data)
#' tTEscanR_obj <- Create_tTEscanR_Object(counts = subset_tRNA_data, assay = "tRNA")
#' tTEscanR_obj <- ComputeAnticodonUsage(object = tTEscanR_obj)

ComputeAnticodonUsage <- function(object, overwrite.assay = FALSE, overwrite.metadata = FALSE, verbose = TRUE){

  ###
  # CALL: User or Compute_tTE()
  # DESCRIPTION: This function takes a tRNA gene expression matrix and extracts from the tRNA genes the anticodons they code for and groups the data accordingly.
  ###

  message("1 . Checking the format of the input data.")
  if (!(inherits(object, 'tTEscanR_Object'))) stop("`object` must be a tTEscanR object.")
  if (verbose) message("- The input contains a proper tTEscanR object.")
  IsIn_tTEscanR_Object(object = object, slot = "assays", section = "tRNA", verbose = FALSE)
  CheckDataFrame(data = object@assays$tRNA)
  if (verbose) message("- The tRNA assay has been properly loaded.")
  CheckNames_tRNA(data = object@assays$tRNA)
  if (verbose) message("- The tRNA data contains tRNA genes with a suitable format.")

  message("  1 . COMPLETED\n", "2 . Extracting anticodons for each tRNA gene.")
  anticodons_for_tRNA_genes <- sapply(strsplit(rownames(object@assays$tRNA), "-"), "[[", 3)
  unique_anticodons <- sort(unique(anticodons_for_tRNA_genes))

  message("  2 . COMPLETED\n", "3 . Pooling counts from each tRNA gene with common anticodons.")
  # Filling the matrix - sum those tRNA genes with the same anticodon
  anticodon_usage_per_condition <- matrix(data = 0, nrow = length(unique_anticodons), ncol = length(colnames(object@assays$tRNA))) # Create an empty matrix (all 0)
  dimnames(anticodon_usage_per_condition) <- list(unique_anticodons, colnames(object@assays$tRNA)) # Dimensions: anticodons retrieved above (rows) and the tRNA data conditions (columns)
  CheckDataFrame(data = anticodon_usage_per_condition)

  for (i in 1:length(unique_anticodons)){ # Iterate over each anticodon retrieved from the tRNA genes
    anticodon_usage_per_condition[i, ] <- Matrix::colSums(object@assays$tRNA[which(anticodons_for_tRNA_genes == unique_anticodons[i]), , drop = FALSE])
  }

  CheckDataFrame(data = anticodon_usage_per_condition)
  message("  3 . COMPLETED\n", "4 . Updating the tTEscanR object.")
  object <- suppressMessages(Update_tTEscanR_Object(object = object, counts =  anticodon_usage_per_condition, assay = "AnticodonUsage", meta.data = anticodons_for_tRNA_genes,
                                                    meta.data.ids = "tRNAsAnticodons", overwrite.assay = overwrite.assay, overwrite.metadata = overwrite.metadata, verbose = FALSE))
  message("  4 . COMPLETED")
  return(object) # The output tTEscanR object has been validated in Update_tTEscanR_Object()
}
