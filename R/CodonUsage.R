#' Compute Codon Usage from mRNA Gene Expression Data
#' @description
#' This function estimates **codon usage profiles** based on gene-level mRNA expression data stored in a \code{tTEscanR_Object}.
#' It optionally accepts pre-computed codon frequency tables or uses internally generated default tables when not provided.
#' When enabled, it can evaluate the correlation between background codon composition and observed mean codon usage.
#' If the additional metrics are to be computed the input \code{tTEscanR_Object} needs to have a **"CorrectionFactor"** stored in the \code{"meta.data"} slot.
#' The default \code{codon_freq} were built using the canonical filter to select one transcript if several were available for the same gene.
#'
#' @param object A \code{tTEscanR_Object} containing a mRNA assay.
#' @param codon_freq Optional; a user-provided codon frequency-per-gene table. If necessary, it can be computed using \code{\link{GetCodonFreq}}.
#' @param species Optional; either \code{"hg38"} (human) or \code{"mm39"} (mouse) to load the default settings. Required if \code{codon_freq} is not provided.
#' @param reduce Numeric; a scaling factor used to normalize large expression values that exceed R's handling capacity. Defaults to 100.
#' @param additional_metrics Logical; if \code{TRUE}, computes: (i) codon exonic background, (ii) mean codon usage, and (iii) correlation between the previous metrics. Defaults to \code{TRUE}.
#' @param corr_method A correlation method accepted by \code{\link{cor}}. Required if \code{additional_metrics} is \code{TRUE}. Defaults to \code{"spearman"}.
#' @param overwrite Logical; if \code{TRUE}, overwrites any existing assay and metadata in the \code{object}. Defaults to \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return An updated \code{tTEscanR_Object} containing a new layer of information \code{"CodonUsage"} in the \code{assays} slot representing the codon usage.
#' Additional computations will be stored in the \code{meta.data} slot as \code{"CodonUsage_AdditionalMetrics"}.
#' @export
#'
#' @examples
#' data(default_tTEscanR_mRNA_data, default_tTEscanR_metadata)
#' tTEscanR_obj <- Create_tTEscanR_Object(counts = default_tTEscanR_mRNA_data, assay = "mRNA",
#'                                        meta.data = list(default_tTEscanR_metadata, "tissue"),
#'                                        meta.data.ids = list("ConditionsLabels", "CorrectionFactor"))
#' tTEscanR_obj <- ComputeCodonUsage(object = tTEscanR_obj, species = "hg38",
#'                                   additional_metrics = FALSE, reduce = 1000)

ComputeCodonUsage <- function(object, codon_freq = NULL, species = NULL, reduce = 100, additional_metrics = TRUE,
                              corr_method = "spearman", overwrite = FALSE, verbose = TRUE){

  ###
  # CALL: User or Compute_tTE()
  # DESCRIPTION: This function takes a mRNA gene expression matrix and a codon frequency per gene table to compute the codon usage.
  # The codon usage represents how much of each codon (row) is needed for translation in a given condition (column).
  ###

  message("\n--- Computation of the codon usage ---", "\n1 . Checking the format of the input data.")
  if (!inherits(object, 'tTEscanR_Object')) stop("'object' must be a tTEscanR object.")
  if (verbose) message("- The input consists of a proper tTEscanR object.")

  # Checking that the mRNA data is in the object and has a proper format
  IsIn_tTEscanR_Object(object = object, slot = "assays", section = "mRNA", verbose = FALSE)
  mRNA_data <- object@assays$mRNA
  CheckDataFrame(data = mRNA_data)
  if (verbose) message("- The mRNA assay has been properly loaded.")

  if (additional_metrics){
    IsIn_tTEscanR_Object(object = object, slot = "meta.data", section = "ConditionsLabels", verbose = FALSE)
    metadata_table <- object@meta.data[["ConditionsLabels"]]
    CheckDataFrame(data = metadata_table)

    IsIn_tTEscanR_Object(object = object, slot = "meta.data", section = "CorrectionFactor", verbose = FALSE)
    batch <- object@meta.data[["CorrectionFactor"]]
    if (!(batch %in% colnames(metadata_table))) stop("The correction factor was not found in the metadata.")
    if (verbose) message("- The 'ConditionsLabels' and 'CorrectionFactor' have been properly loaded.")
  }

  # Loading the codon frequency and assess consistency in gene annotation mRNA_data and codon_freq
  if (!is.null(codon_freq)) {
    if (!any(class(codon_freq) %in% c("dgCMatrix", "data.frame", "matrix"))) {
      stop("Wrong 'codon_freq' format.\n Supported formats: dgCMatrix, data.frame, matrix")
    }
  }
  # if (!is.null(codon_freq) && !(class(codon_freq) %in% c("dgCMatrix", "data.frame", "matrix"))) stop("Wrong 'codon_freq' format.\n Supported formats: dgCMatrix, data.frame, matrix")
  codon_frequency_per_gene_table <- suppressMessages(ConsistencyWithCodonFreq(data = mRNA_data, codon_freq = codon_freq, species = species, verbose = FALSE))
  if (verbose) message("- The codon frequency per gene table has been properly loaded.")

  message("1 . COMPLETED\n", "2 . Retrieving common mRNAs (gene expression and codon frequency table).")
  mRNAs_in_common <- intersect(colnames(codon_frequency_per_gene_table), rownames(mRNA_data))
  if (length(mRNAs_in_common) == 0) stop("No mRNAs in common found between the mRNA data and the codon frequency table.")

  message("2 . COMPLETED\n", "3 . Computing the codon usage matrix.")
  if (verbose) message("- Filtering the datasets based on the common mRNAs.") # Used for computation but not saved
  filtered_codon_frequency_per_gene_table <- as.matrix(codon_frequency_per_gene_table[, mRNAs_in_common])
  codon_usage_per_condition <- filtered_codon_frequency_per_gene_table %*% as.matrix(mRNA_data[mRNAs_in_common, ])
  CheckDataFrame(data = as.matrix(codon_usage_per_condition))

  if ( max(codon_usage_per_condition) >.Machine$integer.max) { # It can happen that the computed numbers are larger than those accepted by R
    if (verbose) message("- Scaling down the matrix to prevent integer overflow.")
    codon_usage_per_condition <- round(codon_usage_per_condition / reduce) # The reduce parameter will be used to divide all the values in the matrix.
    CheckDataFrame(data = codon_usage_per_condition)

    message("- Some values of the computed matrix exceed the maximum value accepted by R.\n", paste("- The matrix has been divided by the factor, 'reduce = ", reduce,"'."))
  }
  message("3 . COMPLETED")
  count <- 4

  if (isTRUE(additional_metrics)){ # There are 3 additional metrics: (i) codon exonic background, (ii) mean codon usage, and (iii) correlation between (i) and (ii).
    message("4 . Computing the additional metrics.")
    available_id_col <- IsIn_tTEscanR_Object(object = object, slot = "meta.data", section = "DataMetadataIndex", verbose = FALSE, compute_assay = TRUE)
    if (isFALSE(available_id_col)) id_col <- NULL else id_col <- object@meta.data[["DataMetadataIndex"]]

    additional_metrics <- ComputeAdditionalMetrics_CodonUsage(codon_usage = codon_usage_per_condition, codon_freq = filtered_codon_frequency_per_gene_table,
                                                              metadata = metadata_table, id_col = id_col, corr_method = corr_method, batch = batch, verbose = verbose)

    if (verbose) message("- Adding the CodonUsage_AdditionalMetrics to the meta.data of the tTEscanR object.") # Nested list
    object <- suppressMessages(Update_tTEscanR_Object(object = object, main_name = "CodonUsage_AdditionalMetrics", meta.data = additional_metrics, overwrite = overwrite, verbose = FALSE))
    message("4 . COMPLETED")
    count <- 5
  }

  message(paste(as.character(count), ". Updating the tTEscanR object."))
  object <- suppressMessages(Update_tTEscanR_Object(object = object, counts = codon_usage_per_condition, assay = "CodonUsage",
                                                    meta.data = list(mRNAs_in_common), meta.data.ids = list("mRNAsInCommon"),
                                                    overwrite = overwrite, verbose = FALSE))
  message(paste(as.character(count), ". COMPLETED\n"), "--- The codon usage has been successfully computed ---\n")
  return(object) # The output tTEscanR object has been validated in Update_tTEscanR_Object()
}

#' Compute Anticodon Usage from tRNA Gene Expression Data
#' @description
#' This function calculates **anticodon usage profiles** from tRNA gene expression data stored in a \code{tTEscanR_Object}.
#' It summarizes the expression of tRNAs by their anticodon identity, which can be used to estimate the tRNA supply landscape.
#' The tRNA gene names need to be properly annotated for proper recognition. Expected format: tRNA-Asn-GTT-5-1.
#'
#' @param object A \code{tTEscanR_Object} containing a tRNA assay.
#' @param overwrite Logical; if \code{TRUE}, overwrites any existing assay and metadata in the \code{object}. Defaults to \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return An updated \code{tTEscanR_Object} containing a new layer of information \code{"AnticodonUsage"} in the \code{assays} slot representing the anticodon usage.
#' @export
#'
#' @examples
#' data(default_tTEscanR_tRNA_data)
#' tTEscanR_obj <- Create_tTEscanR_Object(counts = default_tTEscanR_tRNA_data, assay = "tRNA")
#' tTEscanR_obj <- ComputeAnticodonUsage(object = tTEscanR_obj)

ComputeAnticodonUsage <- function(object, overwrite = FALSE, verbose = TRUE){

  ###
  # CALL: User or Compute_tTE()
  # DESCRIPTION: This function takes a tRNA gene expression matrix and extracts from the tRNA genes the anticodons they code for and groups the data accordingly.
  ###

  message("\n--- Computation of the anticodon usage ---", "\n1 . Checking the format of the input data.")
  if (!(inherits(object, 'tTEscanR_Object'))) stop("'object' must be a tTEscanR object.")
  if (verbose) message("- The input contains a proper tTEscanR object.")
  IsIn_tTEscanR_Object(object = object, slot = "assays", section = "tRNA", verbose = FALSE)
  tRNA_data <- object@assays$tRNA
  tRNA_data <- as.matrix(tRNA_data)
  CheckDataFrame(data = tRNA_data)

  tRNA_genes <- rownames(tRNA_data)
  CheckNames_tRNA(gene_names = tRNA_genes)
  if (verbose) message("- The tRNA assay has been properly loaded and contains tRNA genes with a suitable format.")

  message("1 . COMPLETED\n", "2 . Extracting the anticodons of each tRNA gene.")
  anticodons_for_tRNA_genes <- sapply(strsplit(tRNA_genes, "-"), "[[", 3)
  unique_anticodons <- sort(unique(anticodons_for_tRNA_genes))

  message("2 . COMPLETED\n", "3 . Pooling the counts from each tRNA gene with common anticodons.")
  # Initialize an mapping matrix - 1 to add the gene and 0 to ignore it
  map_factor <- factor(anticodons_for_tRNA_genes, levels = unique_anticodons)
  M <- Matrix:: sparse.model.matrix(~ 0 + map_factor)
  M <- Matrix:: t(M) # Set the anticodons as rows
  rownames(M) <- gsub("map_factor", "", rownames(M)) # Clean the rownames

  # Filling the matrix - sum those tRNA genes with the same anticodon using matrix multiplication
  anticodon_usage_per_condition <- as.matrix(M %*% tRNA_data)
  rownames(anticodon_usage_per_condition) <- unique_anticodons
  colnames(anticodon_usage_per_condition) <- colnames(tRNA_data)
  CheckDataFrame(data = anticodon_usage_per_condition) # Dimensions: anticodons retrieved above (rows) and the tRNA data conditions (columns)

  message("3 . COMPLETED\n", "4 . Updating the tTEscanR object.")
  object <- suppressMessages(Update_tTEscanR_Object(object = object, counts =  anticodon_usage_per_condition, assay = "AnticodonUsage", meta.data = anticodons_for_tRNA_genes,
                                                    meta.data.ids = "tRNAsAnticodons", overwrite = overwrite, verbose = FALSE))
  message("4 . COMPLETED\n", "--- The anticodon usage has been successfully computed ---\n")
  return(object) # The output tTEscanR object has been validated in Update_tTEscanR_Object()
}
