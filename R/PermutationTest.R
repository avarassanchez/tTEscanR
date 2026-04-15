#' Run Codon Usage Permutation Test (Background)
#' @description
#' This function performs a **permutation test** to compare a mRNA dataset to the reference codon frequency-per-gene matrix.
#'
#' @param n_permut Numeric; number of permutations to perform. Defaults to 1000.
#' @param n_features Numeric; number of features to select in each permutation. Defaults to 100. If \code{target_data} is given this parameter will take its length.
#' @param target_data Optional; a mRNA expression count matrix with features as rows and conditions as columns.
#' @param codon_freq Optional; a user-provided codon frequency per gene table. If necessary, it can be computed using \code{\link{GetCodonFreq}}.
#' @param species Optional, a \code{character} string specifying the species reference genome version (used if \code{codon_freq} is not provided or \code{translate} is \code{TRUE}). Supported values include \code{"hg38"} (human) and \code{"mm39"} (mouse).
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return A table with the codons and their frequencies after computing all the permutations.
#' @export
#'
#' @examples
#' data(default_tTEscanR_mRNA_data)
#'
#' genes <- default_tTEscanR_mRNA_data[1:20, ]
#' permut <- GetPermutationDist(n_permut = 100, target_data = genes, species = "hg38")

GetPermutationDist <- function(n_permut = 1000, n_features = 100, target_data = NULL, codon_freq = NULL, species = NULL, verbose = TRUE){

  ###
  # CALL: User
  # DESCRIPTION: This function performs a permutation test with respect to the reference codon_freq genes.
  # In case of having a target_data input, those features are removed from the data used in the permutations.
  ###

  if (n_permut < 2) stop("Error: 'n_permut' must be at least 2.")

  message("--- Computing the permutation analysis ---", "\n1 . Loading the codon frequency per gene table.")
  if (!is.null(target_data)){ # TARGETED APPROACH
    codon_frequency_per_gene_table <- ConsistencyWithCodonFreq(data = target_data, codon_freq = codon_freq, species = species, verbose = verbose)

    target_names <- rownames(target_data)
    is_target <- colnames(codon_frequency_per_gene_table) %in% target_names
    if (!any(is_target)) stop("Error: None of the genes in 'target_data' are present in the codon frequency table.")

    n_features <- sum(is_target)
    codon_frequency_per_gene_table <- codon_frequency_per_gene_table[, !is_target, drop = FALSE] # Isolate the background by removing the target genes

  } else { # CONTROL APPROACH - ONLY TAKES THE codon_freq
    codon_frequency_per_gene_table <- CheckCodonFreqTable(data = codon_freq, species = species, verbose = verbose)
  }

  codon_frequency_per_gene_table <- as.matrix(codon_frequency_per_gene_table)
  n_total_genes <- ncol(codon_frequency_per_gene_table)
  message("1 . COMPLETED\n", "2 . Starting permutation test...")

  perm_results <- replicate(n_permut, {
    index <- sample.int(n_total_genes, n_features) # Sampling of indices
    sub_sums <- rowSums(codon_frequency_per_gene_table[, index, drop = FALSE]) # Subset and sum rows
    return(sub_sums / sum(sub_sums)) # Normalized relative contribution
  })

  codons <- rownames(codon_frequency_per_gene_table)
  perm_matrices <- data.frame(codon = rep(codons, times = n_permut), freq = as.numeric(perm_results), stringsAsFactors = FALSE)
  perm_matrices <- perm_matrices[order(perm_matrices$codon, perm_matrices$freq), ]
  rownames(perm_matrices) <- NULL

  message("2 . COMPLETED\n", "--- The permutation analysis has been successfully executed ---")
  return(perm_matrices) # Returns permutation matrix
}

#' Assess Significance (P-value) & Corrects for Multiple Hypothesis Testing
#'
#' @param dist A table with the codons and their frequencies after completing a permutation test. Output from \code{\link{GetPermutationDist}}.
#' @param value A \code{list} of \code{data.frame} of the codon exonic background of a mRNA gene expression matrix. The codon exonic background can be computed in \code{\link{ComputeCodonUsage}}.
#' @param padj_threshold Numeric; p-value threshold used for highlighting significant features in the volcano plot. Defaults to 0.05.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return A table with the codon exonic background and their significance level before (p-value) and after the correction (p-adjusted value).
#' @export
#'
#' @examples
#' data(default_tTEscanR_mRNA_data, default_tTEscanR_metadata)
#' selected_genes <- default_tTEscanR_mRNA_data[1:20, ]
#' permutation_test <- GetPermutationDist(n_permut = 100, target_data = selected_genes,
#'                                        species = "hg38")
#' tTEscanR_obj <- Create_tTEscanR_Object(counts = default_tTEscanR_mRNA_data, assay = "mRNA",
#'                                        meta.data = list(default_tTEscanR_metadata, "tissue"),
#'                                        meta.data.ids = list("ConditionsLabels", "CorrectionFactor"))
#' tTEscanR_obj <- ComputeCodonUsage(object = tTEscanR_obj, species = "hg38",
#'                                   additional_metrics = FALSE, reduce = 1000)
#'
#' codon_usage <- tTEscanR_obj@assays$CodonUsage
#' codon_background <- rowSums(codon_usage) / sum(rowSums(codon_usage))
#' codons_to_AA <- FeaturesToAA(data_to_translate = names(codon_background),
#'                              notation_from = "codon", notation_to = "aa")
#' codon_background <- data.frame(group = codons_to_AA, codon = names(codon_background),
#'                                freq = as.numeric(codon_background), row.names = NULL)
#' significance <- ObtainSignificance(dist = permutation_test, value = codon_background)

ObtainSignificance <- function(dist, value, padj_threshold = 0.05, verbose = TRUE){

  message("--- Computing the statistical significance ---", "\n1 . Checking the input data.")
  value <- as.data.frame(value)
  dist <- as.data.frame(dist)

  if (ncol(value) != 3 || ncol(dist) != 2) stop("Error: Incorrect number of columns.\n",
                                                "'value' must have 3 columns: (i) group, (ii) codon, (iii) frequency.\n",
                                                "'dist' must have 2 columns: (i) codon, (ii) frequency.")

  colnames(value) <- c("group", "codon", "freq")
  colnames(dist) <- c("codon", "freq")
  if (!is.numeric(value$freq) || !is.numeric(dist$freq)) stop("Error: The frequency values in both 'value' and 'dist' must be numeric characters.")

  message("1 . COMPLETED\n", "2 . Calculating empirical p-values.")
  dist_list <- split(dist$freq, dist$codon)
  dist_medians <- vapply(dist_list, stats::median, numeric(1)) # Calculate the median for all the codons to use as a reference

  results_list <- lapply(seq_len(nrow(value)), function(i){
    current_codon <- value$codon[i]
    obs_freq <- as.numeric(value$freq[i])

    d_freq <- dist_list[[current_codon]]
    if (is.null(d_freq)) return(data.frame(p_val = NA, tail = NA))

    # Compute left or right tail based on the median
    if (obs_freq < dist_medians[current_codon]) {
      p <- mean(d_freq <= obs_freq)
      t <- "left"
    } else {
      p <- mean(d_freq >= obs_freq)
      t <- "right"
    }
    return(data.frame(p_val = p, tail = t))
  })

  res_df <- do.call(rbind, results_list)
  value_results <- cbind(value, res_df)

  message("2 . COMPLETED\n", "3 . Performing FDR multiple test correction.")
  value_results$p_val_adj <- stats:: p.adjust(value_results$p_val, method = "BH")
  value_results$p_val_adj <- round(value_results$p_val_adj, 4)
  value_results$sig_adj <- value_results$p_val_adj < padj_threshold
  # value_results <- value_results %>%
  #   dplyr::group_by(.data$codon) %>%
  #   dplyr::mutate(p_val_adj = stats::p.adjust(.data$p_val, method = "BH"), p_val_adj = round(.data$p_val_adj, 4), sig_adj = .data$p_val_adj < 0.05) %>%
  #   dplyr::ungroup()

  message("3 . COMPLETED\n", "--- The statistical significnace has been successfully computed ---")
  return(value_results)
}
