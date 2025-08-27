#' Run Codon Usage Permutation Test (Background)
#' @description
#' This function performs a **permutation test** to compare a mRNA dataset to the reference codon frequency-per-gene matrix.
#'
#' @param n_permutations Numeric; number of permutations to perform. Defaults to 1000.
#' @param n_features Numeric; number of features to select in each permutation. Defaults to 100. If \code{target_data} is given this parameter will take its length.
#' @param target_data Optional; a mRNA expression count matrix with features as rows and conditions as columns.
#' @param codon_freq Optional; a user-provided codon frequency per gene table. If necessary, it can be computed using \code{\link{ObtainCodonFreqPerGene}}.
#' @param species Optional, a \code{character} string specifying the species reference genome version (used if \code{codon_freq} is not provided or \code{translate} is \code{TRUE}). Supported values include \code{"hg38"} (human) and \code{"mm39"} (mouse).
#' @param filter Optional; a \code{character} string specifying how to choose among multiple transcripts per gene either \code{"canonical"} (default) or \code{"length"} (longest transcript).
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return A table with the codons and their frequencies after computing all the permutations.
#' @export
#'
#' @examples
#' data(subset_mRNA_data)
#'
#' genes <- subset_mRNA_data[1:20, ]
#' permut <- GetPermutationDist(n_permutations = 100, target_data = genes, species = "hg38")

GetPermutationDist <- function(n_permutations = 1000, n_features = 100, target_data = NULL, codon_freq = NULL, species = NULL, filter = "canonical", verbose = TRUE){

  ###
  # CALL: User
  # DESCRIPTION: This function performs a permutation test with respect to the reference codon_freq genes.
  # In case of having a target_data input, those features are removed from the data used in the permutations.
  ###

  message("1 . Loading the codon frequency per gene table.")
  if (!is.null(target_data)){ # TARGETED APPROACH
    codon_frequency_per_gene_table <- ConsistencyWithCodonFreq(data = target_data, codon_freq = codon_freq, species = species, filter = filter, verbose = verbose)

    target_indexes <- which(colnames(codon_frequency_per_gene_table) %in% rownames(target_data))
    if (length(target_indexes) == 0) stop("None of the genes in `target_data` are present in the codon frequency table.")

    codon_frequency_per_gene_table <- codon_frequency_per_gene_table[, -target_indexes]
    n_features <- length(target_indexes)
  } else { # CONTROL APPROACH - ONLY TAKES THE codon_freq
    codon_frequency_per_gene_table <- CheckCodonFreqTable(data = codon_freq, species = species, filter = filter, verbose = verbose)
  }
  message("  1 . COMPLETED\n", "2 . Starting permutation test.")

  perm_matrices <- data.frame(codon = character(), freq = numeric(), stringsAsFactors = FALSE) # Store permutation data into a dataframe

  for (i in 1:n_permutations){
    if (verbose) message(paste("- Permutation", i))

    # Random subsampling of features
    subset_labels <- sample(colnames(codon_frequency_per_gene_table), n_features) # without duplicates - COLUMNS
    codon_usage_perm <- codon_frequency_per_gene_table[ , subset_labels]
    # Relative contribution of each row to the overall total sum of all values
    codon_usage_perm_relative_cont <- rowSums(codon_usage_perm) / sum(rowSums(codon_usage_perm))
    # Organize the data into a dataframe
    codon_usage_perm_relative_cont <- data.frame(codon = names(codon_usage_perm_relative_cont), freq = as.numeric(codon_usage_perm_relative_cont), row.names = NULL)
    perm_matrices <- rbind(perm_matrices, codon_usage_perm_relative_cont)
  }

  perm_matrices <- perm_matrices %>% dplyr::arrange(.data$codon, .data$freq)
  message("  2 . COMPLETED")
  return(perm_matrices) # Returns permutation matrix
}

#' Assess Significance (P-value) & Corrects for Multiple Hypothesis Testing
#'
#' @param dist A table with the codons and their frequencies after completing a permutation test. Output from \code{\link{GetPermutationDist}}.
#' @param value A \code{list} of \code{data.frame} of the codon exonic background of a mRNA gene expression matrix. The codon exonic background can be computed in \code{\link{ComputeCodonUsage}}.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return A table with the codon exonic background and their significance level before (p-value) and after the correction (p-adjusted value).
#' @export
#'
#' @examples
#' data(subset_mRNA_data, metadata)
#' selected_genes <- subset_mRNA_data[1:20, ]
#' permutation_test <- GetPermutationDist(n_permutations = 100, target_data = selected_genes,
#'                                        species = "hg38")
#' tTEscanR_obj <- Create_tTEscanR_Object(counts = subset_mRNA_data, assay = "mRNA",
#'                                        meta.data = list(metadata, "tissue"),
#'                                        meta.data.ids = list("ConditionsLabels", "CorrectionFactor"))
#' tTEscanR_obj <- ComputeCodonUsage(object = tTEscanR_obj, species = "hg38",
#'                                   additional.metrics = FALSE, reduce = 1000)
#'
#' codon_usage <- tTEscanR_obj@assays$CodonUsage
#' codon_background <- rowSums(codon_usage) / sum(rowSums(codon_usage))
#' codons_to_AA <- FeaturesToAA(data_to_translate = names(codon_background),
#'                              notation.from = "codon", notation.to = "aa")
#' codon_background <- data.frame(group = codons_to_AA, codon = names(codon_background),
#'                                freq = as.numeric(codon_background), row.names = NULL)
#' significance <- ObtainSignificance(dist = permutation_test, value = codon_background)

ObtainSignificance <- function(dist, value, verbose = TRUE){

  message("1 . Checking the input data.")
  if (verbose) message("- Analyzing the content in `value`.")
  if (ncol(value) != 3 || (!(is.character(value[[1]])) && !(is.character(value[[2]])) && !(is.numeric(value[[3]])))){
    stop("Wrong format of the `value` parameter.\n", "Please, make sure that it contains 3 columns (strict order).\n",
         "(i) group of the features - character\n", "(ii) features - character\n", "(iii) codon background frequencies - numeric")
  }

  if (verbose) message("- Analyzing the content in `dist`.")
  if (ncol(dist) != 2 || (!(is.character(dist[[1]])) && !(is.numeric(dist[[2]])))){
    stop("Wrong format of the `dist` parameter.\n", "Please, make sure that it contains 2 columns (strict order).\n",
         "(i) features\n", "(ii) codon background frequencies")
  }
  message("  1 . COMPLETED\n", "2 . Iterate over each group.")

  groups <- unique(value[[1]]) # Extract the unique conditions in the data
  p_val <- sig <- test <- list() # Create empty lists to store the results

  for (g in groups) { # Iterate over each condition in the data
    group_data <- value[value[[1]] == g, ] # Get a data subset - based on the group selected in each iteration

    for (c in unique(group_data[[2]])) { # Iterate over the codons present in the data
      observed_row <- group_data[group_data[[2]] == c, ]
      distribution_rows <- dist[dist[[1]] == c, ]

      if (nrow(observed_row) > 0 && nrow(distribution_rows) > 0) {
        obs_freq <- observed_row[[3]]
        dist_freq <- distribution_rows[[2]]
        if (obs_freq < stats::median(dist_freq)){
          result <- length(dist_freq[dist_freq <= obs_freq])/length(dist_freq)
          test <- "left"
        } else {
          result <- length(dist_freq[dist_freq >= obs_freq])/length(dist_freq)
          test <- "right"
        }

        # Store the results
        p_val[[paste(g, c, sep = "-")]] <- result
        test[[paste(g, c, sep = "-")]] <- test
        sig[[paste(g, c, sep = "-")]] <- (result < 0.05)
      }
    }
  }

  if(is.null(p_val) || is.null(sig)) stop("The significance test could not be executed.\n", "Please revise the input parameters.")

  message("  2 . COMPLETED\n", "3 . Performing FDR multiple test correction.")
  results <- data.frame(name = names(p_val), p_val = unlist(p_val), sig = unlist(sig), row.names = NULL) # Define a data frame with the results of all the iterations
  CheckDataFrame(results) # Add checkpoint
  results <- tidyr::separate(results, .data$name, into = c("group", "codon"), sep = "-")

  # Perform FDR correction for multiple testing
  results <- results %>% dplyr::group_by(.data$codon) %>% dplyr::mutate(p_val_adj = stats::p.adjust(p_val, method = "BH")) %>% dplyr::ungroup()

  # Transform the format of the data
  results$p_val_adj <- round(results$p_val_adj, 4)
  results$sig_adj <- (results$p_val_adj < 0.05)
  value_results <- value %>% dplyr::left_join(results, by = c("group", "codon"))
  message("  3 . COMPLETED")
  return(value_results)
}
