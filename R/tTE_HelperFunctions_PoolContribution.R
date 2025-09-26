ComputeIndividualGeneCorrelation <- function(codon_usage, mean_codon_usage, corr_method){

  ###
  # CALL: ExaminePoolContribution()
  # DESCRIPTION: Identification of outliers: calculate correlation of each condition's codon usage to the mean codon usage across conditions.
  # This function takes as input the size corrected codon usage matrix and the mean codon usage.
  # It correlates each individual codon usage value with the mean codon usage (reference).
  ###

  # Generates a vector with the conditions included in data
  condition_correlations_to_mean_codon_usage <- numeric(length(colnames(codon_usage)))
  names(condition_correlations_to_mean_codon_usage) <- colnames(codon_usage)

  # Normalize the data and convert it into long format
  codon_usage_size_corrected <- DataToLongFormat(data = codon_usage, normalize = TRUE, rownames_to_column = "codon", names_to = "condition", values_to = "usage")

  for (i in 1:length(colnames(codon_usage))) { # Iterates over the conditions in data
    # Retrieve the codon usage of the condition and correlate the mean codon usage and the current codon usage
    select_condition_codon_usage <- dplyr::filter(codon_usage_size_corrected, .data$condition == colnames(codon_usage)[i])$usage
    condition_correlations_to_mean_codon_usage[i] <- stats::cor(mean_codon_usage$mean_usage_across_conditions, select_condition_codon_usage, method = corr_method)
  }
  return(condition_correlations_to_mean_codon_usage) # Returns a matrix that stores how similar the codon usage in a particular condition is to the mean codon usage across all conditions
}

ComputeTopNGenes <- function(data, N){

  ###
  # CALL: ExaminePoolContribution()
  # DESCRIPTION: This function takes each condition of the data and sorts the codons by their contribution and extracts the top N contributors.
  ###

  # Define an empty matrix based on the conditions in the codon_pool_contribution
  codon_pool_contribution_of_top_N_genes <- numeric(length(colnames(data)))
  top_N_codon_pool_contributors_per_condition <- matrix(data = 0, nrow = N, ncol = length(colnames(data)))

  for (i in 1:length(colnames(data))) { # Iterate over each condition in the input data
    top_N_genes <- sort(data[, i], decreasing = TRUE)[1:N] # Retrieves the top genes by sorting them based on their contribution
    codon_pool_contribution_of_top_N_genes[i] <- sum(top_N_genes) # Sum the contribution of the top N genes
    top_N_codon_pool_contributors_per_condition[, i] <- names(top_N_genes)
  }

  colnames(top_N_codon_pool_contributors_per_condition) <- colnames(data)
  rownames(top_N_codon_pool_contributors_per_condition) <- paste("top", seq_len(nrow(top_N_codon_pool_contributors_per_condition)), "gene", sep = "")

  codon_pool_diversity_vs_correlation <- tibble::tibble(condition = colnames(data), codon_diversity = 1 - codon_pool_contribution_of_top_N_genes)
  return(list(codon_pool_diversity_vs_correlation, top_N_codon_pool_contributors_per_condition)) # Returns a list with two matrices
}

ComputeWithoutTopNGenes <- function(data, codon_freq, mean_codon_usage, extract_topN_genes, corr_method){

  ###
  # CALL: ExaminePoolContribution()
  # DESCRIPTION: This function computes the correlation between the codon usage matrix computed
  # without the top N gene codon pool contributors and the mean codon usage across conditions (baseline, considering all genes).
  ###

  for (i in 1:length(colnames(data))) { # Remove expression of top N codon pool contributors (so their codon pool contribution becomes 0)
    select_condition_top_N_contributors <- extract_topN_genes[[2]][, i]
    data[select_condition_top_N_contributors, ] <- 0
  }

  # Computes the codon usage without the top N codon pool contributors
  removed_top_N_genes_codon_usage <- as.matrix(codon_freq) %*% as.matrix(data)
  removed_top_N_genes_codon_usage <- sweep(removed_top_N_genes_codon_usage, 2, colSums(removed_top_N_genes_codon_usage), "/") # Normalize data

  # Define an empty matrix based on the conditions in the data
  removed_top_N_genes_correlation_to_mean_codon_usage <- numeric(length(colnames(data)))
  names(removed_top_N_genes_correlation_to_mean_codon_usage) <- colnames(data)

  for (i in 1:length(colnames(data))) { # Fills the matrix with the correlation values
    removed_top_N_genes_correlation_to_mean_codon_usage[i] <- stats::cor(removed_top_N_genes_codon_usage[, i], mean_codon_usage$mean_usage_across_conditions, method = corr_method)
  }

  codon_pool_diversity_vs_correlation <- tibble::tibble(condition = colnames(data),
                                                        codon_diversity = 1 - extract_topN_genes[[1]]$codon_diversity,
                                                        correlation = removed_top_N_genes_correlation_to_mean_codon_usage)
  return(codon_pool_diversity_vs_correlation) # Returns the correlation matrix
}
