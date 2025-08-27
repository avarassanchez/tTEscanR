ComputeStatisticalSignificance <- function(level, tTE_scores, data_mRNA, data_tRNA, tRNA_exp, corr_method, verbose){

  ###
  # CALL: Compute_tTE()
  # DESCRIPTION: This function takes computes a statistical test to obtain the significance of tTE scores  by comparing the mRNA data (codon usage or AA demand) to the tRNA expression data.
  # The tRNA expression data is necessary to increase the number of features to consider to perform the randomization.
  # If performed with the anticodon usage or AA supply less features could be interrogated, and therefore we will achieve less power.
  ###

  # Extract the features and the conditions to consider - `data_mRNA` and `data_tRNA` have been filtered for common conditions
  tTE_conditions <- colnames(data_tRNA)
  tTE_factors <- rownames(data_tRNA) # Anticodons
  tTE_scores <- as.numeric(tTE_scores) # Extract the numeric vector
  tRNAsAnticodons <- sapply(strsplit(rownames(tRNA_exp), "-"), "[[", 3) # Extract the anticodons from the tRNA data

  if (level == "aa") { # Extract the AA from the anticodons based on the genetic code
    level_for_tRNA_genes <- as.vector(Biostrings::translate(Biostrings::reverseComplement(Biostrings::DNAStringSet(tRNAsAnticodons)), no.init.codon = TRUE))
  } else { # Directly performing the analysis over the anticodons
    level_for_tRNA_genes <- tRNAsAnticodons
  }

  tTE_p_values <- numeric(length(tTE_conditions)) # Generate another empty vector (all values = 0), this time to store the significance results

  for (i in 1:length(tTE_conditions)) { # Iterate over the conditions present in the data
    select_condition <- tTE_conditions[i] # Retrieve the condition interrogated each time
    N = 1000 # Number of null distribution to sample

    # Select the mRNA and tRNA usage of the selected condition (vector)
    mRNA_usage <- data_mRNA[, select_condition] # observed (true) AA supply usage for this condition
    true_tRNA_usage <- tRNA_exp[, select_condition] # observed (true) tRNA gene usage for this condition: this will be shuffled N times to determine statistical significance

    null_TEs_level = numeric(N) # Generate another empty vector (all values = 0), this time to store the values of each permutation

    for (j in 1:N) { # Iterates over the name of iterations defined by parameter N
      shuffled_tRNA_usage <- sample(true_tRNA_usage) # Randomizes the order of the tRNA genes
      random_tibble <- tibble::tibble(sums = shuffled_tRNA_usage, level = level_for_tRNA_genes) # Sums the counts per each tRNA gene
      random_usage <- dplyr::summarize(dplyr::group_by(random_tibble, .data$level), sums = sum(.data$sums)) # Sums together the counts per each anticodon
      random_usage <- random_usage[which(random_usage$level %in% tTE_factors), ] # Checks which anticodons are in the data

      # In order to perform the correlation both data pieces need to have matching features
      # Filter the mRNA data (mRNA_usage) based on the features in random_usage - specific for the codon level
      if (level == "codon") {
        keep_indexes <- which(as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(names(mRNA_usage)))) %in% random_usage$level)
        mRNA_usage_filtered <- mRNA_usage[keep_indexes]
        names(mRNA_usage_filtered) <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(names(mRNA_usage_filtered))))
        random_usage <- random_usage[random_usage$level %in% names(mRNA_usage_filtered), ]
      } else {
        mRNA_usage_filtered <- mRNA_usage
      }

      # Extract the sums vector and divide it by the total usage
      random_usage <- random_usage$sums
      random_usage <- random_usage / sum(random_usage)

      # Fills the value of the matrix with the tTE score (correlation)
      null_TEs_level[j] <- stats::cor(random_usage, mRNA_usage_filtered, method = corr_method) # Correlation between shuffled tRNA data and unshuffled mRNA data (AA demand or codon usage)
    }

    # Fit parameters (mu and sigma for null distribution)
    null_TEs_level <- data.frame(null_dist = null_TEs_level)
    mean_random_correlations <- mean(null_TEs_level$null_dist)
    sd_random_correlations <- stats::sd(null_TEs_level$null_dist)

    # Calculate p-value for the actual tTE
    tTE_p_values[i] <- 1 - stats::pnorm(tTE_scores[i], mean = mean_random_correlations, sd = sd_random_correlations)
  }

  # Generate summary table
  tTE_results_tibble <- tibble::tibble(condition = tTE_conditions, tTE = tTE_scores, p_value = tTE_p_values, neg_log10_tTE_p_value = -log10(tTE_p_values))
  return(tTE_results_tibble) # The output consists of a table with the input tTE scores and their corresponding p-values
}

ComputeCorrelation <- function(data_mRNA, data_tRNA, corr_method, verbose){

  ###
  # CALL: Compute_tTE()
  # DESCRIPTION: This function takes two count matrices and performs the correlation between their entries.
  ###

  # Extract the features and the conditions to consider - both datasets have been filtered for common conditions
  tTE_conditions <- colnames(data_mRNA)
  tTE_factors <- rownames(data_mRNA)
  tTEs_at_level <- numeric(length(tTE_conditions)) # Generate empty vector with as many elements as in tTE_conditions
  names(tTEs_at_level) <- tTE_conditions # Include as names of the vector the elements in tTE_conditions

  if (verbose) message(paste("- Computing the", corr_method, "correlation."))
  for (i in 1:length(tTE_conditions)) { # Iterates over the shared conditions
    # Selects the condition that we are analyzing in each iteration and filters the data matrices accordingly
    select_condition <- tTE_conditions[i]
    values_from_tRNA <- data_tRNA[, select_condition]
    values_from_mRNA <- data_mRNA[tTE_factors, select_condition]

    tTEs_at_level[i] <- stats::cor(values_from_tRNA, values_from_mRNA, method = corr_method) # The tTE is computed as the correlation coefficient
  }

  if (is.null(tTEs_at_level)) stop("The tTE scores could not be computed.")
  return(tTEs_at_level) # The output consists of vector with the features and the correlation value, the theoretical translation efficiency (tTE) score
}

FilterMatrix <- function(data_mRNA, data_tRNA, level, verbose){

  ###
  # CALL: Compute_tTE()
  # DESCRIPTION: This function takes two count matrices (in this case mRNA and tRNA data) and selects the factors (rows) and conditions (columns) that are present in both data sources.
  # This step is crucial during Compute_tTE() as to compute the tTE (correlation) scores matching entries are required.
  ###

  tTE_conditions <- intersect(colnames(data_mRNA), colnames(data_tRNA)) # Extracting the matching conditions

  # Extracting the matching factors
  if (level == "codon"){ # Perform the reverse complement of the codons to find the proper codon-anticodon pairs
    translated_codons <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(rownames(data_mRNA))))
    tTE_factors <- intersect(translated_codons, rownames(data_tRNA))
  } else { # Dealing with AA, perform a direct intersection
    tTE_factors <- intersect(rownames(data_mRNA), rownames(data_tRNA))
  }

  # NO INTERSECTION
  if (length(tTE_conditions) == 0) stop("No intersecting conditions were found between the matrices.\n",
                                        "If needed check that the annotation (columns) between the matrices follows the same format.\n",
                                        "Check the correct use of the `name_sep` input parameter.")
  if (verbose) message(paste("- A total of", length(tTE_conditions), "intersecting conditions were found."))
  if (length(tTE_factors) == 0) stop("No intersecting features were found between the matrices.\n",
                                     "If needed check that the anotation (rows) between the matrices follows the same format.")

  if (verbose) message("- Filtering the datasets by the intersecting conditions.") # In order to filter the mRNA data we need to undo the reverse-complement translation
  filtered_mRNA_data <-  if (level == "codon") data_mRNA[as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(tTE_factors))), tTE_conditions] else data_mRNA[tTE_factors, tTE_conditions]
  filtered_tRNA_data <- data_tRNA[tTE_factors, tTE_conditions]
  CheckDataFrame(data = filtered_mRNA_data)
  CheckDataFrame(data = filtered_tRNA_data)

  return(list(filtered_mRNA_data, filtered_tRNA_data)) # The output consists of a list with the two filtered count matrices
}
