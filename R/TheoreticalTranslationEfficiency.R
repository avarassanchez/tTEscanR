#' Compute the Theoretical Translation Efficiency (tTE) score
#' @description
#' This function calculates the theoretical translation efficiency (tTE) score by integrating codon-anticodon usage and/or amino acid demand-supply across conditions.
#'
#' @param object A \code{tTEscanR_Object} containing a mRNA and a codon usage assay and/or a tRNA and an anticodon usage assay.
#' @param level Either \code{"codon"}, \code{"aa"} or \code{"both"} to indicate which analysis to perform.
#' @param genetic_code A \code{character} string to specify the genetic code to be used. Defaults to \code{"Standard"}.
#' @param corr_method A correlation method accepted by \code{\link{cor}}. Defaults to \code{"spearman"}.
#' @param compute_significance Logical; if \code{TRUE}, computes the statistical significance (p-value) of the tTE scores. Defaults to \code{TRUE}.
#' @param overwrite Logical; if \code{TRUE}, overwrites any existing tTE table in the \code{tTEscanR_Object}. Defaults to \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return An updated \code{tTEscanR_Object} containing a new layer of information representing the translation efficiency table for the matching conditions in the mRNA and tRNA data.
#' @export
#'
#' @examples
#' data(default_tTEscanR_mRNA_data, default_tTEscanR_tRNA_data, default_tTEscanR_metadata)
#' tTEscanR_obj <- Create_tTEscanR_Object(counts = list(mRNA = default_tTEscanR_mRNA_data,
#'                                                      tRNA = default_tTEscanR_tRNA_data),
#'                                        meta.data = list(default_tTEscanR_metadata, "tissue"),
#'                                        meta.data.ids = list("ConditionsLabels", "CorrectionFactor"))
#' tTEscanR_obj <- ComputeCodonUsage(object = tTEscanR_obj, species = "hg38",
#'                                   additional_metrics = FALSE, reduce = 10000)
#' tTEscanR_obj <- ComputeAnticodonUsage(object = tTEscanR_obj)
#' tTEscanR_obj <- Compute_tTE(object = tTEscanR_obj, level = "codon", compute_significance = FALSE)

Compute_tTE <- function(object, level = c("codon", "aa", "both"), genetic_code = "Standard", corr_method = c("spearman", "pearson", "kendall"), compute_significance = TRUE, overwrite = FALSE, verbose = TRUE){

  ###
  # CALL: User or Run_tTEscanR_pipeline()
  # DESCRIPTION: This function computes the theoretical translation efficiency (tTE) score by correlating codon-anticodon usage and/or amino acid demand-supply across conditions.
  ###

  message("\n--- Computation of the theoretical translation efficiency (tTE) ---", "\n1 . Checking the format of the input data.")
  if (!(inherits(object, 'tTEscanR_Object'))) stop("'object' must be a tTEscanR object.")
  if (verbose) message("- The input contains a proper tTEscanR object.")

  corr_method <- match.arg(corr_method)
  level <- match.arg(level)
  # if (!(corr_method %in% c("spearman", "pearson", "kendall"))) stop("Please specify a suitable `corr_method`.\n", "For further details check the cor() documentation.")

  target_levels <- if (level == "both") names(assay_map_TE) else level
  if (is.null(target_levels) || !all(target_levels %in% names(assay_map_TE))) stop("Please specify a valid 'level' input parameter: codon, aa, both.")

  # Evaluate that the required data assays are present in the object
  IsIn_tTEscanR_Object(object = object, slot = "assays", section = "tRNA", update_assay = FALSE, overwrite = FALSE, verbose = FALSE)

  # Evaluate that the required metadata variables are present in the object
  IsIn_tTEscanR_Object(object = object, slot = "meta.data", section = "ConditionsLabels", verbose = FALSE)
  IsIn_tTEscanR_Object(object = object, slot = "meta.data", section = "CorrectionFactor")
  meta <- object@meta.data$ConditionsLabels
  batch <- object@meta.data$CorrectionFactor
  if (!(batch %in% colnames(meta))) stop("The correction factor was not found in the metadata.")

  final_results <- list() # Empty list to store the results
  message("1 . COMPLETED\n", "2 . Computing the theoretical translation efficiency (tTE) analysis.")
  for (i in target_levels){
    if (verbose) message("\nProcessing level: ", i)

    config <- assay_map_TE[[i]]
    mRNA_raw <- object@assays[[config$mRNA]]
    tRNA_raw <- object@assays[[config$tRNA]]

    CheckDataFrame(data = mRNA_raw)
    CheckDataFrame(data = tRNA_raw)

    filter_results <- FilterMatrix(data_mRNA = mRNA_raw, data_tRNA = tRNA_raw, level = i, verbose = verbose) # Retain matching conditions in the mRNA and tRNA data

    if (verbose) message("- Filtering the metadata by the intersecting conditions.") # Filter the data and metadata based on the shared content
    filtered_data <- FilterByMetadata(data = filter_results$mRNA, metadata = meta, verbose = FALSE) # Evaluating the mRNA data

    if (verbose) message("- Size-correcting the data.") # The metadata is just retrieved once as it has been previously checked the consistency between mRNA and tRNA datasets
    mRNA_filtered <- ComputeSizeCorrection(data = filtered_data$data, metadata = filtered_data$metadata, batch = batch, verbose = FALSE)
    tRNA_filtered <- ComputeSizeCorrection(data = filter_results$tRNA, metadata = filtered_data$metadata, batch = batch, verbose = FALSE)

    # Compute the correlation between the mRNA and the tRNA data
    tTE_results_tibble <- ComputeCorrelation(data_mRNA = mRNA_filtered, data_tRNA = tRNA_filtered, corr_method = corr_method, verbose = verbose)

    if (isTRUE(compute_significance)){
      if (verbose) message("- Computing the statistical significance (this may take time)...")

      # Filter the tRNA expression data based on the shared conditions extracted above
      tRNA_exp_data <- object@assays$tRNA[, colnames(filter_results$mRNA), drop = FALSE] # Usage of tRNA expression matrix to increment the number of features considered in the shuffling analysis
      CheckDataFrame(data = tRNA_exp_data)
      tRNA_exp_data <- suppressMessages(ComputeSizeCorrection(data = tRNA_exp_data, metadata = filtered_data$metadata, batch = batch, verbose = FALSE))

      tTE_results_tibble <- ComputeStatisticalSignificance(level = i, tTE_scores = tTE_results_tibble, data_mRNA = mRNA_filtered, genetic_code = genetic_code,
                                                           data_tRNA = tRNA_filtered, tRNA_exp = tRNA_exp_data, corr_method = corr_method, verbose = verbose)
    }

    final_results[[config$id]] <- tTE_results_tibble
    if (verbose) message("\nCompleted the processing level: ", i)
  }

  message("\n2 . COMPLETED\n", "--- The theoretical translation efficiency (tTE) was properly computed ---\n")
  if (verbose) message("- Updating the tTEscanR object.")
  object <- Update_tTEscanR_Object(object = object, meta.data = final_results, meta.data.ids = names(final_results), overwrite = overwrite, verbose = FALSE)
  return(object) # The output tTEscanR object has been validated in Update_tTEscanR_Object()
}

ComputeStatisticalSignificance <- function(level, tTE_scores, data_mRNA, data_tRNA, tRNA_exp, corr_method, genetic_code, N_perm = 1000, verbose){

  ###
  # CALL: Compute_tTE()
  # DESCRIPTION: This function takes computes a statistical test to obtain the significance of tTE scores  by comparing the mRNA data (codon usage or AA demand) to the tRNA expression data.
  # The tRNA expression data is necessary to increase the number of features to consider to perform the randomization.
  # If performed with the anticodon usage or AA supply less features could be interrogated, and therefore we will achieve less power.
  ###

  tTE_conditions <- colnames(data_tRNA)
  n_cond <- length(tTE_conditions)
  tRNA_names <- rownames(tRNA_exp)
  tRNAsAnticodons <- sapply(strsplit(tRNA_names, "-"), "[[", 3) # Extract the anticodons from the tRNA data

  if (level == "aa") {
    if (!(genetic_code %in% colnames(final_matrix_genetic_code))) stop("The genetic code '", genetic_code, "' is not a valid idnetifier.")
    aa_dictionary <- stats::setNames(final_matrix_genetic_code[[genetic_code]], final_matrix_genetic_code$Anticodon)
    tRNA_map_level <- unname(aa_dictionary[tRNAsAnticodons])
    if (any(is.na(tRNA_map_level))) warning("Some anticodons mapped to codons missing from the genetic code table. They will be ignored.")

  } else {
    codon_to_anti_dict <- stats::setNames(final_matrix_genetic_code$Anticodon, final_matrix_genetic_code$Codon)
    tRNA_map_level <- tRNAsAnticodons
    rownames(data_mRNA) <- unname(codon_to_anti_dict[rownames(data_mRNA)])
  }

  common_features <- intersect(rownames(data_mRNA), unique(tRNA_map_level))
  data_mRNA <- as.matrix(data_mRNA[common_features, , drop = FALSE])

  keep_tRNA_indices <- which(tRNA_map_level %in% common_features)
  tRNA_exp_filtered <- as.matrix(tRNA_exp[keep_tRNA_indices, , drop = FALSE])
  tRNA_map_level_filtered <- tRNA_map_level[keep_tRNA_indices]

  # Initialize the output matrix
  map_factor <- factor(tRNA_map_level_filtered , levels = common_features)
  M <- Matrix::sparse.model.matrix(~ 0 + map_factor)
  M <- Matrix::t(M) # Features as rows and columns as tRNA genes

  tTE_p_values <- numeric(n_cond) # Generate another empty vector (all values = 0), this time to store the significance results

  pb <- utils::txtProgressBar(min = 0, max = n_cond, style = 3) # Start a time counter to track the progress of the function execution
  for (i in seq_len(n_cond)) { # Iterate over the conditions present in the data

    curr_cond <- tTE_conditions[i] # Retrieve the condition interrogated each time
    mRNA_vec <- data_mRNA[, curr_cond]  # Observed (true) usage for this condition
    tRNA_vec_raw <- tRNA_exp_filtered[, curr_cond] # Observed (true) tRNA gene usage for this condition: this will be shuffled N_perm times to determine statistical significance
    n_tRNA_genes <- length(tRNA_vec_raw)
    null_correlations <- numeric(N_perm) # Generate another empty vector (all values = 0), this time to store the values of each permutation
    # observed_score <- as.numeric(tTE_scores[i])

    for (j in seq_len(N_perm)) { # Iterates over the number of iterations defined by parameter N_perm
      shuffled_tRNA <- tRNA_vec_raw[sample.int(n_tRNA_genes)] # Randomizes the order of the tRNA genes
      agg_usage <- as.numeric(M %*% shuffled_tRNA) # Aggregated usage - sum of tRNA genes per feature
      # agg_usage_rel <- agg_usage / sum(agg_usage) # Normalization of the relative usage

      # Fills the value of the matrix with the tTE score (correlation)
      null_correlations[j] <- stats::cor(agg_usage, mRNA_vec, method = corr_method) # Correlation between shuffled tRNA data and unshuffled mRNA data (AA demand or codon usage)
    }

    # Fit parameters (mu and sigma for null distribution)
    mean_random_correlations <- mean(null_correlations)
    sd_random_correlations <- stats::sd(null_correlations)

    # Calculate p-value for the actual tTE
    tTE_p_values[i] <- stats::pnorm(tTE_scores[i], mean = mean_random_correlations, sd = sd_random_correlations, lower.tail = FALSE)

    utils::setTxtProgressBar(pb, i) # Increase the progress bar
  }

  close(pb)

  # Generate summary table - table with the input tTE scores and their corresponding p-values
  return(tibble::tibble(condition = tTE_conditions, tTE = as.numeric(tTE_scores), p_value = tTE_p_values, neg_log10_tTE_p_value = -log10(tTE_p_values)))
}

ComputeCorrelation <- function(data_mRNA, data_tRNA, corr_method, verbose){

  ###
  # CALL: Compute_tTE()
  # DESCRIPTION: This function takes two count matrices and performs the correlation between their entries.
  ###

  if (verbose) message(paste("- Computing the", corr_method, "correlation."))
  tTE_conditions <- colnames(data_mRNA) # Extract the features and the conditions to consider - both datasets have been filtered for common conditions

  # Selects the condition that we are analyzing in each iteration and filters the data matrices accordingly
  tTEs_at_level <- vapply(seq_along(tTE_conditions), function(i) {
    stats::cor(data_mRNA[, i], data_tRNA[, i], method = corr_method)
  }, FUN.VALUE = numeric(1))
  names(tTEs_at_level) <- tTE_conditions # Include as names of the vector the elements in tTE_conditions

  if (any(is.na(tTEs_at_level))) warning("Some conditions produced NA correlation scores. Check for zero variance in those columns.")
  if (is.null(tTEs_at_level)) stop("The tTE scores could not be computed.")
  return(tTEs_at_level) # The output consists of vector with the features and the correlation value, the theoretical translation efficiency (tTE) score
}

FilterMatrix <- function(data_mRNA, data_tRNA, level = c("codon", "aa"), verbose){

  ###
  # CALL: Compute_tTE()
  # DESCRIPTION: This function takes two count matrices (in this case mRNA and tRNA data) and selects the factors (rows) and conditions (columns) that are present in both data sources.
  # This step is crucial during Compute_tTE() as to compute the tTE (correlation) scores matching entries are required.
  ###

  level <- match.arg(level)
  tTE_conditions <- intersect(colnames(data_mRNA), colnames(data_tRNA)) # Extracting the matching conditions
  if (length(tTE_conditions) == 0) stop("No intersecting conditions were found between the matrices.\n",
                                        "If needed check that the annotation (columns) between the matrices follows the same format.")

  # Extracting the matching factors
  row_mRNA <- rownames(data_mRNA)
  row_tRNA <- rownames(data_tRNA)

  if (level == "codon"){ # Perform the reverse complement of the codons to find the proper codon-anticodon pairs
    codon_to_anti_dict <- stats::setNames(final_matrix_genetic_code$Anticodon, final_matrix_genetic_code$Codon)
    translated_codons <- unname(codon_to_anti_dict[row_mRNA])
    # translated_codons <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(row_mRNA)))
    tTE_factors <- intersect(translated_codons, row_tRNA)
    if (length(tTE_factors) == 0) stop("No intersecting features found for codons.")

    valid_index <- translated_codons %in% tTE_factors
    mRNA_filtered_rows <- row_mRNA[valid_index]
    tRNA_filtered_rows <- translated_codons[valid_index]

  } else { # Dealing with AA, perform a direct intersection
    tTE_factors <- intersect(row_mRNA, row_tRNA)
    if (length(tTE_factors) == 0) stop("No intersecting features found for amino acids")
    mRNA_filtered_rows <- tTE_factors
    tRNA_filtered_rows <- tTE_factors
  }

  if (verbose) message(paste("- Found", length(tTE_conditions), "conditions and", length(tRNA_filtered_rows), "features."))
  if (verbose) message("- Filtering the datasets by the intersecting conditions.")
  filtered_mRNA <- data_mRNA[mRNA_filtered_rows, tTE_conditions, drop = FALSE]
  filtered_tRNA <- data_tRNA[tRNA_filtered_rows, tTE_conditions, drop = FALSE]

  CheckDataFrame(data = filtered_mRNA)
  CheckDataFrame(data = filtered_tRNA)

  # filtered_row_mRNA <- rownames(filtered_mRNA)
  # if (level == "codon") {
  #   reordered_tRNA_names <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(filtered_row_mRNA)))
  #   filtered_tRNA <- filtered_tRNA[reordered_tRNA_names, , drop = FALSE]
  # } else {
  #   filtered_tRNA <- filtered_tRNA[filtered_row_mRNA, , drop = FALSE]
  # }

  return(list(mRNA = filtered_mRNA, tRNA = filtered_tRNA)) # The output consists of a list with the two filtered count matrices
}
