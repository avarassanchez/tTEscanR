#' @importFrom magrittr %>%
#' @importFrom rlang .data :=
NULL

#' Compute the Exonic Background of the Codon Usage
#' @description
#' This function calculates the **codon/anticodon usage background** based solely on exonic sequence composition, independent of expression levels.
#' It provides a reference distribution of codon/anticodon frequencies across conditions, used to normalize/compare against observed usage patterns derived from expression data.
#'
#' @param data A codon usage matrix with codons as rows and conditions or samples as columns.
#'
#' @return A \code{matrix} with the codon exonic background.
#' @export
#'
#' @examples
#' data(default_tTEscanR_mRNA_data)
#' tTEscanR_obj <- Create_tTEscanR_Object(counts = default_tTEscanR_mRNA_data, assay = "mRNA")
#' tTEscanR_obj <- ComputeCodonUsage(object = tTEscanR_obj, species = "hg38",
#'                                   additional_metrics = FALSE, reduce = 1000)
#' exonic_background <- ComputeExonicBackground(data = tTEscanR_obj@assays$CodonUsage)

ComputeExonicBackground <- function(data){

  ###
  # CALL: User and ComputeAdditionalMetrics_CodonUsage()
  # DESCRIPTION: This function calculates the codon exonic background of a count matrix.
  ###

  row_totals <- rowSums(data)
  total <- sum(row_totals)

  if (total == 0) stop ("- All values in the matrix are 0.") # Check in order to avoid denominator = 0
  exonic_background <- row_totals / total # Relative contribution of each row to the overall total sum of all values
  return (exonic_background) # Returns a numeric vector with an element per row
}

#' Compute the Correlation Between Mean Usage and Exonic Background
#' @description
#' This function calculates the **correlation** between observed **mean usage** and the **exonic background**.
#' It provides a metric for evaluating how much usage is driven by underlying sequence composition versus condition-specific expression.
#'
#' @param mean A \code{matrix} representing mean usage across condition. Can be computed using \code{\link{ComputeMeanUsage}}.
#' @param background A \code{matrix} or table representing frequencies in exonic regions. Can be computed using \code{\link{ComputeExonicBackground}}.
#' @param corr_method A correlation method accepted by \code{\link{cor}}. Defaults to \code{"spearman"}.
#'
#' @return Integer; correlation information between \code{mean} and \code{background}.
#' @export
#'
#' @examples
#' data(default_tTEscanR_mRNA_data, default_tTEscanR_metadata)
#' tTEscanR_obj <- Create_tTEscanR_Object(counts = default_tTEscanR_mRNA_data, assay = "mRNA")
#' tTEscanR_obj <- ComputeCodonUsage(object = tTEscanR_obj, species = "hg38",
#'                                   additional_metrics = FALSE, reduce = 1000)
#' exonic_background <- ComputeExonicBackground(data = tTEscanR_obj@assays$CodonUsage)
#' # Input: expression count matrix, need to provide metadata and batch parameters
#' mean_codon_usage <- ComputeMeanUsage(data = tTEscanR_obj@assays$CodonUsage,
#'                                      mode = "raw", metadata = default_tTEscanR_metadata,
#'                                      batch = "tissue")
#' correlation_background <- ComputeCorrelationBackground(mean = mean_codon_usage,
#'                                                        background = exonic_background)

ComputeCorrelationBackground <- function(mean, background, corr_method = "spearman"){

  ###
  # CALL: User and ComputeAdditionalMetrics_CodonUsage()
  # DESCRIPTION: This function takes a mean and background vector and computes their correlation.
  # The default correlation method can be modified by the user based on the documentation of stats::cor().
  ###

  mean_vector <- stats::setNames(mean[[2]], mean[[1]]) # Convert the `mean` input from a tibble to a named vector
  common_features <- intersect(names(mean_vector), names(background)) # Extract the common features between mean and background

  # If needed filter the vectors accordingly
  if (length(common_features) != length(mean_vector) || length(common_features) != length(background)){
    message("- The correlation will be computed over shared features.\n", "- Filtering steps applied to 'mean' and 'background' variables.")
    mean_vector <- mean_vector[common_features]
    background <- background[common_features]
  }

  correlation_to_exonic_background <- round(stats::cor(mean_vector, background, method = corr_method), 3) # Retain 3 decimals
  return(correlation_to_exonic_background)
}

#' Compute Usage Across Conditions (Mean Usage)
#' @description
#' This function computes the **average usage** of codon, anticodons, or amino acids across conditions, useful for summarizing feature usage trends across sample groups.
#' It supports direct input of a count matrix and extraction from a \code{tTEscanR_Object}.
#' When the input \code{data} is a \code{tTEscanR_Object} the parameters \code{metadata} and \code{batch} will be extracted from the object, and ignored if specified as input parameters.
#' Therefore, variables \code{assay} and \code{metadata} need to be coherent with the rules described in \code{\link{Create_tTEscanR_Object}}.
#'
#' @param data A \code{tTEscanR_Object} or expression count \code{matrix} (with codons, anticodons or amino acids as features).
#' @param assay Optional; a character string specifying the name of the assay to retrieve from the \code{tTEscanR_Object}.
#' @param mode Either \code{"raw"}, \code{"size-corrected"} or \code{"long-format"} to specify the format the input data belongs to. Defaults to \code{"raw"}.
#' @param metadata Optional; a \code{data.frame} with the meta-information related with the conditions in \code{data}. There has to be one column with the same labels as the column names.
#' @param id_col Optional; a factor based on \code{metadata} columns to define the variable to use to link it with the \code{data}. If \code{NULL} the column with the highest agreement will be automatically selected.
#' @param batch Optional; a factor based on \code{metadata} columns to define the variable to correct for. Required if \code{mode} is \code{"raw"} and \code{data} is not a \code{tTEscanR_Object}.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return An updated \code{tTEscanR_Object} if \code{data} is a \code{tTEscanR_Object}. A \code{data.frame} containing a new layer of information representing the mean codon usage if \code{data} is an expression count matrix.
#' @export
#'
#' @examples
#' data(default_tTEscanR_tRNA_data, default_tTEscanR_metadata)
#' tTEscanR_obj <- Create_tTEscanR_Object(counts = default_tTEscanR_tRNA_data, assay = "tRNA",
#'                                        meta.data = list(default_tTEscanR_metadata, "tissue"),
#'                                        meta.data.ids = list("ConditionsLabels", "CorrectionFactor"))
#' # Input: tTEscanR object containing metadata and batch parameters
#' tTEscanR_obj <- ComputeAnticodonUsage(object = tTEscanR_obj)
#' anticodon_mean_usage <- ComputeMeanUsage(data = tTEscanR_obj, assay = "AnticodonUsage")

ComputeMeanUsage <- function(data, assay = NULL, mode = c("raw", "size-corrected", "long_format"), metadata = NULL, id_col = NULL, batch = NULL, verbose = TRUE){

  ###
  # CALL: User and ComputeAdditionalMetrics_CodonUsage()
  # DESCRIPTION: This function takes an object or a count matrix to compute the mean usage of their features (codons, anticodons, AA).
  # The data can de input in different formats that need to be properly specified using the parameter mode.
  # The input data can be size-corrected, normalized and transformed into a long format.
  # The output consists of a tibble with the features and their corresponding mean value.
  # If the input is a tTEscanR object this will be updated.
  ###

  # if (!mode %in% c("raw", "size-corrected", "long-format")) stop("Please specify a suitable `mode` parameter.\n", "Accepted formats: raw, size-corrected, long-format")
  mode <- match.arg(mode)

  is_object <- inherits(data, 'tTEscanR_Object')
  assay_map <- list(CodonUsage = "codon", AnticodonUsage = "anticodon", AADemand = "AA", AASupply = "AA") # Links each 'assay' to the data to be retrieved

  if (verbose) message("A . Evaluating the input.")
  if (is_object) { # Dealing with a tTEscanR object
    if (verbose) message("- The input is a tTEscanR object.")
    if (!assay %in% names(assay_map)) stop("Invalid 'assay'.\n", "Please specify a valid 'assay' input parameter: CodonUsage, AnticodonUsage, AADemand, AASupply.")

    raw_mat <- data@assays[[assay]]
    metadata <- data@meta.data$ConditionsLabels
    batch <- data@meta.data$CorrectionFactor
    var_name <- assay_map[[assay]]

    IsIn_tTEscanR_Object(object = data, slot = "assays", section = assay, verbose = verbose)
    IsIn_tTEscanR_Object(object = data, slot = "meta.data", section = "CorrectionFactor", verbose = FALSE)
    IsIn_tTEscanR_Object(object = data, slot = "meta.data", section = "ConditionsLabels", verbose = FALSE)

  } else { # Dealing with a dataset
    raw_mat <- data
    var_name <- "feature"
  }

  if (verbose) message("A . COMPLETED\n", "B . Calculating the mean usage across conditions.")

  if (mode == "raw"){ # THE DATA HAS NOT BEEN SIZE-CORRECTED NOT NORMALIZED
    if (verbose) message("- Size correcting the matrix to account for sequencing depth.")
    if (!is.null(id_col)) {
      if (!id_col %in% colnames(metadata)) is_col <- NULL
      message("The 'id_col' was not found in the 'metadata'. Automatic detection will be implemented.")
    }
    filtered <- FilterByMetadata(data = raw_mat, metadata = metadata, id_col = id_col) # Filter the data and the metadata based on their matching entries
    raw_mat <- suppressMessages(ComputeSizeCorrection(data = filtered[[1]], metadata = filtered[[2]], batch = batch, verbose = FALSE))
    mode <- "size-corrected"
  }

  if (mode == "size-corrected"){ # THE DATA HAS NOT BEEN NORMALIZED
    if (verbose) message("- Normalizing the usage matrix.")
    # raw_mat <- TransformFormat(data = raw_mat, normalize = TRUE, rownames_to_column = var_name, names_to = "conditions", values_to = "usage")

    norm_mat <- t(t(raw_mat)  / colSums(raw_mat)) # Normalization
    means <- rowMeans(norm_mat, na.rm = TRUE) # Mean calculation
    usage_results <- tibble::tibble(!!var_name := names(means), mean_usage_across_conditions = as.numeric(means))

  } else if (mode == "long-format"){ # THE DATA HAS BEEN SIZE-CORRECTED AND NORMALIZED
    # if (!var_name %in% names(raw_mat) || !utils::hasName(raw_mat, "usage")){ # Check that the required columns exist in data - Very sensitive to input `data` column's names
    #   stop("Missing column or wrong column labeling in the `data` input parameter.\n",
    #        "Input a suitable long format `data` or specify the correct `mode`.\n",
    #        paste("Expected (long format) column names: ", var_name, ", conditions, usage."))
    # }
    usage_results <- data %>% dplyr::group_by(.data[[var_name]]) %>%
      dplyr::summarize(mean_usage_across_conditions = mean(.data[["usage"]], na.rm = TRUE)) # Group the data by features (codons, anticodons or AAs) and compute the mean
    if (verbose) message("- The data is in a proper long format.")
  }

  CheckDataFrame(usage_results, required_names = FALSE)
  if (verbose) message("B . COMPLETED")
  if (is_object) { # Update the input tTEscanR object
    data <- suppressMessages(Update_tTEscanR_Object(object = data, meta.data = usage_results, meta.data.ids = paste(assay, "MeanUsage", sep = "_"), verbose = FALSE))
    return(data)
  } else {
    return(usage_results) # Returns a tibble with two columns: (i) feature, and (ii) mean_usage_across_conditions
  }
}

ComputeAdditionalMetrics_CodonUsage <- function(codon_usage, codon_freq, metadata, id_col = NULL, batch = NULL, corr_method, verbose){

  ###
  # CALL: ComputeCodonUsage()
  # DESCRIPTION: This function takes a codon per condition count matrix and a suitable codon frequency per gene matrix as inputs.
  # It computes 3 additional metrics: (i) codon exonic background, (ii) mean codon usage, and (iii) correlation between (i) and (ii)
  ###

  # CODON EXONIC BACKGROUND
  if (verbose) message("- Computing the codon exonic background.")
  codon_exonic_background <- ComputeExonicBackground(data = codon_freq)
  if (is.null(codon_exonic_background)) stop("The codon exonic background could not be computed.\n", "Failure in ComputeExonicBackground().")

  # MEAN CODON USAGE
  if (verbose) message("- Computing the mean codon usage.")
  mean_codon_usage <- ComputeMeanUsage(data = codon_usage, metadata = metadata, id_col = id_col, batch = batch, verbose = FALSE)
  if (is.null(mean_codon_usage)) stop("The mean codon usage could not be computed.\n", "Failure in ComputeMeanUsage().")

  # CORRELATION BACKGROUND-MEAN
  if (verbose) message("- Computing the correlation between the mean codon usage and the codon exonic background.")
  correlation_background_mean <- ComputeCorrelationBackground(mean = mean_codon_usage, background = codon_exonic_background, corr_method = corr_method)
  if (is.null(correlation_background_mean)) stop("The correlation mean-background could not be computed.\n", "Failure in ComputeCorrelationBackground().")

  return(list(CodonExonicBackground = codon_exonic_background, MeanCodonUsage = mean_codon_usage, MeanCodonCorr = correlation_background_mean)) # Store all the metrics in a named list
}

ConsistencyWithCodonFreq <- function(data, codon_freq, species, verbose){

  ###
  # CALL: ComputeCodonUsage()
  # DESCRIPTION: This function checks the format of a mRNA gene expression count matrix and a codon per gene matrix.
  # The codon per gene matrix can be defined by the user or loaded as default from the tTEscanR memory.
  ###

  if (verbose) message("\n------------------------------\n", "A . Assessing the 'codon_freq' matrix.")
  codon_frequency_per_gene_table <- CheckCodonFreqTable(data = codon_freq, species = species) # Checks the user-defined codon_freq or loads a default table if possible
  if (verbose) message("A . COMPLETED\n", "B . Evaluating the consistency across gene annotations.", "- Vector 1: codon frequency per gene table\n", "- Vector 2: mRNA gene expression data")

  gene_annotation <- CheckGeneAnnotation(vector1 = colnames(codon_frequency_per_gene_table), vector2 = rownames(data), verbose = verbose) # The user will need to input an already translated data if there are inconsistencies in the gene annotation
  if (verbose) message("B . COMPLETED\n", "------------------------------\n")

  return(codon_frequency_per_gene_table) # Return the codon (rows) frequency per gene (columns) table
}

CheckNames_tRNA <- function(gene_names){

  ###
  # CALL: ComputeAnticodonUsage()
  # DESCRIPTION: This function checks the format of the tRNA genes annotated (position and annotation format).
  # There is no output, unless an error is reported.
  ###

  if (is.null(gene_names)) stop("No tRNA genes found in row names.")

  pattern <- "^tRNA-[A-Za-z]{3,4}-[ATGC]{3}(-[0-9]+-[0-9]+)?$" # Check the actual content of the tRNA gene label - accepts both tRNA gene or isoacceptor names
  # pattern <- "^([Mm][Tt]-)?(tRNA-|T|trn)?[A-Za-z]{1,3}[0-9]?(-[ATGC]{3})?(-[0-9]+(-[0-9]+)?)?$"
  invalid_rows <- grep(pattern, gene_names, invert = TRUE)

  # Report if there is any row name that does not follow the requirements above
  if (length(invalid_rows) > 0) stop("Inconsistent tRNA gene format.\n", "Expected format: tRNA-Asn-GTT-5-1\n",
                                     "1st invalid format: ", gene_names[invalid_rows[1]], "\n",
                                     "Invalid rows: ", paste(invalid_rows, collapse = ", "))
}
