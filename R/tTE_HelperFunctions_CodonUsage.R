#' @importFrom magrittr %>%
#' @importFrom rlang .data
NULL

#' Compute the Exonic Background of the Codon Usage
#' @description
#' This function calculates the **codon usage background** based solely on exonic sequence composition, independent of gene expression levels.
#' It provides a reference distribution of codon frequencies across conditions, used to normalize/compare against observed codon usage patterns derived from expression data.
#'
#' @param data A codon usage matrix with codons as rows and conditions or samples as columns.
#'
#' @return A \code{matrix} with the codon exonic background.
#' @export
#'
#' @examples
#' data(subset_mRNA_data)
#' tTEscanR_obj <- Create_tTEscanR_Object(counts = subset_mRNA_data, assay = "mRNA")
#' tTEscanR_obj <- ComputeCodonUsage(object = tTEscanR_obj, species = "hg38",
#'                                   additional.metrics = FALSE, reduce = 1000)
#' exonic_background <- ComputeCodonExonicBackground(data = tTEscanR_obj@assays$CodonUsage)

ComputeCodonExonicBackground <- function(data){

  ###
  # CALL: User and ComputeAdditionalMetrics_CodonUsage()
  # DESCRIPTION: This function calculates the codon exonic background of a count matrix.
  ###

  if (sum(rowSums(data)) == 0) stop ("- All values in the matrix are 0.") # Check in order to avoid denominator = 0
  exonic_background <- rowSums(data) / sum(rowSums(data)) # Relative contribution of each row to the overall total sum of all values
  return (exonic_background) # Returns a numeric vector with an element per row
}

#' Compute the Correlation Between Mean Codon Usage and Codon Exonic Background
#' @description
#' This function calculates the **correlation** between observed **mean codon usage** and the **codon exonic background**.
#' It provides a metric for evaluating how much codon usage is driven by underlying sequence composition versus condition-specific expression.
#'
#' @param mean A \code{matrix} representing mean codon usage across condition. Can be computed using \code{\link{ComputeMeanUsage}}.
#' @param background A \code{matrix} or table representing codon frequencies in exonic regions. Can be computed using \code{\link{ComputeCodonExonicBackground}}.
#' @param corr_method A correlation method accepted by \code{\link{cor}}. Defaults to \code{"spearman"}.
#'
#' @return Integer; correlation information between \code{mean} and \code{background}.
#' @export
#'
#' @examples
#' data(subset_mRNA_data, metadata)
#' tTEscanR_obj <- Create_tTEscanR_Object(counts = subset_mRNA_data, assay = "mRNA")
#' tTEscanR_obj <- ComputeCodonUsage(object = tTEscanR_obj, species = "hg38",
#'                                   additional.metrics = FALSE, reduce = 1000)
#' exonic_background <- ComputeCodonExonicBackground(data = tTEscanR_obj@assays$CodonUsage)
#' # Input: expression count matrix, need to provide metadata and corr_factor parameters
#' mean_codon_usage <- ComputeMeanUsage(data = tTEscanR_obj@assays$CodonUsage,
#'                                      mode = "raw", metadata = metadata,
#'                                      corr_factor = "tissue")
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
    message("- The correlation will be computed over shared features.\n", "- Filtering steps applied to `mean` and `background` variables.")
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
#' When the input \code{data} is a \code{tTEscanR_Object} the parameters \code{metadata} and \code{corr_factor} will be extracted from the object, and ignored if specified as input parameters.
#' Therefore, variables \code{assay} and \code{metadata} need to be coherent with the rules described in \code{\link{Create_tTEscanR_Object}}.
#'
#' @param data A \code{tTEscanR_Object} or expression count \code{matrix} (with codons, anticodons or amino acids as features).
#' @param assay Optional; a character string specifying the name of the assay to retrieve from the \code{tTEscanR_Object}.
#' @param mode Either \code{"raw"}, \code{"size-corrected"} or \code{"long-format"} to specify the format the input data belongs to. Defaults to \code{"raw"}.
#' @param metadata Optional; a \code{data.frame} with the meta-information related with the conditions in \code{data}. There has to be one column with the same labels as the column names.
#' @param corr_factor Optional; a factor based on \code{metadata} columns to define the variable to correct for. Required if \code{mode} is \code{"raw"} and \code{data} is not a \code{tTEscanR_Object}.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return An updated \code{tTEscanR_Object} if \code{data} is a \code{tTEscanR_Object}. A \code{data.frame} containing a new layer of information representing the mean codon usage if \code{data} is an expression count matrix.
#' @export
#'
#' @examples
#' data(subset_tRNA_data, metadata)
#' tTEscanR_obj <- Create_tTEscanR_Object(counts = subset_tRNA_data, assay = "tRNA",
#'                                        meta.data = list(metadata, "tissue"),
#'                                        meta.data.ids = list("ConditionsLabels", "CorrectionFactor"))
#' # Input: tTEscanR object containing metadata and corr_factor parameters
#' tTEscanR_obj <- ComputeAnticodonUsage(object = tTEscanR_obj)
#' anticodon_mean_usage <- ComputeMeanUsage(data = tTEscanR_obj, assay = "AnticodonUsage")

ComputeMeanUsage <- function(data, assay = NULL, mode = "raw", metadata = NULL, corr_factor = NULL, verbose = TRUE){

  ###
  # CALL: User and ComputeAdditionalMetrics_CodonUsage()
  # DESCRIPTION: This function takes an object or a count matrix to compute the mean usage of their features (codons, anticodons, AA).
  # The data can de input in different formats that need to be properly specified using the parameter mode.
  # The input data can be size-corrected, normalized and transformed into a long format.
  # The output consists of a tibble with the features and their corresponding mean value.
  # If the input is a tTEscanR object this will be updated.
  ###

  # This table links each `assay` parameter to the different layers of data that will need to be retrieved
  assay_map <- list(CodonUsage = list(var_name = "codon", slot = "CodonUsage"),
                    AnticodonUsage = list(var_name = "anticodon", slot = "AnticodonUsage"),
                    AADemand = list(var_name = "AA", slot = "AADemand"),
                    AASupply = list(var_name = "AA", slot = "AASupply"))

  message("A . Evaluating the input.")
  if (inherits(data, 'tTEscanR_Object')){ # Dealing with tTEscanR object
    if (verbose) message("- The input is a tTEscanR object.\n", paste("- Checking if the", assay, "assay is in the object."))
    if (!assay %in% names(assay_map) || is.null(assay)) stop("Please specify a valid `assay` input parameter: CodonUsage, AnticodonUsage, AADemand, AASupply.")

    # Checking the input data
    IsIn_tTEscanR_Object(object = data, slot = "assays", section = assay, verbose = verbose)
    var_name <- assay_map[[assay]]$var_name
    object <- data # Create copy of the object to be able to update it at the end
    data <- data@assays[[assay_map[[assay]]$slot]]
    CheckDataFrame(data)
    if (verbose) message(paste("The assay", assay, "has been properly loaded."))

    # Checking the metadata
    IsIn_tTEscanR_Object(object = object, slot = "meta.data", section = "ConditionsLabels", verbose = FALSE)
    CheckDataFrame(data = object@meta.data$ConditionsLabels)
    metadata <- object@meta.data$ConditionsLabels
    if (verbose) message("- The ConditionsLabels metadata has been properly loaded.")

    # Checking the correction factor
    IsIn_tTEscanR_Object(object = object, slot = "meta.data", section = "CorrectionFactor", verbose = FALSE)
    corr_factor <- object@meta.data$CorrectionFactor

  } else { # Dealing with a dataset
    if (verbose) message("- The input is a count matrix.")

    # Checking the input data and metadata
    CheckDataFrame(data = data)
    if (verbose) message(paste("The assay", assay, "has been properly loaded."))
    CheckDataFrame(data = metadata)
    if (verbose) message("The metadata has been properly loaded.")

    var_name <- "feature"
  }

  if (!(corr_factor %in% colnames(metadata))) stop("The correction factor was not found in the metadata")
  if (verbose) message("- The correction factor has been properly loaded.")
  message("  A . COMPLETED\n", "B . Transforming (if applicable) the input data based on the `mode` parameter.")

  # To compute the mean usage the data needs to be size-corrected (using DESeq2) and normalized (change format).
  if (!(mode %in% c("raw", "size-corrected", "long-format"))) stop("Please specify a suitable `mode` parameter.\n", "Accepted formats: raw, size-corrected, long-format")

  # THE DATA HAS NOT BEEN SIZE-CORRECTED NOT NORMALIZED
  if (mode == "raw"){

    if (verbose) message("- Size correcting the usage matrix to account for sequencing depth.")
    filtered <- FilterByMetadata(data = data, metadata = metadata) # Filter the data and the metadata based on their matching entries
    data <- suppressMessages(ComputeSizeCorrection(data = filtered[[1]], metadata = filtered[[2]], corr_factor = corr_factor, verbose = FALSE))
    mode <- "size-corrected"
  }

  # THE DATA HAS NOT BEEN NORMALIZED
  if (mode == "size-corrected"){
    # Using the same function we normalize the data, and transform it from a wide format into a long format (each row is a different entry)
    if (verbose) message("- Normalizing and transforming to long format the usage matrix.")
    data <- DataToLongFormat(data = data, normalize = TRUE, rownames_to_column = var_name, names_to = "conditions", values_to = "usage")
    mode <- "long-format"
  }

  # THE DATA HAS BEEN SIZE-CORRECTED AND NORMALIZED
  if (mode == "long-format"){
    # Check that the required columns exist in data - Very sensitive to input `data` column's names
    if(!var_name %in% names(data) || !utils::hasName(data, "usage")){
      stop("Missing column or wrong column labelling in the `data` input parameter.\n",
           "Input a suitable long format `data` or specify the correct `mode`.\n",
           paste("Expected (long format) column names: ", var_name, ", conditions, usage."))
    }
    if (verbose) message("- The data is in a proper long format.")
  }

  # The data is size-corrected, normalized and in a long format
  message("  B . COMPLETED\n","C . Calculating the mean usage across conditions.")
  usage_across_conditions <- dplyr::group_by(data, .data[[var_name]]) %>% # Group the data by features (codons, anticodons or AAs)
    dplyr::summarize(mean_usage_across_conditions = mean(.data[["usage"]])) # Compute the mean

  CheckDataFrame(usage_across_conditions, names = FALSE)
  message("  C . COMPLETED")

  if (inherits(data, 'tTEscanR_Object')) { # Update the input tTEscanR object
    object <- suppressMessages(Update_tTEscanR_Object(object = object, meta.data = usage_across_conditions, meta.data.ids = paste(assay, "MeanUsage", sep = "_"), verbose = FALSE))
    return(object)
  } else {
    return(usage_across_conditions) # Returns a tibble with two columns: (i) feature, and (ii) mean_usage_across_conditions
  }
}

ComputeAdditionalMetrics_CodonUsage <- function(codon_usage, codon_freq, metadata, corr_factor = NULL, corr_method, verbose){

  ###
  # CALL: ComputeCodonUsage()
  # DESCRIPTION: This function takes a codon per condition count matrix and a suitable codon frequency per gene matrix as inputs.
  # It computes 3 additional metrics: (i) codon exonic background, (ii) mean codon usage, and (iii) correlation between (i) and (ii)
  ###

  # CODON EXONIC BACKGROUND
  if (verbose) message("- Computing the codon exonic background.")
  codon_exonic_background <- ComputeCodonExonicBackground(data = codon_freq)
  if (is.null(codon_exonic_background)) stop("The codon exonic background could not be computed.")

  # MEAN CODON USAGE
  if (verbose) message("- Computing the mean codon usage.")
  mean_codon_usage <- ComputeMeanUsage(data = codon_usage, metadata = metadata, corr_factor = corr_factor, verbose = FALSE) # verbose = FALSE - Try to avoid too messy outputs
  if (is.null(mean_codon_usage)) stop("The mean codon usage could not be computed.")

  # CORRELATION BACKGROUND-MEAN
  if (verbose) message("- Computing the correlation between the mean codon usage and the codon exonic background.")
  correlation_background_mean <- ComputeCorrelationBackground(mean = mean_codon_usage, background = codon_exonic_background, corr_method = corr_method)
  if (is.null(correlation_background_mean)) stop("The correlation mean-background could not be computed.")

  # Store all the metrics in a named list
  additional_metrics <- list(CodonExonicBackground = codon_exonic_background, MeanCodonUsage = mean_codon_usage, MeanCodonCorr = correlation_background_mean)
  return(additional_metrics)
}

ConsistencyWithCodonFreq <- function(data, codon_freq, filter, species, verbose){

  ###
  # CALL: ComputeCodonUsage()
  # DESCRIPTION: This function checks the format of a mRNA gene expression count matrix and a codon per gene matrix.
  # The codon per gene matrix can be defined by the user or loaded as default from the tTEscanR memory.
  ###

  message("\n------------------------------\n", "A . Assessing the `codon_freq` matrix.")
  codon_frequency_per_gene_table <- CheckCodonFreqTable(data = codon_freq, species = species, filter = filter) # Check the user-defined codon_freq or loads a default table if possible

  message("  A . COMPLETED\n", "B . Evaluating the consistency across gene annotations.")
  if (verbose) message("- Vector 1: codon frequency per gene table\n", "- Vector 2: mRNA gene expression data")

  # The user would needs to input an already translated data if there are inconsistencies in the gene annotation
  gene_annotation <- CheckGeneAnnotation(vector1 = colnames(codon_frequency_per_gene_table), vector2 = rownames(data), verbose = verbose)
  message("  B . COMPLETED\n", "------------------------------\n")

  return(codon_frequency_per_gene_table) # Return the codon (rows) frequency per gene (columns) table
}

CheckNames_tRNA <- function(data){

  ###
  # CALL: ComputeAnticodonUsage()
  # DESCRIPTION: This function checks the format of the tRNA genes annotated.
  # Evaluates their position in the input data as well as that they follow the standard annotation format.
  # There is no output, unless an error is reported.
  ###

  # The tRNA genes need to be located in the rows
  if (is.null(rownames(data))) stop("No tRNA genes found in row names.")

  # Check the actual content of the tRNA gene label
  tRNA_label <- grepl("^tRNA-[A-Za-z]{3,4}-[ATGC]{3}-[0-9]+-[0-9]+$", rownames(data))
  invalid_rows <- which(!tRNA_label)

  # Report if there is any row name that does not follow the requirements above
  if (length(invalid_rows) > 0) stop("Inconsistent tRNA gene format.\n", "Expected format: tRNA-Asn-GTT-5-1\n",
                                     "1st invalid format: ", rownames(data)[invalid_rows[1]], "\n",
                                     "Invalid rows: ", paste(invalid_rows, collapse = ", "))
}
