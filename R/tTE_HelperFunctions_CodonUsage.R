#' @importFrom magrittr %>%
#' @importFrom rlang .data :=
NULL

#' Compute the Exonic Background of the Codon/Anticodon Usage
#' @description
#' This function calculates the \strong{codon/anticodon usage background} based
#' solely on exonic sequence composition, independent of expression levels.
#' It provides a reference distribution of codon/anticodon frequencies across
#' conditions, used to normalize/compare against observed usage patterns derived
#' from expression data.
#'
#' @param data A codon usage matrix with codons as rows and conditions or
#' samples as columns.
#'
#' @return A \code{matrix} with the codon/anticodon background.
#' @export
#'
#' @examples
#' data(default_tTEscanR_mRNA_data)
#' tTEscanR_obj <- CreateObject(
#'     counts = default_tTEscanR_mRNA_data,
#'     assay = "mRNA"
#' )
#' tTEscanR_obj <- ComputeCodonUsage(
#'     object = tTEscanR_obj, species = "hg38",
#'     additional_metrics = FALSE, reduce = 1000
#' )
#' exonic_background <- ComputeExonicBackground(data = getAssay(
#'     tTEscanR_obj,
#'     "CodonUsage"
#' ))
ComputeExonicBackground <- function(data) {
    row_totals <- rowSums(data)
    total <- sum(row_totals)

    # Check in order to avoid denominator = 0
    if (total == 0) stop("- All values in the matrix are 0.")

    # Relative contribution of each row to the overall total sum of all values
    exonic_background <- row_totals / total

    # Returns a numeric vector with an element per row
    return(exonic_background)
}

#' Compute the Correlation Between Mean Usage and Exonic Background
#' @description
#' This function calculates the \strong{correlation} between observed
#' \strong{mean usage} and the \strong{exonic background}. It provides a metric
#' for evaluating how much usage is driven by underlying sequence composition
#' versus condition-specific expression.
#'
#' @param mean A \code{matrix} containing the mean usage across conditions. Can
#' be computed using \code{\link{ComputeMeanUsage}}.
#' @param background A \code{matrix} or table containing the frequencies in
#' exonic regions. Can be computed using \code{\link{ComputeExonicBackground}}.
#' @param corr_method A correlation method accepted by \code{\link{cor}}.
#' Defaults to \code{"spearman"}.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#' Defaults to \code{TRUE}.
#'
#' @return Integer; correlation information between \code{mean} and
#' \code{background}.
#' @export
#'
#' @examples
#' data(default_tTEscanR_mRNA_data, default_tTEscanR_metadata)
#' tTEscanR_obj <- CreateObject(
#'     counts = default_tTEscanR_mRNA_data,
#'     assay = "mRNA"
#' )
#' tTEscanR_obj <- ComputeCodonUsage(
#'     object = tTEscanR_obj, species = "hg38",
#'     additional_metrics = FALSE, reduce = 1000
#' )
#' codon_usage <- getAssay(tTEscanR_obj, "CodonUsage")
#' exonic_background <- ComputeExonicBackground(data = codon_usage)
#' # Input: expression count matrix, need to provide metadata & batch parameters
#' mean_codon_usage <- ComputeMeanUsage(
#'     data = codon_usage,
#'     mode = "raw",
#'     metadata = default_tTEscanR_metadata,
#'     batch = "tissue"
#' )
#' corr_back <- ComputeCorrelationBackground(
#'     mean = mean_codon_usage,
#'     background = exonic_background
#' )
ComputeCorrelationBackground <- function(
    mean, background, corr_method = "spearman", verbose = TRUE
) {
    #  Convert the 'mean' input from a tibble to a named vector
    mean_vector <- stats::setNames(mean[[2]], mean[[1]])

    # Extract the common features between mean and background
    common_features <- intersect(names(mean_vector), names(background))

    # If needed filter the vectors accordingly
    if (length(common_features) != length(mean_vector) ||
        length(common_features) != length(background)) {
        if (verbose) {
            message(
                "- The correlation will be computed over shared features.\n",
                "- Filtering steps applied to 'mean' and 'background'."
            )
        }
        mean_vector <- mean_vector[common_features]
        background <- background[common_features]
    }

    corr_to_exonic_back <- round(stats::cor(mean_vector, background,
        method = corr_method
    ), 3) # 3 decimals
    return(corr_to_exonic_back)
}

#' Compute Usage Across Conditions (Mean Usage)
#' @description
#' This function computes the \strong{average usage} of codons, anticodons, or
#' amino acids across conditions, useful for summarizing feature usage trends
#' across sample groups. It supports direct input of a count matrix and
#' extraction from a \code{tTEscanR_Object}. When the input \code{data} is a
#' \code{tTEscanR_Object} the parameters \code{metadata} and \code{batch} will
#' be extracted from the object, and ignored if specified as input parameters.
#' Therefore, variables \code{assay} and \code{metadata} need to be coherent
#' with the rules described in \code{\link{CreateObject}}.
#'
#' @param data A \code{tTEscanR_Object} or expression count \code{matrix}
#' (with codons, anticodons or amino acids as features).
#' @param assay Optional; a character string specifying the name of the assay
#' to retrieve from the \code{tTEscanR_Object}.
#' @param mode Either \code{"raw"}, \code{"size-corrected"} or
#' \code{"long-format"} to specify the format that the input data belongs to.
#' Defaults to \code{"raw"}.
#' @param metadata Optional; a \code{data.frame} with the meta-information
#' related with the conditions in \code{data}. There has to be one column with
#' the same labels as the column names.
#' @param id_col Optional; a factor based on \code{metadata} columns to define
#' the variable to use to link it with the \code{data}. If \code{NULL} the
#' column with the highest agreement will be automatically selected.
#' @param batch Optional; a factor based on \code{metadata} columns to define
#' the variable to correct for. Required if \code{mode} is \code{"raw"} and
#' \code{data} is not a \code{tTEscanR_Object}.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#' Defaults to \code{TRUE}.
#'
#' @return An updated \code{tTEscanR_Object} if \code{data} is a
#' \code{tTEscanR_Object}. A \code{data.frame} containing a new layer of
#' information representing the mean codon usage if \code{data} is an expression
#' count matrix.
#' @export
#'
#' @examples
#' data(default_tTEscanR_tRNA_data, default_tTEscanR_metadata)
#' tTEscanR_obj <- CreateObject(
#'     counts = default_tTEscanR_tRNA_data,
#'     assay = "tRNA",
#'     meta.data = list(default_tTEscanR_metadata, "tissue"),
#'     meta.data.ids = list(
#'         "ConditionsLabels",
#'         "CorrectionFactor"
#'     )
#' )
#' # Input: tTEscanR object containing metadata and batch parameters
#' tTEscanR_obj <- ComputeAnticodonUsage(object = tTEscanR_obj)
#' anticodon_mean_usage <- ComputeMeanUsage(
#'     data = tTEscanR_obj,
#'     assay = "AnticodonUsage"
#' )
ComputeMeanUsage <- function(
    data, assay = NULL, metadata = NULL, id_col = NULL, batch = NULL,
    mode = c("raw", "size-corrected", "long_format"), verbose = TRUE
) {
    mode <- match.arg(mode)
    is_object <- inherits(data, "tTEscanR_Object")
    extract_data <- GeneralChecksMeanUsage( # Performs step A
        data = data, assay = assay, metadata = metadata, batch = batch,
        verbose = verbose
    )
    mat <- extract_data$raw_mat
    var_name <- extract_data$var_name

    if (verbose) message("B . Calculating the mean usage across conditions.")
    if (mode == "raw") { # Data has not been size-corrected not norm.
        mat <- ModeRaw(
            id_col = id_col, metadata = extract_data$metadata,
            raw_mat = mat, batch = extract_data$batch, verbose = verbose
        )
        mode <- "size-corrected"
    }
    if (mode == "size-corrected") { # Data has not been norm.
        if (verbose) message("- Normalizing the usage matrix.")
        norm_mat <- t(t(mat) / colSums(mat)) # Normalization
        means <- rowMeans(norm_mat, na.rm = TRUE) # Mean calculation
        usage_results <- tibble::tibble(
            !!var_name := names(means),
            mean_usage_across_conditions = as.numeric(means)
        )
    } else if (mode == "long-format") { # Data has been size-corrected & norm.
        usage_results <- data %>%
            dplyr::group_by(.data[[var_name]]) %>% # Group by features
            dplyr::summarize(mean_usage_across_conditions = mean( # Compute mean
                .data[["usage"]],
                na.rm = TRUE
            ))
        if (verbose) message("- The data is in a proper long format.")
    }
    CheckDataFrame(usage_results, required_names = FALSE)
    if (verbose) message("B . COMPLETED")
    if (is_object) { # Update the input tTEscanR object
        data <- UpdateObject(
            object = data, meta.data = usage_results, verbose = FALSE,
            meta.data.ids = paste(assay, "MeanUsage", sep = "_")
        )
        return(data)
    } else {
        return(usage_results) # Returns a tibble with: feature & mean usage
    }
}

GeneralChecksMeanUsage <- function(data, assay, metadata, batch, verbose) {
    assay_map <- list( # Links data to retrieve
        CodonUsage = "codon", AnticodonUsage = "anticodon",
        AADemand = "AA", AASupply = "AA"
    )
    if (verbose) message("A . Evaluating the input.")
    is_object <- inherits(data, "tTEscanR_Object")
    if (is_object) { # Dealing with a tTEscanR object

        if (verbose) message("- The input is a tTEscanR object.")
        if (!assay %in% names(assay_map)) {
            stop(
                "Invalid 'assay'.\n", "Please specify a valid 'assay' input ",
                "parameter: CodonUsage, AnticodonUsage, AADemand, AASupply."
            )
        }

        IsInObject(
            object = data, slot = "assays", section = assay, verbose = FALSE
        )
        IsInObject(
            object = data, slot = "meta.data",
            section = "CorrectionFactor", verbose = FALSE
        )
        IsInObject(
            object = data, slot = "meta.data",
            section = "ConditionsLabels", verbose = FALSE
        )

        raw_mat <- getAssay(data, assay)
        metadata <- getMetadata(data, "ConditionsLabels")
        batch <- getMetadata(data, "CorrectionFactor")
        var_name <- assay_map[[assay]]
    } else { # Dealing with a dataset
        raw_mat <- data
        var_name <- "feature"
    }

    if (verbose) message("A . COMPLETED")
    return(list(
        raw_mat = raw_mat, metadata = metadata, batch = batch,
        var_name = var_name
    ))
}

ModeRaw <- function(id_col, metadata, raw_mat, batch, verbose) {
    if (verbose) {
        message(
            "- Size correcting the matrix to account for sequencing depth."
        )
    }

    if (!is.null(id_col)) {
        if (!id_col %in% colnames(metadata)) is_col <- NULL
        if (verbose) {
            message(
                "The 'id_col' was not found in the 'metadata'. ",
                "Automatic detection will be implemented."
            )
        }
    }

    # Filter the data and the metadata based on their matching entries
    filtered <- FilterByMetadata(
        data = raw_mat, metadata = metadata, id_col = id_col
    )
    raw_mat <- ComputeSizeCorrection(
        data = filtered[[1]], verbose = FALSE,
        metadata = filtered[[2]], batch = batch
    )

    return(raw_mat)
}

ComputeMetricsCodonUsage <- function(
    codon_usage, codon_freq, metadata, id_col = NULL, batch = NULL,
    corr_method, verbose
) {
    # CODON EXONIC BACKGROUND
    if (verbose) message("- Computing the codon exonic background.")
    codon_exonic_back <- ComputeExonicBackground(data = codon_freq)

    if (is.null(codon_exonic_back)) {
        stop(
            "The codon exonic background could not be computed.\n",
            "Failure in ComputeExonicBackground()."
        )
    }

    # MEAN CODON USAGE
    if (verbose) message("- Computing the mean codon usage.")
    mean_codon_usage <- ComputeMeanUsage(
        data = codon_usage, metadata = metadata,
        id_col = id_col, batch = batch, verbose = FALSE
    )
    if (is.null(mean_codon_usage)) {
        stop(
            "The mean codon usage could not be computed.\n",
            "Failure in ComputeMeanUsage()."
        )
    }
    # CORRELATION BACKGROUND-MEAN
    if (verbose) {
        message(
            "- Computing the correlation between the mean codon ",
            "usage and the codon exonic background."
        )
    }
    corr_back_mean <- ComputeCorrelationBackground(
        mean = mean_codon_usage, background = codon_exonic_back,
        corr_method = corr_method, verbose = verbose
    )
    if (is.null(corr_back_mean)) {
        stop(
            "The correlation mean-background could not be computed.\n",
            "Failure in ComputeCorrelationBackground()."
        )
    }
    return(list( # Store all the metrics in a named list
        CodonExonicBackground = codon_exonic_back,
        MeanCodonUsage = mean_codon_usage,
        MeanCodonCorr = corr_back_mean
    ))
}

ConsistencyWithCodonFreq <- function(
    data, codon_freq, species, verbose = FALSE
) {
    if (verbose) {
        message(
            "\n------------------------------\n",
            "A . Assessing the 'codon_freq' matrix."
        )
    }

    # Checks the user-defined codon_freq or loads a default table if possible
    codon_freq <- CheckCodonFreqTable(data = codon_freq, species = species)
    if (verbose) {
        message(
            "A . COMPLETED\n",
            "B . Evaluating the consistency across gene annotations."
        )
    }

    gene_annot <- CheckGeneAnnotation(
        vector1 = colnames(codon_freq),
        vector2 = rownames(data), verbose = verbose
    )
    if (verbose) message("B . COMPLETED\n", "------------------------------\n")

    # Return the codon (rows) frequency per gene (columns) table
    return(codon_freq)
}

CheckNames_tRNA <- function(gene_names) {
    if (is.null(gene_names)) stop("No tRNA genes found in row names.")

    # Check the actual content of the tRNA gene label
    # Accepts both tRNA gene or isoacceptor names
    pattern <- "^tRNA-[A-Za-z]{3,4}-[ATGC]{3}(-[0-9]+-[0-9]+)?$"
    invalid_rows <- grep(pattern, gene_names, invert = TRUE)

    # Report if there is any rowname that does not follow the requirements above
    if (length(invalid_rows) > 0) {
        stop(
            "Inconsistent tRNA gene format.\n",
            "Expected format: tRNA-Asn-GTT-5-1\n",
            "1st invalid format: ", gene_names[invalid_rows[1]], "\n",
            "Invalid rows: ", paste(invalid_rows, collapse = ", ")
        )
    }
}
