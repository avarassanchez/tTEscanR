#' Compute Codon Usage from mRNA Gene Expression Data
#' @description
#' Estimates \strong{codon usage profiles} based on gene-level mRNA expression
#' data. It optionally accepts pre-computed codon frequency tables or uses
#' internally generated default tables when not provided. When enabled, it can
#' evaluate the correlation between background codon composition and observed
#' mean codon usage.
#'
#' @param object A \code{\link[MultiAssayExperiment]{MultiAssayExperiment}}
#'     containing a mRNA assay.
#' @param codon_freq Optional; a user-provided codon frequency-per-gene table.
#'     If necessary, it can be computed using \code{\link{getCodonFreq}}.
#' @param species Optional; either \code{"hg38"} (human) or \code{"mm39"}
#'     (mouse) to load the default settings. Required if \code{codon_freq} is
#'     not provided.
#' @param reduce Numeric; a scaling factor used to normalize large expression
#'     values that exceed R's handling capacity. Defaults to 100.
#' @param additional_metrics Logical; if \code{TRUE}, computes: (i) codon exonic
#'     background, (ii) mean codon usage, and (iii) correlation between the
#'     previous metrics. Defaults to \code{TRUE}.
#' @param corr_method A correlation method accepted by \code{\link{cor}}.
#'     Required if \code{additional_metrics} is \code{TRUE}. Defaults to
#'     \code{"spearman"}.
#' @param overwrite Logical; if \code{TRUE}, overwrites any existing assay and
#'     metadata in the \code{object}. Defaults to \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#'     Defaults to \code{TRUE}.
#'
#' @return An updated \code{\link[MultiAssayExperiment]{MultiAssayExperiment}}
#'     containing a new layer of information in the \code{assays}
#'     representing the \code{"CodonUsage"}. Additional computations will be
#'     stored in the \code{param} slot as \code{"CodonUsage_AdditionalMetrics"}.
#' @export
#'
#' @examples
#' data("default_tTEscanR_mRNA_data", package = "tTEscanR")
#' data("default_tTEscanR_metadata", package = "tTEscanR")
#'
#' tTEscanR_obj <- createObject(
#'     counts = default_tTEscanR_mRNA_data,
#'     assay = "mRNA",
#'     meta_data = default_tTEscanR_metadata,
#'     sample_id = "conditions",
#'     params = list("CorrectionFactor" = "tissue")
#' )
#' tTEscanR_obj <- computeCodonUsage(
#'     object = tTEscanR_obj, species = "hg38",
#'     additional_metrics = FALSE, reduce = 1000
#' )
computeCodonUsage <- function(object, codon_freq = NULL, species = NULL,
    additional_metrics = TRUE, reduce = 100, corr_method = "spearman",
    overwrite = FALSE, verbose = TRUE) {
    get_data <- generalChecksUsage(
        step = "codon", section = "mRNA", object = object,
        codon_freq = codon_freq, species = species, verbose = verbose
    )
    data <- get_data$mRNA
    codon_freq <- get_data$codon_freq
    if (verbose) {
        message("1 . COMPLETED")
        message("2 . Retrieving common mRNAs (mRNA and codon frequency data).")
    }
    common <- intersect(colnames(codon_freq), rownames(data))
    if (length(common) == 0) stop("No common mRNAs in mRNA & codon frequency.")
    if (verbose) {
        message("2 . COMPLETED\n", "3 . Computing the codon usage matrix.")
        message("- Filtering the datasets based on the common mRNAs.")
    }
    filt_codon_freq <- as.matrix(codon_freq[, common, drop = FALSE])
    codon_usage <- filt_codon_freq %*% as.matrix(data[common, , drop = FALSE])
    checkDataFrame(data = as.matrix(codon_usage))
    codon_usage <- checkIntegerLength(
        data = codon_usage, reduce = reduce, verbose = verbose
    )
    if (verbose) message("3 . COMPLETED")
    count <- 4
    if (isTRUE(additional_metrics)) {
        object <- getAdditionalMetrics(
            object = object, overwrite = overwrite, codon_usage = codon_usage,
            codon_freq = filt_codon_freq, corr = corr_method, verbose = verbose
        )
        count <- 5
    }
    if (verbose) message(count, ". Updating the tTEscanR object.")
    object <- updateObject(
        object = object, counts = codon_usage, assay = "CodonUsage",
        params = list("mRNAsInCommon" = common), verbose = FALSE,
        overwrite = overwrite
    )
    if (verbose) {
        message(count, ". COMPLETED")
        message("--- The codon usage has been successfully computed ---\n")
    }
    return(object)
}

getAdditionalMetrics <- function(object, overwrite, codon_usage, codon_freq,
    corr, verbose) {
    ## Check the correction factor
    checkParamPresent(object = object, param_name = "CorrectionFactor")
    batch <- S4Vectors::metadata(object)[["CorrectionFactor"]]

    ## Check the metadata in the object
    meta <- MultiAssayExperiment::colData(object)
    if (S4Vectors::isEmpty(meta)) {
        stop("The 'object' does not contain metadata.")
    }
    checkDataFrame(data = meta)

    if (!(batch %in% colnames(meta))) {
        stop("The correction factor was not found in the metadata.")
    }
    if (verbose) {
        message(
            "- The metadata and correction factor have",
            " been properly loaded.\n", "4 . Computing the additional metrics."
        )
    }

    additional_metrics <- computeMetricsCodonUsage(
        codon_usage = codon_usage, codon_freq = codon_freq, metadata = meta,
        corr_method = corr, batch = batch, verbose = verbose
    )
    if (verbose) {
        message(
            "- Adding the CodonUsage_AdditionalMetrics to the ",
            "meta.data of the tTEscanR object meta.data."
        )
    } # Nested list
    object <- updateObject(
        object = object,
        params = list("CodonUsage_AdditionalMetrics" = additional_metrics),
        overwrite = overwrite, verbose = FALSE
    )
    if (verbose) message("4 . COMPLETED")
    return(object)
}

#' Compute Anticodon Usage from tRNA Gene Expression Data
#' @description
#' Calculates \strong{anticodon usage profiles} from tRNA gene
#' abundance data. It summarizes the abundance of tRNAs by their anticodon
#' identity, which can be used to estimate the tRNA supply landscape. The tRNA
#' gene names need to be properly annotated for proper recognition.
#'
#' @param object A \code{\link[MultiAssayExperiment]{MultiAssayExperiment}}
#'     containing a tRNA assay.
#' @param overwrite Logical; if \code{TRUE}, overwrites any existing assay and
#'     metadata in the \code{object}. Defaults to \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#'     Defaults to \code{TRUE}.
#'
#' @return An updated \code{\link[MultiAssayExperiment]{MultiAssayExperiment}}
#'     containing a new layer of information in the \code{assays}
#'     representing the \code{"AnticodonUsage"}.
#' @export
#'
#' @examples
#' data("default_tTEscanR_tRNA_data", package = "tTEscanR")
#'
#' tTEscanR_obj <- createObject(
#'     counts = default_tTEscanR_tRNA_data,
#'     assay = "tRNA"
#' )
#' tTEscanR_obj <- computeAnticodonUsage(object = tTEscanR_obj)
computeAnticodonUsage <- function(object, overwrite = FALSE, verbose = TRUE) {
    extract_tRNA_data <- generalChecksUsage(
        step = "anticodon", section = "tRNA", object = object, verbose = verbose
    )
    tRNA_data <- extract_tRNA_data$tRNA_data
    tRNA_genes <- extract_tRNA_data$tRNA_genes
    if (verbose) {
        message("1 . COMPLETED")
        message("2 . Extracting the anticodons of each tRNA gene.")
    }
    anticodons <- sub("^[^-]+-[^-]+-([^-]+).*", "\\1", tRNA_genes)
    unique_anticodons <- sort(unique(anticodons))
    if (verbose) {
        message(
            "2 . COMPLETED\n",
            "3 . Pooling the counts from each tRNA gene with common anticodons."
        )
    }
    ## Initialize mapping matrix
    map_factor <- factor(anticodons, levels = unique_anticodons)
    i_indices  <- as.integer(map_factor)  # Row positions (anticodons)
    j_indices  <- seq_along(anticodons)   # Column positions (tRNA genes)
    M <- Matrix::sparseMatrix(
        i = i_indices, j = j_indices,
        x = 1,
        dims = c(length(unique_anticodons), length(anticodons)),
        dimnames = list(unique_anticodons, colnames(tRNA_genes))
    )
    anticodon_usage <- as.matrix(M %*% tRNA_data) # Multiply tRNA same anticodon

    checkDataFrame(data = anticodon_usage) # anticodons x conditions
    if (verbose) message("3 . COMPLETED\n", "4 . Updating the tTEscanR object.")
    object <- updateObject(
        object = object, counts = anticodon_usage, assay = "AnticodonUsage",
        verbose = FALSE, params = list("tRNAsAnticodons" = anticodons),
        overwrite = overwrite
    )
    if (verbose) {
        message(
            "4 . COMPLETED\n --- The anticodon usage has been successfully ",
            "computed ---\n"
        )
    }
    return(object)
}

generalChecksUsage <- function(step, section, object, codon_freq = NULL,
    species = NULL, verbose) {
    if (verbose) message("\n--- Computation of the ", step, " usage ---")
    if (verbose) message("1 . Checking the format of the input data.")
    checkObject(object = object, verbose = verbose)

    ## Check that the data is in object with proper format
    checkAssayPresent(object = object, assay_name = section)
    data <- SummarizedExperiment::assay(object, section)

    if (step == "anticodon") {
        data <- as.matrix(data)
        tRNA_genes <- rownames(data)
        checkNamestRNA(gene_names = tRNA_genes)
    }
    checkDataFrame(data = data)
    if (verbose) message("- The ", section, " assay has been properly loaded.")

    if (step == "codon") {
        codon_freq <- generalChecksCodonUsage(
            data = data, codon_freq = codon_freq, species = species,
            verbose = verbose
        )
    }
    if (step == "codon") {
        if (verbose) {
            message("- The codon frequency table has been properly loaded.")
        }
        return(list(mRNA = data, codon_freq = codon_freq))
    }

    if (step == "anticodon") {
        if (verbose) message("- The tRNA genes with a suitable format.")
        return(list(tRNA_data = data, tRNA_genes = tRNA_genes))
    }
}

generalChecksCodonUsage <- function(data, codon_freq, species, verbose) {
    ## Loading codon_freq & check consistency with mRNA_data
    available_formats <- c("dgCMatrix", "data.frame", "matrix")
    if (!is.null(codon_freq)) {
        if (!any(class(codon_freq) %in% available_formats)) {
            stop(
                "Wrong 'codon_freq' format.\n Supported formats: dgCMatrix, ",
                "data.frame, matrix"
            )
        }
    }

    codon_freq_table <- consistencyWithCodonFreq(
        data = data, codon_freq = codon_freq, species = species,
        verbose = verbose
    )

    return(codon_freq_table)
}
