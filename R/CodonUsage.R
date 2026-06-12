#' Compute Codon Usage from mRNA Gene Expression Data
#' @description
#' This function estimates \strong{codon usage profiles} based on gene-level
#' mRNA expression data stored in a \code{tTEscanR_Object}. It optionally
#' accepts pre-computed codon frequency tables or uses internally generated
#' default tables when not provided. When enabled, it can evaluate the
#' correlation between background codon composition and observed mean codon
#' usage. If the additional metrics are to be computed the input
#' \code{tTEscanR_Object} needs to The default \code{codon_freq} were built
#' using the canonical filter to select one transcript if several were
#' available for the same gene.
#'
#' @param object A \code{tTEscanR_Object} containing a mRNA assay.
#' @param codon_freq Optional; a user-provided codon frequency-per-gene table.
#' If necessary, it can be computed using \code{\link{GetCodonFreq}}.
#' @param species Optional; either \code{"hg38"} (human) or \code{"mm39"}
#' (mouse) to load the default settings. Required if \code{codon_freq} is not
#' provided.
#' @param reduce Numeric; a scaling factor used to normalize large expression
#' values that exceed R's handling capacity. Defaults to 100.
#' @param additional_metrics Logical; if \code{TRUE}, computes: (i) codon exonic
#' background, (ii) mean codon usage, and (iii) correlation between the previous
#' metrics. Defaults to \code{TRUE}.
#' @param corr_method A correlation method accepted by \code{\link{cor}}.
#' Required if \code{additional_metrics} is \code{TRUE}. Defaults to
#' \code{"spearman"}.
#' @param overwrite Logical; if \code{TRUE}, overwrites any existing assay and
#' metadata in the \code{object}. Defaults to \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#' Defaults to \code{TRUE}.
#'
#' @return An updated \code{tTEscanR_Object} containing a new layer of
#' information \code{"CodonUsage"} in the \code{assays} slot representing the
#' codon usage. Additional computations will be stored in the \code{meta.data}
#' slot as \code{"CodonUsage_AdditionalMetrics"}.
#' @export
#'
#' @examples
#' data(default_tTEscanR_mRNA_data, default_tTEscanR_metadata)
#' tTEscanR_obj <- CreateObject(
#'     counts = default_tTEscanR_mRNA_data,
#'     assay = "mRNA",
#'     meta.data = list(default_tTEscanR_metadata, "tissue"),
#'     meta.data.ids = list("ConditionsLabels", "CorrectionFactor")
#' )
#' tTEscanR_obj <- ComputeCodonUsage(
#'     object = tTEscanR_obj, species = "hg38",
#'     additional_metrics = FALSE, reduce = 1000
#' )
ComputeCodonUsage <- function(
    object, codon_freq = NULL, species = NULL, additional_metrics = TRUE,
    reduce = 100, corr_method = "spearman", overwrite = FALSE, verbose = TRUE
) {
    get_data <- GeneralChecksUsage(
        step = "codon", section = "mRNA", object = object,
        codon_freq = codon_freq, species = species, verbose = verbose
    )
    data <- get_data$mRNA
    codon_freq <- get_data$codon_freq
    if (verbose) {
        message(
            "1 . COMPLETED\n", "2 . Retrieving common mRNAs (gene expression ",
            "and codon frequency table)."
        )
    }
    common <- intersect(colnames(codon_freq), rownames(data))
    if (length(common) == 0) stop("No common mRNAs in mRNA & codon frequency.")
    if (verbose) {
        message("2 . COMPLETED\n", "3 . Computing the codon usage matrix.")
        message("- Filtering the datasets based on the common mRNAs.")
    }
    filt_codon_freq <- as.matrix(codon_freq[, common, drop = FALSE])
    codon_usage <- filt_codon_freq %*% as.matrix(data[common, , drop = FALSE])
    CheckDataFrame(data = as.matrix(codon_usage))
    codon_usage <- CheckIntegerLength(
        data = codon_usage, reduce = reduce, verbose = verbose
    )
    if (verbose) message("3 . COMPLETED")
    count <- 4
    if (isTRUE(additional_metrics)) {
        object <- GetAdditionalMetrics(
            object = object, overwrite = overwrite, codon_usage = codon_usage,
            codon_freq = filt_codon_freq, corr = corr_method, verbose = verbose
        )
        count <- 5
    }
    if (verbose) message(as.character(count), ". Updating the tTEscanR object.")
    object <- UpdateObject(
        object = object, counts = codon_usage, assay = "CodonUsage",
        meta.data = list(common), verbose = FALSE,
        meta.data.ids = list("mRNAsInCommon"), overwrite = overwrite
    )
    if (verbose) message(
        as.character(count), ". COMPLETED\n",
        "--- The codon usage has been successfully computed ---\n"
    )
    return(object) # tTEscanR object was validated in UpdateObject()
}

GetAdditionalMetrics <- function(
    object, overwrite, codon_usage, codon_freq, corr, verbose
) {
    IsInObject(
        object = object, slot = "meta.data",
        section = "ConditionsLabels", verbose = FALSE
    )
    meta <- getMetadata(object, "ConditionsLabels")
    CheckDataFrame(data = meta)
    IsInObject(
        object = object, slot = "meta.data",
        section = "CorrectionFactor", verbose = FALSE
    )
    batch <- getMetadata(object, "CorrectionFactor")
    if (!(batch %in% colnames(meta))) {
        stop("The correction factor was not found in the metadata.")
    }
    if (verbose) {
        message(
            "- The 'ConditionsLabels' and 'CorrectionFactor' have",
            " been properly loaded.\n", "4 . Computing the additional metrics."
        )
    }
    id_col <- IsInObject(
        object = object, slot = "meta.data", section = "DataMetadataIndex",
        verbose = FALSE, compute_assay = TRUE
    )
    id <- if (isFALSE(id_col)) NULL else {
        getMetadata(object, "DataMetadataIndex")
    }
    additional_metrics <- ComputeMetricsCodonUsage(
        codon_usage = codon_usage, codon_freq = codon_freq, metadata = meta,
        id_col = id, corr_method = corr, batch = batch, verbose = verbose
    )
    if (verbose) {
        message(
            "- Adding the CodonUsage_AdditionalMetrics to the ",
            "meta.data of the tTEscanR object meta.data."
        )
    } # Nested list
    object <- UpdateObject(
        object = object, main_name = "CodonUsage_AdditionalMetrics",
        meta.data = additional_metrics, overwrite = overwrite, verbose = FALSE
    )
    if (verbose) message("4 . COMPLETED")
    return(object)
}

#' Compute Anticodon Usage from tRNA Gene Expression Data
#' @description
#' This function calculates \strong{anticodon usage profiles} from tRNA gene
#' expression data stored in a \code{tTEscanR_Object}. It summarizes the
#' expression of tRNAs by their anticodon identity, which can be used to
#' estimate the tRNA supply landscape. The tRNA gene names need to be properly
#' annotated for proper recognition. Expected format: tRNA-Asn-GTT-5-1.
#'
#' @param object A \code{tTEscanR_Object} containing a tRNA assay.
#' @param overwrite Logical; if \code{TRUE}, overwrites any existing assay and
#' metadata in the \code{object}. Defaults to \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#' Defaults to \code{TRUE}.
#'
#' @return An updated \code{tTEscanR_Object} containing a new layer of
#' information \code{"AnticodonUsage"} in the \code{assays} slot representing
#' the anticodon usage.
#' @export
#'
#' @examples
#' data(default_tTEscanR_tRNA_data)
#' tTEscanR_obj <- CreateObject(
#'     counts = default_tTEscanR_tRNA_data,
#'     assay = "tRNA"
#' )
#' tTEscanR_obj <- ComputeAnticodonUsage(object = tTEscanR_obj)
ComputeAnticodonUsage <- function(object, overwrite = FALSE, verbose = TRUE) {
    extract_tRNA_data <- GeneralChecksUsage(
        step = "anticodon", section = "tRNA", object = object, verbose = verbose
    )
    tRNA_data <- extract_tRNA_data$tRNA_data
    tRNA_genes <- extract_tRNA_data$tRNA_genes
    if (verbose){
        message(
            "1 . COMPLETED\n",
            "2 . Extracting the anticodons of each tRNA gene."
        )
    }
    anticodons <- vapply(strsplit(tRNA_genes, "-"), "[[", 3,
        FUN.VALUE = character(1)
    )
    unique_anticodons <- sort(unique(anticodons))
    if (verbose) {
        message(
            "2 . COMPLETED\n",
            "3 . Pooling the counts from each tRNA gene with common anticodons."
        )
    }
    # Initialize an mapping matrix - 1 to add the gene and 0 to ignore it
    map_factor <- factor(anticodons, levels = unique_anticodons)
    M <- Matrix::sparse.model.matrix(~ 0 + map_factor)
    M <- Matrix::t(M) # Set the anticodons as rows
    rownames(M) <- gsub("map_factor", "", rownames(M)) # Clean the rownames
    anticodon_usage <- as.matrix(M %*% tRNA_data) # Multiply tRNA same anticodon
    rownames(anticodon_usage) <- unique_anticodons
    colnames(anticodon_usage) <- colnames(tRNA_data)
    CheckDataFrame(data = anticodon_usage) # anticodons x conditions
    if (verbose) message("3 . COMPLETED\n", "4 . Updating the tTEscanR object.")
    object <- UpdateObject(
        object = object, counts = anticodon_usage, assay = "AnticodonUsage",
        verbose = FALSE, meta.data = anticodons,
        meta.data.ids = "tRNAsAnticodons", overwrite = overwrite
    )
    if (verbose) {
        message(
            "4 . COMPLETED\n --- The anticodon usage has been successfully ",
            "computed ---\n"
        )
    }
    return(object) # tTEscanR object was validated in UpdateObject()
}

GeneralChecksUsage <- function(
    step, section, object, codon_freq = NULL, species = NULL, verbose
) {
    if (verbose) message("\n--- Computation of the ", step, " usage ---")
    if (verbose) message("\n1 . Checking the format of the input data.")
    if (!inherits(object, "tTEscanR_Object")) {
        stop("'object' must be a tTEscanR object.")
    }
    if (verbose) message("- The input consists of a proper tTEscanR object.")

    IsInObject( # Check that the data is in object with proper format
        object = object, slot = "assays", section = section, verbose = FALSE
    )

    data <- getAssay(object, section)
    if (step == "anticodon") {
        data <- as.matrix(data)
        tRNA_genes <- rownames(data)
        CheckNames_tRNA(gene_names = tRNA_genes)
    }
    CheckDataFrame(data = data)
    if (verbose) message("- The ", section, " assay has been properly loaded.")

    if (step == "codon"){
        codon_freq <- GeneralChecksCodonUsage(
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

GeneralChecksCodonUsage <- function(data, codon_freq, species, verbose) {
    # Loading codon_freq & check consistency with mRNA_data
    available_formats <- c("dgCMatrix", "data.frame", "matrix")
    if (!is.null(codon_freq)) {
        if (!any(class(codon_freq) %in% available_formats)) {
            stop(
                "Wrong 'codon_freq' format.\n Supported formats: dgCMatrix, ",
                "data.frame, matrix"
            )
        }
    }

    codon_freq_table <- ConsistencyWithCodonFreq(
        data = data, codon_freq = codon_freq, species = species,
        verbose = verbose
    )

    return(codon_freq_table)
}
