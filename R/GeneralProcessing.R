#' Selection of the Optimal tRNA Cut Cutoff
#'
#' @param data tRNA gene expression count \code{matrix} with tRNA genes as rows
#' and conditions as columns.
#' @param num_iter Numeric; value to select the number of iterations to perform
#' in order to determine the optimal cutoff. Defaults to 1000.
#' @param cutoffs_limits Minimum and maximum values to test to search for the
#' optimal tRNA cuts threshold. Defaults to c(50, 10000).
#' @param compute_aa Logic; if \code{TRUE}, computes the amino acid supply,
#' otherwise only considers the anticodon usage. Defaults to \code{FALSE}.
#' @param generate_plot Logic; if \code{TRUE}, generates a correlation plot.
#' Defaults to \code{TRUE}.
#' @param slope_threshold Numeric; value to consider for the determination of
#' the correlation stability. Defaults to 0.001.
#' @param rho_threshold Numeric; value to consider for the determination of the
#' correlation strength. Defaults to 0.95.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#' Defaults to \code{TRUE}.
#'
#' @returns Table with the optimal cutoff at the anticodon isoacceptor and amino
#' acid isotype.
#' @export
#'
#' @examples
#' data(default_tTEscanR_tRNA_data)
#' optimal_tRNA_cutoffs <- tRNASetCutoff(
#'     data = default_tTEscanR_tRNA_data,
#'     generate_plot = FALSE, num_iter = 100,
#'     cutoffs_limits = c(3000, 5000)
#' )
tRNASetCutoff <- function(
    data, num_iter = 1000, cutoffs_limits = c(50, 10000), generate_plot = TRUE,
    slope_threshold = 0.001, rho_threshold = 0.95, compute_aa = FALSE,
    verbose = TRUE
) {
    if (verbose) message("1 . Computing reference tTEscanR object.")
    ref <- ReferenceObject(data = data, compute_aa = compute_aa)
    anticodon <- ref$anticodon
    supply <- ref$supply
    col_sums <- colSums(data)
    data <- TransformCounts(data)
    if (verbose) message("1 . COMPLETED\n2 . Extracting the potential cutoffs.")
    cuts <- c()
    retrieved_conditions <- c()
    for (i in seq(cutoffs_limits[1], cutoffs_limits[2], by = 1)) {
        num_cond <- sum(col_sums > i)
        if (!num_cond %in% retrieved_conditions) {
            retrieved_conditions <- c(retrieved_conditions, num_cond)
            cuts <- c(cuts, i)
        }
    }
    if (verbose) {
        message(
            "- Cutoffs retrieved: ", paste(cuts, collapse = ", "), "\n2 . ",
            "COMPLETED\n3 . Compute correlations and determine optimal cutoff."
        )
    }
    cutoff_res <- RunIterations(
        num_iter = num_iter, data = data, anticodon = anticodon, cuts = cuts,
        supply = supply, slope = slope_threshold, rho = rho_threshold)
    if (verbose) message("3 . COMPLETED")
    cutoff_res <- dplyr::bind_rows(cutoff_res, .id = "iteration")
    if (isTRUE(generate_plot)) {
        if (verbose) message("4 . Generating plots.")
        cor_results <- ComputeCorrelations(
            data = data, ref_anticodon = anticodon, ref_supply = supply,
            cutoffs = cuts
        )
        corr_plot <- CorrelationCutoffPlot(data = cor_results)
        hist_plot <- SelectionCutoffPlot(data = cutoff_res)
        if (verbose) message("4 . COMPLETED")
        return(list(
            optimal_cutoff = cutoff_res,
            correlation_plot = corr_plot, histogram_plot = hist_plot
        ))
    } else {
        return(list(optimal_cutoff = cutoff_res))
    }
}

#' Generate a tRNA expression matrix
#'
#' @param data \code{SummarizedExperiment}, \code{ChromatinAssay}, or
#' \code{SeuratObject}.
#' @param confidence_set Either a file path to the tRNA annotations
#' (confidence set file from gtRNAdb), or a GRanges object. Contains the set of
#' high confidence tRNA genes
#' @param tRNA_name_map Optional; a \code{data.frame} with the tRNA gene names
#' linked to \code{confidence_set}. Contains two columns: tRNAscan-SE and
#' GtRNAdb gene ids.
#' @param species Optional; either \code{"hg38"} (human) or \code{"mm39"}
#' (mouse) to load the default \code{confidence_set}
#' @param flanking_region Integer; number of nucleotides that form the flanking
#' region of each tRNA. Defaults to 100.
#' @param assay Optional; a character string specifying the name of the assay to
#' retrieve from \code{data} if it is a \code{SeuratObject}. Defaults to
#' \code{"peaks"}.
#' @param name_sep A string delimiter to format the tRNA gene names in the
#' output matrix. Defaults to \code{c("-", "-")}.
#' @param save Logical; if \code{TRUE} stores the generated tRNA matrix into a
#' file.
#' @param out_name Optional; name for the saved plot (if \code{save} specified).
#' @param out_directory Optional; path to the directory where the plot will be
#' saved (if \code{save} specified).
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#' Defaults to \code{TRUE}.
#'
#' @return Sparse matrix of tRNA counts (tRNAs x cells)
#' @export
#'
tRNAGetMatrix <- function(
    data, assay = "peaks", confidence_set = NULL, tRNA_name_map = NULL,
    species = NULL, flanking_region = 100, name_sep = c("-", "-"), save = TRUE,
    out_name = NULL, out_directory = NULL, verbose = TRUE
) {
    if (verbose) message("1 . Importing the high-confidence tRNA annotations.")
    tRNA_granges <- Get_tRNARanges(
        conf = confidence_set, flanking_region = flanking_region,
        map = tRNA_name_map, species = species, verbose = verbose
    )
    if (verbose) message(
        "1 . COMPLETED\n", "2 . Filtering tRNA genes with unknown annotations."
    )
    tRNA_granges <- Filter_tRNARanges(
        granges = tRNA_granges, flanking_region = flanking_region,
        map = tRNA_name_map, species = species
    )
    if (verbose) message(
        "2 . COMPLETED\n", "3 . Finding overlaps and aggregating tRNA counts."
    )
    tRNA_matrix <- GenerateMatrix(
        data = data, tRNA_granges = tRNA_granges, assay = assay,
        name_sep = name_sep, verbose = verbose
    )
    if (verbose) message("3 . COMPLETED\n")
    if (isTRUE(save)) {
        if (verbose) message("4 . Exporting tRNA expression matrix.")
        output_file <- GetOutputName(
            action = "file", out_name = out_name, out_directory = out_directory,
            save_format = "rds", verbose = verbose
        )
        saveRDS(tRNA_matrix, output_file)
        if (verbose) message("4 . COMPLETED\n")
    }
    return(tRNA_matrix) # tRNA gene matrix to us as input in tTEscanR
}

Get_tRNARanges <- function(conf, flanking_region, map, species, verbose) {
    # No annotations provided, using default data (if possible)
    if (is.null(conf)) {
        if (is.null(species)) stop(
            "No 'confidence_set' provided and no 'species' specified.\n",
            "Specify 'hg38' for human or 'mm39' for mouse default data."
        )
        granges <- SelectDefaultData(
            species = species, action = "confidence_set"
        )
    # The conf parameter is already a GRanges object
    } else if (inherits(conf, "GRanges")) {
        if (verbose) message("The input is already a GRanges object.")
        granges <- conf
    # The conf parameter is the file to the confidence-set file
    } else if (is.character(conf)) { # Input is a path
        if (file.exists(conf) || dir.exists(conf)) {
            if (verbose) message("Importing tRNA annotations from: ", conf)
            granges <- tRNAscanImport::import.tRNAscanAsGRanges(conf)
        } else stop("The provided path does not exist: ", conf)
    # Incorrect format
    } else {
        stop(
            "The 'confidence_set' argument must be a valid file path ",
            "(string) or a GRanges object."
        )
    }
    return(granges)
}

Filter_tRNARanges <- function(granges, flanking_region, map, species) {
    granges <- granges[granges$tRNA_anticodon != "NNN"] # Filter unknown
    granges <- GenomicRanges::trim(granges + flanking_region)

    if (is.null(map)) {
        if (is.null(species)) stop(
            "No 'tRNA_name_map' provided and no 'species' specified.\n",
            "Specify 'hg38' for human or 'mm39' for mouse default data."
        )
        map <- SelectDefaultData(species = species, action = "tRNA_map")

    }

    tRNA_id <- paste0(GenomicRanges::seqnames(granges), ".trna", granges$no)
    map_vec <- stats::setNames(map$GtRNAdb_id, map$tRNAscan.SE_id)
    granges$gene_name <- map_vec[tRNA_id]

    granges$gene_biotype <- 'tRNA'

    return(granges)
}

GenerateMatrix <- function(data, tRNA_granges, assay, name_sep, verbose) {
    if (inherits(data, "SummarizedExperiment")) {

        counts_peak_matrix <- SummarizedExperiment::assay(data, "counts")
        peak_ranges <- SummarizedExperiment::rowRanges(data)
        match <- GenomicRanges::findOverlaps(tRNA_granges, peak_ranges)
        M <- Matrix::sparseMatrix(
            i = S4Vectors::queryHits(match), j = S4Vectors::subjectHits(match),
            x = 1, dims = c(length(tRNA_granges), length(peak_ranges))
        )
        tRNA_matrix <- M %*% as.matrix(counts_peak_matrix)
        rownames(tRNA_matrix) <- tRNA_granges$gene_name
        colnames(tRNA_matrix) <- colnames(counts_peak_matrix)

    } else if (inherits(data, c("Seurat", "ChromatinAssay"))) {
        chrom_assay <- if (inherits(data, "Seurat")) {
            data[[assay %||% "peaks"]] # Check the peaks identifier
        } else data

        if (is.null(chrom_assay)) stop("Requested assay not found in object.")
        frags <- Signac::Fragments(chrom_assay)

        if (length(frags) == 0) stop("No fragments found in the input assay.")
        all_tRNA_coords <- Signac::GRangesToString(tRNA_granges, sep = name_sep)

        tRNA_matrix <- Signac::FeatureMatrix(
            fragments = frags, features = tRNA_granges,
            cells = colnames(chrom_assay), sep = name_sep, verbose = verbose
        )

        matched_indices <- match(rownames(tRNA_matrix), all_tRNA_coords)
        rownames(tRNA_matrix) <- tRNA_granges$gene_name[matched_indices]

    } else {
        stop(
            "Unsupported 'data' type. Use SummarizedExperiment, Seurat, ",
            "or ChromatinAssay."
        )
    }

    return(tRNA_matrix)
}

#' Filter Out Conditions With Low tRNA Cuts
#' @description
#' This function filters a tRNA expression matrix by removing conditions
#' (columns) that fall below a specific total read count (\code{cutoff}). It
#' is useful for eliminating low-quality or poorly sequenced conditions that
#' may bias downstream analyses.
#'
#' @param data A \code{matrix} or \code{data.frame} of tRNA gene expression
#' data, with tRNA genes as rows and conditions as columns.
#' @param cutoff Numeric; minimum total number of tRNA cuts required to retain
#' a condition in \code{data}. Defaults to 5000.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#' Defaults to \code{TRUE}.
#'
#' @return A filtered \code{matrix} or \code{data.frame} with tRNAs below the
#' cutoff removed.
#' @export
#'
#' @examples
#' data(default_tTEscanR_tRNA_data)
#' tRNA_data_filtered <- tRNAFilterCuts(
#'     data = default_tTEscanR_tRNA_data,
#'     cutoff = 5000
#' )
tRNAFilterCuts <- function(data, cutoff = 5000, verbose = TRUE) {
    if (verbose) {
        message(
            "--- Filtering the tRNA gene abundance count matrix ---",
            "\n- Initial number of samples: ", ncol(data)
        )
    }

    keep_samples <- Matrix::colSums(data) >= cutoff
    if (!any(keep_samples)) {
        stop("No samples pass the cutoff (", cutoff, "). Adjust the parameter.")
    }

    data <- data[, keep_samples, drop = FALSE]
    if (verbose) {
        message("- Samples retained: ", ncol(data),
                "\n--- The tRNA gene abundance count matrix has been ",
                "successfully filtered ---"
        )
    }

    return(data) # Returns the filtered tRNA gene expression count matrix.
}

#' Annotate the tRNA genes from tRNA tags
#'
#' @param data A \code{matrix} with tRNA tags as rows.
#' @param tRNA_bed Path to the directory that contains the .bed file
#' @param flanking_region Numeric; number of bases to include expand the region
#' interrogated. Defaults to 100.
#'
#' @returns A \code{data} with the translated tRNA gene names.
#' @export
#'
tRNASetGenes <- function(data, tRNA_bed, flanking_region = 100) {
    tRNA_table <- utils::read.delim(
        tRNA_bed,
        header = FALSE,
        stringsAsFactors = FALSE
    )

    tRNA_table$V2 <- tRNA_table$V2 + 1 - flanking_region
    tRNA_table$V3 <- tRNA_table$V3 + flanking_region

    ref_regions <- paste(tRNA_table$V1, tRNA_table$V2, tRNA_table$V3, sep = "-")
    mapping <- stats::setNames(tRNA_table$V4, ref_regions)

    new_names <- mapping[rownames(data)]
    rownames(data) <- ifelse(is.na(new_names), rownames(data), new_names)

    return(data)
}

#' Combine large matrices
#' @description
#' This function efficiently combines individual matrices.
#'
#' @param ... A variable number of \code{matrix}.
#'
#' @return A single sparse matrix with as a combination of all the input
#' matrices.
#' @export
#'
#' @examples
#' df1 <- matrix(c(1, 0, 0, 2), nrow = 2, dimnames = list(
#'     c("geneA", "geneB"),
#'     c("s1", "s2")
#' ))
#' df2 <- matrix(c(3, 0, 0, 4), nrow = 2, dimnames = list(
#'     c("geneB", "geneC"),
#'     c("s2", "s3")
#' ))
#' merged_matrix <- MergeMatrices(df1, df2)
MergeMatrices <- function(...) {
    cnnew <- character()
    rnnew <- character()
    x <- numeric()
    i <- numeric()
    j <- numeric()

    # Iterates over each matrices passed
    for (df in list(...)) {
        df <- as.matrix(df)
        storage.mode(df) <- "numeric"

        cnold <- colnames(df)
        rnold <- rownames(df)

        cnnew <- union(cnnew, cnold)
        rnnew <- union(rnnew, rnold)

        cindnew <- match(cnold, cnnew)
        rindnew <- match(rnold, rnnew)

        non_zero <- which(df != 0, arr.ind = TRUE)
        i <- c(i, rindnew[non_zero[, 1]])
        j <- c(j, cindnew[non_zero[, 2]])
        x <- c(x, df[non_zero])
    }

    # Create sparse matrix - only the non-zero values are stored
    result_df <- Matrix::sparseMatrix(
        i = i, j = j, x = x, dims = c(length(rnnew), length(cnnew)),
        dimnames = list(rnnew, cnnew)
    )

    return(result_df)
}

#' Aggregates data by group
#' @description
#' This function calculates the row sums of a given matrix to combine columns
#' that share the same group.
#'
#' @param data A \code{matrix} with the features to group for as columns.
#' @param group_labels A \code{vector} with the metadata features to group
#' the columns in \code{data}
#'
#' @return A \code{matrix} with the conditions merged based on the metadata.
#' @export
#'
#' @examples
#' data <- data.frame(
#'     sample_1 = c(10, 5, 20), sample_2 = c(15, 8, 25),
#'     sample_3 = c(12, 6, 22), sample_4 = c(1, 2, 3),
#'     sample_5 = c(4, 5, 6), sample_6 = c(7, 8, 9)
#' )
#' rownames(data) <- c("gene_1", "gene_2", "gene_3")
#' groups <- c("cond_A", "cond_A", "cond_A", "cond_B", "cond_B", "cond_B")
#' data_combined <- GroupConditions(data = data, group_labels = groups)
GroupConditions <- function(data, group_labels) {
    # Check the dimensions
    if (ncol(data) != length(group_labels)) {
        stop(
            "Dimension mismatch: The number of columns in 'data' (",
            ncol(data), ") must match the length of 'group_labels' (",
            length(group_labels), ")."
        )
    }

    # Check for missing (NA) labels
    if (any(is.na(group_labels))) {
        stop(
            "NAs found in 'group_labels'. Rremove or handle missing grouping ",
            "information before running."
        )
    }

    # Check the order of the labels - Ensure group_labels contains factors
    if (is.factor(group_labels)) {
        groups <- droplevels(group_labels)
    } else {
        groups <- factor(group_labels)
    }

    if (inherits(data, "Matrix")) {
        M <- Matrix::sparse.model.matrix(~ 0 + groups)
        res <- data %*% M
        res <- as.matrix(res)
    } else {
        if (!is.matrix(data)) mat_data <- as.matrix(data) else mat_data <- data
        res <- t(rowsum(t(mat_data), group = groups, reorder = FALSE))
    }

    rownames(res) <- rownames(data)
    colnames(res) <- levels(groups)

    return(as.data.frame(res))
}

#' Transform the format of a table
#' @description
#' This function converts a \code{matrix} or \code{data.frame} into a tidy,
#' long-format \code{tibble}. Optionally normalizes the values.
#'
#' @param data A table to be converted. Supported formats: \code{matrix} or
#' \code{data.frame}.
#' @param normalize Logical; if \code{TRUE}, values are converted to relative
#' frequencies. Defaults to \code{FALSE}.
#' @param rownames_to_column A character string specifying the name of the new
#' column that will hold the former row names in \code{data}.
#' @param names_to A character string specifying the name of the new column
#' that will hold the former column names in \code{data}.
#' @param values_to A character string specifying the name of the new column
#' that will hold the corresponding values from the pivoted columns in
#' \code{data}.
#'
#' @return A tibble of the input \code{data}
#' @export
#'
#' @examples
#' data(default_tTEscanR_tRNA_data)
#' tRNA_long_format <- TransformFormat(
#'     data = default_tTEscanR_tRNA_data,
#'     normalize = FALSE,
#'     rownames_to_column = "tRNA_genes",
#'     names_to = "condition",
#'     values_to = "abundance"
#' )
#' tTEobj <- CreateObject(
#'     counts = default_tTEscanR_tRNA_data,
#'     assay = "tRNA"
#' )
#' tTEobj <- ComputeAnticodonUsage(object = tTEobj)
#' anticodon_long_format <- TransformFormat(
#'     data = getAssay(tTEobj, "AnticodonUsage"),
#'     normalize = TRUE,
#'     rownames_to_column = "anticodons",
#'     names_to = "condition", values_to = "usage"
#' )
TransformFormat <- function(
    data, normalize, rownames_to_column, names_to, values_to
) {
    # Vectorized normalization (if required)
    if (isTRUE(normalize)) data <- t(t(data) / colSums(data))

    # Ensure data is a matrix
    if (!is.matrix(data)) data <- as.matrix(data)

    num_row <- nrow(data)
    num_col <- ncol(data)

    # Transform the data into a long format
    long_format <- data.frame(rep.int(rownames(data), num_col),
        rep.int(colnames(data), rep.int(num_row, num_col)),
        as.vector(data),
        stringsAsFactors = FALSE
    )
    colnames(long_format) <- c(rownames_to_column, names_to, values_to)

    CheckDataFrame(long_format) # Evaluate that the data is properly defined

    # Returns the processed data (long-format and normalized if applicable).
    return(long_format)
}
