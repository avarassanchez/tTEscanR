#' Selection of the Optimal tRNA Cut Cutoff
#'
#' @param data tRNA gene expression count \code{matrix} with tRNA genes as rows and conditions as columns.
#' @param num_iter Numeric; value to select the number of iterations to perform in order to determine the optimal cutoff. Defaults to 1000.
#' @param cutoffs_limits Minimum and maximum values to test to search for the optimal tRNA cuts threshold. Defaults to c(50, 10000).
#' @param compute_aa Logic; if \code{TRUE}, computes the amino acid supply, otherwise only considers the anticodon usage. Defaults to \code{FALSE}.
#' @param generate_plot Logic; if \code{TRUE}, generates a correlation plot. Defaults to \code{TRUE}.
#' @param slope_threshold Numeric; value to consider for the determination of the correlation stability. Defaults to 0.001.
#' @param rho_threshold Numeric; value to consider for the determination of the correlation strength. Defaults to 0.95.
#'
#' @returns Table with the optimal cutoff at the anticodon isoacceptor and amino acid isotype.
#' @export
#'
#' @examples
#' data(default_tTEscanR_tRNA_data)
#' optimal_tRNA_cutoffs <- Set_tRNACutoff(data = default_tTEscanR_tRNA_data,
#'                                        generate_plot = FALSE, num_iter = 100,
#'                                        cutoffs_limits = c(3000, 5000))

Set_tRNACutoff <- function(data, num_iter = 1000, cutoffs_limits = c(50, 10000), generate_plot = TRUE, slope_threshold = 0.001, rho_threshold = 0.95, compute_aa = FALSE) {

  message("1 . Computing reference tTEscanR object.")
  object_ref <- suppressMessages(Create_tTEscanR_Object(counts = data, assay = "tRNA", verbose = FALSE))
  object_ref <- suppressMessages(ComputeAnticodonUsage(object_ref, verbose = FALSE))
  if (isTRUE(compute_aa)) object_ref <- suppressMessages(ComputeAAUsage(object_ref, level = "supply", verbose = FALSE))

  calculate_reference <- function(assay_slot){
    value <- rowSums(assay_slot)
    return(value / sum(value))
  }

  ref_anticodon <- calculate_reference(object_ref@assays$AnticodonUsage)
  if (isTRUE(compute_aa)) ref_supply <- calculate_reference(object_ref@assays$AASupply) else ref_supply <- NULL

  data_long <- TransformCounts(data)
  # if (is.null(cutoffs)) cutoffs <- 1:seq_len(colnames(object_ref@assays$mRNA))
  message("1 . COMPLETED\n", "2 . Extracting the potential cutoffs.")

  cutoffs_tRNA <- seq(cutoffs_limits[1], cutoffs_limits[2], by = 1)
  cutoffs <- c()
  retrieved_conditions <- c()

  col_sums <- colSums(data)
  for (i in cutoffs_tRNA){
    num_cond <- sum(col_sums > i)
    if (!num_cond %in% retrieved_conditions){
      retrieved_conditions <- c(retrieved_conditions, num_cond)
      cutoffs <- c(cutoffs, i)
    }
  }

  message(paste("- Cutoffs retrieved:", paste(cutoffs, collapse = ", ")), "\n2 . COMPLETED\n", "3 . Compute correlations and determine optimal cutoff.")

  tRNA_cutoff_results_list <- list() # Store the optimal results of each iteration
  pb <- utils::txtProgressBar(min = 0, max = num_iter, style = 3) # Start a time counter to track the progress of the function execution
  for (i in 1:num_iter){ # Obtain the optimal cutoff per iteration
    out <- Iterate_tRNACutoff(data = data_long, anticodon = ref_anticodon, supply = ref_supply,
                              slope_threshold = slope_threshold, rho_threshold = rho_threshold, cutoffs = cutoffs)
    tRNA_cutoff_results_list[[i]] <- out
    utils::setTxtProgressBar(pb, i) # Increase the progress bar
  }
  close(pb)
  message("3 . COMPLETED")

  results_tRNA_cutoff_long <- dplyr::bind_rows(tRNA_cutoff_results_list, .id = "iteration") # Format the output table of optimal codons per each iteration

  if (isTRUE(generate_plot)){
    message("4 . Generating plots.")
    cor_results <- ComputeCorrelations(data = data_long, ref_anticodon = ref_anticodon, ref_supply = ref_supply, cutoffs = cutoffs)
    correlation_plot <- CorrelationCutoffPlot(data = cor_results)
    histogram_plot <- SelectionCutoffPlot(data = results_tRNA_cutoff_long)

    message("4 . COMPLETED")
    return(list(optimal_cutoff = results_tRNA_cutoff_long, correlation_plot = correlation_plot, histogram_plot = histogram_plot))
  } else {
    return(list(optimal_cutoff = results_tRNA_cutoff_long))
  }
}

#' Generate a tRNA expression matrix
#'
#' @param data \code{SummarizedExperiment}, \code{ChromatinAssay}, or \code{SeuratObject}.
#' @param confidence_set A directory path to tRNA annotations (confidence set) in .ss or .bed format.
#' @param tRNA_name_map Optional; a \code{data.frame} with the tRNA gene names linked to \code{tRNA_annotations}.
#' @param flanking_region Integer; number of nucleotides that form the flanking region of each tRNA. Defaults to 100.
#' @param assay Optional; a character string specifying the name of the assay to retrieve from \code{data} if it is a \code{SeuratObject}. Defaults to \code{"peaks"}.
#' @param name_sep A string delimiter to format the tRNA gene names in the output matrix. Defaults to \code{c("-", "-")}.
#' @param save Logical; if \code{TRUE} stores the generated tRNA matrix into a file.
#' @param out_name Optional; name for the saved plot (if \code{save} specified).
#' @param out_directory Optional; path to the directory where the plot will be saved (if \code{save} specified).
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return Sparse matrix of tRNA counts (tRNAs x cells)
#' @export
#'

Get_tRNAMatrix <- function(data, assay = "peaks", confidence_set, tRNA_name_map, flanking_region = 100, name_sep = c("-", "-"),
                           save = TRUE, out_name = NULL, out_directory = NULL, verbose = TRUE){

  ###
  # DESCRIPTION: This function takes a Seurat or Signac object from scATAC-seq data and the corresponding tRNA
  # reference gene annotations to generate a tRNA gene expression count matrix using the proper Signac built-in functions.
  ###

  # Common tRNA annotation preprocessing
  if (verbose) message("1 . Importing and filtering tRNA annotations refernece file.")
  tRNA_granges <- tRNAscanImport::import.tRNAscanAsGRanges(confidence_set)
  tRNA_granges <- tRNA_granges[tRNA_granges$tRNA_anticodon != "NNN"] # Remove tRNAs which have unknown/undefined anticodon
  tRNA_granges <- GenomicRanges::trim(tRNA_granges + flanking_region)
  tRNA_granges$tRNA_instance <- stats::ave(seq_along(tRNA_granges), tRNA_granges$tRNA_type, tRNA_granges$tRNA_anticodon, FUN = seq_along) # Obtain the last index of the standard tRNA gene name

  if (!is.null(tRNA_name_map)){
    tRNA_gene_id <- paste0(GenomicRanges::seqnames(tRNA_granges), ".trna", tRNA_granges$no)
    map_vec <- stats::setNames(tRNA_name_map$GtRNAdb_id, tRNA_name_map$tRNAscan.SE_id)
    gene_names <- map_vec[tRNA_gene_id]
    tRNA_granges$gene_name
  } else {
    message("No 'tRNA_name_map' provided, the tRNA coordinates will not be translated.")
  }

  message("1 . COMPLETED\n", "2 . Finding overlaps and aggregating counts by tRNA gene..")
  if (inherits(data, "SummarizedExperiment")){
    counts_peak_matrix <- SummarizedExperiment::assay(data, "counts")
    peak_ranges <- SummarizedExperiment::rowRanges(data)

    overlaps <- GenomicRanges::findOverlaps(tRNA_granges, peak_ranges)
    M <- Matrix::sparseMatrix(i = S4Vectors::queryHits(overlaps), j = S4Vectors::subjectHits(overlaps), x = 1, dims = c(length(tRNA_granges), length(peak_ranges)))
    rownames(tRNA_counts) <- tRNA_granges$gene_name
    colnames(tRNA_counts) <- colnames(counts_peak_matrix)

    tRNA_matrix <- M %*% counts_peak_matrix

  } else if (inherits(data, c("Seurat", "ChromatinAssay"))){
    chrom_assay <- if(inherits(data, "Seurat")) data[[assay %||% "peaks"]] else data
    if (is.null(chrom_assay)) stop(paste("Requested assay not found in object."))

    frags <- Signac::Fragments(chrom_assay)
    if (length(frags) == 0) stop("No fragments found in the input assay.")

    names(tRNA_granges) <- paste0("tRNA_", seq_along(tRNA_granges))
    tRNA_matrix <- Signac::FeatureMatrix(fragments = frags, features = tRNA_granges, cells = colnames(chrom_assay), sep = name_sep, verbose = verbose)

  } else {
    stop("Unsupported data type. Must be SummarizedExperiment, Seurat, or ChromatinAssay.")
  }

  rownames(tRNA_matrix) <- tRNA_granges$gene_name

  message("2 . COMPLETED\n")
  if (isTRUE(save)){
    message("3 . Exporting tRNA expression matrix.")
    output_file <- GetOutputName(action = "file", out_name = out_name, out_directory = out_directory, save_format = "rds")
    saveRDS(tRNA_matrix, output_file)
    message("3 . COMPLETED\n")
  }
  return(tRNA_matrix) # Returns a tRNA gene expression matrix ready to be used as input in tTEscanR main anticodon functions
}

#' Filter Out Conditions With Low tRNA Cuts
#' @description
#' This function filters a tRNA expression matrix by removing conditions (columns) that fall below a specific total read count (\code{cutoff}).
#' It is useful for eliminating low-quality or poorly sequenced conditions that may bias downstream analyses.
#'
#' @param data A \code{matrix} or \code{data.frame} of tRNA gene expression data, with tRNA genes as rows and conditions as columns.
#' @param cutoff Numeric; minimum total number of tRNA cuts required to retain a condition in \code{data}. Defaults to 5000.
#'
#' @return A filtered \code{matrix} or \code{data.frame} with tRNAs below the cutoff removed.
#' @export
#'
#' @examples
#' data(default_tTEscanR_tRNA_data)
#' tRNA_data_filtered <- Filter_tRNACuts(data = default_tTEscanR_tRNA_data, cutoff = 5000)

Filter_tRNACuts <- function(data, cutoff = 5000){

  ###
  # CALL: User
  # DESCRIPTION: This function takes a tRNA gene expression count matrix and filters the tRNA genes
  # based on the cutoff parameter to avoid low-quality data in the downstream analyses.
  ###

  message("--- Filtering the tRNA gene abundance count matrix ---", paste("\n- Initial number of samples:", ncol(data)))
  keep_samples <- Matrix::colSums(data) >= cutoff
  if (!any(keep_samples)) stop("No samples pass the cutoff (", cutoff, "). Adjust the parameter.")
  data <- data[, keep_samples, drop = FALSE]
  message("- Samples retained: ", ncol(data))
  message("--- The tRNA gene abundance count matrix has been successfully filtered ---")

  return(data) # Returns the filtered tRNA gene expression count matrix.
}

#' Annotate the tRNA genes from tRNA tags
#'
#' @param data A \code{matrix} with tRNA tags as rows.
#' @param tRNA_bed Path to the directory that contains the .bed file
#' @param flanking_region Numeric; number of bases to include expand the region interrogated. Defaults to 100.
#'
#' @returns A \code{data} with the translated tRNA gene names.
#' @export
#'

Set_tRNAgenes <- function(data, tRNA_bed, flanking_region = 100){

  tRNA_table <- utils::read.delim(tRNA_bed, header = FALSE, stringsAsFactors = FALSE)

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
#' @return A single sparse matrix with as a combination of all the input matrices.
#' @export
#'
#' @examples
#' df1 <- matrix(c(1, 0, 0, 2), nrow = 2, dimnames = list(c("geneA", "geneB"), c("s1", "s2")))
#' df2 <- matrix(c(3, 0, 0, 4), nrow = 2, dimnames = list(c("geneB", "geneC"), c("s2", "s3")))
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

  # Create sparse matrix - memory-efficient, only the non-zero values are stored
  result_df <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(length(rnnew), length(cnnew)), dimnames = list(rnnew, cnnew))
  return(result_df)
}

#' Aggregates data by group
#' @description
#' This function calculates the row sums of a given matrix to combine columns that share the same group.
#'
#' @param data A \code{matrix} with the features to group for as columns.
#' @param group_labels A \code{vector} with the metadata features to group the columns in \code{data}
#'
#' @return A \code{matrix} with the conditions merged based on the metadata.
#' @export
#'
#' @examples
#' data <- data.frame(sample_1 = c(10, 5, 20), sample_2 = c(15, 8, 25),
#'                    sample_3 = c(12, 6, 22), sample_4 = c(1, 2, 3),
#'                    sample_5 = c(4, 5, 6), sample_6 = c(7, 8, 9))
#' rownames(data) <- c("gene_1", "gene_2", "gene_3")
#' groups <- c("cond_A", "cond_A", "cond_A", "cond_B", "cond_B", "cond_B")
#' data_combined <- GroupConditions(data = data, group_labels = groups)

GroupConditions <- function(data, group_labels) {

  # Check the dimensions
  if (ncol(data) != length(group_labels)) stop("Dimension mismatch: The number of columns in 'data' (", ncol(data), ") must match the length of 'group_labels' (", length(group_labels), ").")

  # Check for missing (NA) labels
  if (any(is.na(group_labels))) stop("NAs found in 'group_labels'. Please remove or handle missing grouping information before running.")

  # Check the order of the labels
  groups <- if (is.factor(group_labels)) droplevels(group_labels) else factor(group_labels) # Ensure group_labels contains factors

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
#' This function converts a \code{matrix} or \code{data.frame} into a tidy, long-format \code{tibble}. Optionally normalizes the values.
#'
#' @param data A table to be converted. Supported formats: \code{matrix} or \code{data.frame}.
#' @param normalize Logical; if \code{TRUE}, values are converted to relative frequencies. Defaults to \code{FALSE}.
#' @param rownames_to_column A character string specifying the name of the new column that will hold the former row names in \code{data}.
#' @param names_to A character string specifying the name of the new column that will hold the former column names in \code{data}.
#' @param values_to A character string specifying the name of the new column that will hold the corresponding values from the pivoted columns in \code{data}.
#'
#' @return A tibble of the input \code{data}
#' @export
#'
#' @examples
#' data(default_tTEscanR_tRNA_data)
#' tRNA_long_format <- TransformFormat(data = default_tTEscanR_tRNA_data, normalize = FALSE,
#'                                     rownames_to_column = "tRNA_genes",
#'                                     names_to = "condition", values_to = "abundance")
#' tTEobj <- Create_tTEscanR_Object(counts = default_tTEscanR_tRNA_data, assay = "tRNA")
#' tTEobj <- ComputeAnticodonUsage(object = tTEobj)
#' anticodon_long_format <- TransformFormat(data = tTEobj@assays$AnticodonUsage,
#'                                          normalize = TRUE, rownames_to_column = "anticodons",
#'                                          names_to = "condition", values_to = "usage")

TransformFormat <- function(data, normalize, rownames_to_column, names_to, values_to){

  ###
  # CALL: User
  # DESCRIPTION: This function transforms the input data from a wide format to a long format. With normalize = TRUE, performs the normalization of the data to account for sequencing depth.
  ###

  if (isTRUE(normalize)) data <- t(t(data) / colSums(data)) # Vectorized normalization (if required)
  if (!is.matrix(data)) data <- as.matrix(data) # Ensure data is a matrix

  num_row <- nrow(data)
  num_col <- ncol(data)

  long_format <- data.frame(rep.int(rownames(data), num_col), rep.int(colnames(data), rep.int(num_row, num_col)), as.vector(data), stringsAsFactors = FALSE) # Transform the data into a long format
  colnames(long_format) <- c(rownames_to_column, names_to, values_to)
  CheckDataFrame(long_format) # Evaluate that the data is properly defined

  return(long_format) # Returns the processed data (long-format and normalized if applicable).
}
