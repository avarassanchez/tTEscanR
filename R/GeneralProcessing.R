
#' Selection of the Optimal tRNA Cut Cutoff
#'
#' @param data tRNA gene expression count \code{matrix} with tRNA genes as rows and conditions as columns.
#' @param cutoffs Set of values to test to search for the optimal tRNA cuts threshold. Defaults to seq(50, 10000, by = 50).
#' @param generate_plot Logic; if \code{TRUE}, generates a correlation plot. Defaults to \code{TRUE}.
#' @param slope_threshold Numeric; value to consider for the determination of the correlation stability.
#'
#' @returns Table with the optimal cutoff at the anticodon isoacceptor and amino acid isotype.
#' @export
#'
#' @examples
#' data(subset_tRNA_data)
#' optimal_tRNA_cutoffs <- tRNACutsCutoff(data = subset_tRNA_data,
#'                                        cutoffs = seq(50, 100, by = 10),
#'                                        generate_plot = FALSE)

tRNACutsCutoff <- function(data, cutoffs = seq(50, 10000, by = 50), generate_plot = TRUE, slope_threshold = 0.001) {


  message("1 . Computing reference tTEscanR object.")
  object_ref <- suppressMessages(Create_tTEscanR_Object(counts = data, assay = "tRNA", verbose = FALSE))
  object_ref <- suppressMessages(ComputeAnticodonUsage(object_ref, verbose = FALSE))
  object_ref <- suppressMessages(ComputeAAUsage(object_ref, level = "supply", verbose = FALSE))

  ref_anticodon <- rowSums(object_ref@assays$AnticodonUsage)
  ref_supply <- rowSums(object_ref@assays$AASupply)
  ref_anticodon <- ref_anticodon / sum(ref_anticodon)
  ref_supply <- ref_supply / sum(ref_supply)
  message("  1 . COMPLETED\n", "2 . Transform data.")

  data_long <- TransformCounts(data)
  message("  2 . COMPLETED\n", "3 . Compute correlations.")

  cor_results <- ComputeCorrelations(data_long, ref_anticodon, ref_supply, cutoffs)
  message("  3 . COMPLETED\n", "4 . Determine optimal cutoff.")

  optimal_cutoff <- Selection_Cutoff(cor_results, slope_threshold = slope_threshold)
  message("  4 . COMPLETED")

  plot <- NULL
  if (isTRUE(generate_plot)){
    message("5 . Generate plot.")
    plot <- CorrelationCutoffPlot(cor_results)
    message("  5 . COMPLETED")
  }

  if (!is.null(plot)) return(list(optimal_cutoff = optimal_cutoff, plot = plot)) else return(optimal_cutoff)
}

#' Convert a Chromatin Assay into a tRNA Expression Matrix
#' @description
#' This function extracts tRNA expression data from a given \code{ChromatinAssay} or \code{SeuratObject} and formats it into a standardized expression matrix.
#' if no \code{tRNA_annotations} is provided, internal annotations for human (\code{"hg38"}) and mouse (\code{"mm39"}) can be used.
#' Optionally, the resulting matrix can be saved to a file for downstream analysis.
#'
#' @param chrom A \code{ChromatinAssay} or \code{SeuratObject} containing tRNA data.
#' @param assay Optional; a character string specifying the name of the assay to retrieve from the \code{chrom} if it is a \code{SeuratObject}.
#' @param tRNA_annotations A directory path to tRNA annotations in .ss format.
#' @param flanking_region Integer; number of nucleotides that form the flanking region of each tRNA. Defaults to 100.
#' @param species Either \code{"hg38"} (human) or \code{"mm39"} (mouse) to specify which default reference tRNA prediction file to use. Required if \code{tRNA_annotations} is not provided.
#' @param name_sep A string delimiter to format the tRNA gene names in the output matrix. Defaults to \code{c("-", "-")}.
#' @param save Logical; if \code{TRUE}, the output tRNA expression matrix will be saved as a \code{.rds} file. Defaults to \code{TRUE}.
#' @param out_name Optional; name for the output file. Required if \code{save} is enabled.
#' @param out_directory Optional; path to the directory where the output file will be saved. Required if \code{save} is enabled.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return A tRNA expression matrix with tRNA genes as rows and samples or conditions as columns.
#' @export

Get_tRNAMatrix <- function(chrom, assay = NULL, tRNA_annotations, flanking_region = 100, species = NULL, name_sep = c("-", "-"),
                           save = TRUE, out_name = NULL, out_directory = NULL, verbose = TRUE){

  ###
  # CALL: User
  # DESCRIPTION: This function takes a Seurat or Signac object from scATAC-seq data and the corresponding tRNA
  # reference gene annotations to generate a tRNA gene expression count matrix using the proper Signac built-in functions.
  ###

  message("1 . Checking the format of the input data.")
  if (verbose) message("- Checking the `chrom_obj.`")
  if (!(inherits(chrom, "Seurat")) && !(inherits(chrom_obj, "ChromatinAssay")))  stop("The `chrom` input needs to be a Seurat or a ChromatinAssay object.")

  if (inherits(chrom, "Seurat")) { # DEALING WITH SEURAT OBJECT
    if (verbose) message("- The input `chrom_obj` is a Seurat object.")
    if (is.null(assay)) assay <- "peaks"
    if (!(assay %in% names(chrom@assays))) stop(paste("The `chrom` Seurat object does not contain a", assay, "assay."))

    chrom_obj <- chrom
  }

  if (inherits(chrom_obj, "ChromatinAssay")) { # DEALING WITH CHROMATIN ASSAY
    if (verbose) message("- The input `chrom_obj` is a ChromatinAssay object.\n", "- Creating a Seurat object from the input ChromatinAssay.")
    chrom_obj <- Seurat::CreateSeuratObject(counts = chrom, assay = assay,)
  }

  if (verbose) message ("- Checking the `tRNA_annotations.`")
  tRNA_granges <- tRNAscanImport::import.tRNAscanAsGRanges(tRNA_annotations)

  message("  1 . COMPLETED\n", "2 . Extracting metadata from the tRNA gene annotation file.")
  tRNA_granges <- tRNA_granges[which(tRNA_granges$tRNA_anticodon != "NNN")] # remove tRNAs which have unknown/undefined anticodon
  GenomicRanges::ranges(tRNA_granges) <- GenomicRanges::ranges(tRNA_granges) + 100 # include the flanking 100 nt of each tRNA

  message("  2 . COMPLETED\n", "3 . Creating tRNA expression matrix.")
  names(tRNA_granges) <- 1:length(tRNA_granges)
  tRNA_granges$gene_name <- paste(names(tRNA_granges), tRNA_granges$tRNA_type, tRNA_granges$tRNA_anticodon, sep = '_')
  tRNA_granges$gene_biotype <- 'tRNA'

  tRNA_expression_matrix <- Signac::FeatureMatrix(fragments = Signac::Fragments(chrom_obj[[assay]]), features = tRNA_granges, cells = colnames(chrom_obj), sep = name_sep, verbose = verbose)

  message("  3 . COMPLETED\n", "4 . Exporting tRNA expression matrix.")
  if (isTRUE(save)){
    if (verbose) message("- Saving the tRNA expression into a file.")
    output_file <- GetOutputName(action = "file", out_name = out_name, out_directory = out_directory, save_format = "rds")
    saveRDS(tRNA_expression_matrix, output_file)
  }
  message("  4 . COMPLETED\n")
  return(tRNA_expression_matrix) # Returns a tRNA gene expression matrix ready to be used as input in tTEscanR main anticodon functions
}

#' Translate chromatin coordinates to tRNA genes
#'
#' @param data A tRNA expression matrix with tRNA genes as rows and samples or conditions as columns. Can be generated using \code{\link{Get_tRNAMatrix}}
#' @param tRNA_bed A string indicating the path to the \code{".bed"} file from the organism's reference genome
#' @param flanking_region Integer; number of nucleotides that form the flanking region of each tRNA. Defaults to 100.
#'
#' @return The \code{data} tRNA expression data with the tRNA genes annotated for the given fragment coordinates.
#' @export
#'

Set_tRNAgenes <- function(data, tRNA_bed, flanking_region = 100){

  tRNA_table <- utils::read.delim(tRNA_bed, header = FALSE)
  tRNA_table <- tRNA_table[ , 1:4]
  tRNA_table$V2 <- tRNA_table$V2 + 1
  tRNA_table$V2 <- tRNA_table$V2 - flanking_region
  tRNA_table$V3 <- tRNA_table$V3 + flanking_region
  tRNA_table$region <- paste(tRNA_table$V1, tRNA_table$V2, tRNA_table$V3, sep = "-")

  translate_table <- tRNA_table[ , 4:5]

  names_row <- rownames(data)
  ref_names <- translate_table[ ,2]
  new_names <- translate_table[ ,1]

  tr <- stats::setNames(new_names, nm = ref_names)
  tRNA_id <- dplyr::recode(names_row, !!!tr)

  rownames(data) <- tRNA_id

  return(data)
}

#' Filter Out Conditions With Low tRNA Cuts
#' @description
#' This function filters a tRNA expression matrix by removing tRNA genes (rows) that fall below a specific total read count (\code{cutoff}).
#' It is useful for eliminating low-quality or poorly sequenced conditions that may bias downstream analyses.
#'
#' @param data A \code{matrix} or \code{data.frame} of tRNA gene expression data, with tRNA genes as rows and conditions as columns.
#' @param cutoff Numeric; minimum total number of tRNA cuts required. Defaults to 5000.
#'
#' @return A filtered \code{matrix} or \code{data.frame} with tRNAs below the cutoff removed.
#' @export
#'
#' @examples
#' data(subset_tRNA_data)
#' tRNA_data_filtered <- tRNACutsFilter(data = subset_tRNA_data, cutoff = 5000)

tRNACutsFilter <- function(data, cutoff = 5000){

  ###
  # CALL: User
  # DESCRIPTION: This function takes a tRNA gene expression count matrix and filters the tRNA genes
  # based on the cutoff parameter to avoid low-quality data in the downstream analyses.
  ###

  column_name_counts <- Matrix::colSums(data)
  unique_column_names <- names(column_name_counts)[column_name_counts < cutoff]
  indices <- which(colnames(data) %in% unique_column_names)
  if(!is.null(indices)) data <- data[, -indices] else stop("No features pass the cutoff. Plase, adjust it accordingly.")

  return(data) # Returns the filtered tRNA gene expression count matrix.
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
#' data(subset_tRNA_data)
#' tRNA_long_format <- DataToLongFormat(data = subset_tRNA_data, normalize = FALSE,
#'                                      rownames_to_column = "tRNA_genes",
#'                                      names_to = "condition", values_to = "expression")
#' tTEobj <- Create_tTEscanR_Object(counts = subset_tRNA_data, assay = "tRNA")
#' tTEobj <- ComputeAnticodonUsage(object = tTEobj)
#' anticodon_long_format <- DataToLongFormat(data = tTEobj@assays$AnticodonUsage,
#'                                          normalize = TRUE, rownames_to_column = "anticodons",
#'                                          names_to = "condition", values_to = "usage")

DataToLongFormat <- function(data, normalize, rownames_to_column, names_to, values_to){

  ###
  # CALL: User and ComputeMeanUsage()
  # DESCRIPTION: This function transforms the input data from a wide format to a long format.
  # Moreover, with normalize = TRUE, performs the normalization of the data to account for sequencing depth.
  ###

  data_processed <- if (isTRUE(normalize)) sweep(data, 2, colSums(data), FUN = "/") else data # Normalize the data if required

  # Transform the data into a long format
  long_format <- data_processed %>% as.data.frame() %>%
    tibble::rownames_to_column(var = rownames_to_column) %>%
    tidyr::pivot_longer(-{{ rownames_to_column }}, names_to = names_to, values_to = values_to)

  CheckDataFrame(long_format) # Evaluate that the data is properly defined
  return(long_format) # Returns the processed data (long-format and normalized if applicable).
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
    df <- as.matrix(df)          # force to matrix
    storage.mode(df) <- "numeric"  # force numeric

    cnold <- colnames(df)
    rnold <- rownames(df)

    # Collect the names of the rows and columns
    cnnew <- union(cnnew, cnold)
    rnnew <- union(rnnew, rnold)

    cindnew <- match(cnold, cnnew)
    rindnew <- match(rnold, rnnew)

    # Finds non-zero values
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
#' @param group.labels A \code{vector} with the metadata features to group the columns in \code{data}
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
#' data_combined <- GroupConditions(data = data, group.labels = groups)

GroupConditions <- function(data, group.labels) {
  # Ensure group.labels is a factor
  groups <- factor(group.labels)
  n_groups <- length(levels(groups))

  # If sparse matrix (from Matrix package)
  if (inherits(data, "dgCMatrix") || inherits(data, "dgTMatrix") || inherits(data, "Matrix")) {
    # Create sparse indicator matrix for groups
    G <- Matrix::sparse.model.matrix(~0 + groups)

    # Sparse matrix multiplication (fast + memory-efficient)
    res <- data %*% G
    res <- as.data.frame(as.matrix(res))

  } else {
    # Dense matrix: use base::rowsum (optimized C code)
    res <- t(rowsum(t(as.matrix(data)), group = groups))
    res <- as.data.frame(res)
  }

  # Add names
  rownames(res) <- rownames(data)
  colnames(res) <- levels(groups)

  return(res)
}

