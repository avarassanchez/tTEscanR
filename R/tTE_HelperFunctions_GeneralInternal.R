###
# General functions used internally by the tTEscanR main functions
# 2 parts: (i) checking functions, (ii) input data control
###

# (i) Checking functions

CheckCodonFreqTable <- function(data, species, verbose = TRUE){

  ###
  # CALL: ConsistencyWithCodonFreq()
  # DESCRIPTION: This function checks the codon frequency per gene table needed to compute the codon usage.
  # If no user-defined codon frequency table is given, it loads (if enabled) the codon frequency matrix from the tTEscanR memory.
  ###

  if (is.null(data)){ # No user-defined codon frequency table has been given
    if (is.null(species)) stop("No 'codon_freq' provided and no 'species' specified.\n", "Specify 'hg38' for human or 'mm39' for mouse default data.")
    return(SelectDefaultData(species = species)) # Extract from the tTEscanR memory the default codon frequency table
  }

  # A user-defined codon-frequency table has been given
  CheckDataFrame(data = data, required_names = TRUE) # Check that the format of the table is suitable
  actual_codons <- rownames(data)

  return(data) # Returns the codon frequency per gene table
}

CheckDataFrame <- function(data, required_names = TRUE){

  ###
  # CALL: CheckCodonFreqTable()
  # DESCRIPTION: This function evaluates the format of the input data. There is no return value, the execution will stop upon an error.
  ###

  if (is.null(data)) stop("The dataset is empty or could not be found.") # Control - in case the data object is empty
  dims <- dim(data)

  if (is.null(dims) || dims[1] < 2 || dims[2] < 2){ # We use at least 2 columns and 2 rows to define a table - evaluate that it is not empty
    row_val <- ifelse(is.null(dims), 0, dims[1])
    col_val <- ifelse(is.null(dims), 0, dims[2])
    stop(sprintf("Incorrect dimensions. Expected at least 2x2. Found: %d x %d", row_val, col_val))
  }

  if (required_names){ # The required_names parameter is used to specify if the data needs to contain row and column names
    row_names <- rownames(data)
    col_names <- colnames(data)
    if (is.null(row_names) || is.null(col_names)) stop("Dataset labels missing. Both rownames and colnames are required.")
    if (any(row_names == "", na.rm = TRUE) || any(col_names == "", na.rm = TRUE)) stop("Dataset labels are incomplete. Empty strings detected in rownames or colnames.")
  }
  return(invisible(NULL))
}

CheckGeneAnnotation <- function(vector1, vector2, verbose = TRUE){

  ###
  # CALL: ConsistencyWithCodonFreq()
  # DESCRIPTION: This function compares the gene annotation of two vectors to ensure consistency.
  ###

  # Check in which format each vector of genes is
  format1 <- IsEnsemblID(vector1)
  format2 <- IsEnsemblID(vector2)

  if (format1 != format2) stop("Annotation mismatch detected!\n", "- Codon frequency table: ", format1, "\n", "- mRNA data: ", format2, ".\n", "Please ensure both datasets use the same gene naming convention.")
  if (verbose) message("Both vectors are in the same format: ", format1, ".")
}

# (ii) Input data control

SelectDefaultData <- function(species, action = "codon_freq"){

  ###
  # CALL: Multiple internal functions
  # DESCRIPTION: This function enables to retrieve from the tTEscanR memory the required codon frequency-per-gene table.
  ###

  if (action == "codon_freq"){
    valid_species <- c("hg38", "mm39")
    if (!(species %in% valid_species)) stop("Incorrect 'species'.\n", "Supported formats: hg38 for human and mm39 for mouse.")

    message(paste("- The default dataset:", species, "will be used as 'codon_freq'."))

    # Construct the name of the object to retrieve
    data_name <- paste("codon_freq_table_canonical", species, sep = "_")
    codon_freq <- get(data_name, envir = asNamespace("tTEscanR"))

    return(codon_freq) # Return the retrieved codon frequency-per-gene table
  }

}

IdentifyInputFormat <- function(data, mode = c("fix", "flexible")) {

  ###
  # CALL: Multiple internal functions
  # DESCRIPTION: Check the format of the input data and return a string with the format of the input data.
  ###

  mode <- match.arg(mode) # Validate the string arguments

  is_table_like <- function(x) {
    is.data.frame(x) || is.matrix(x) || inherits(x, "dgCMatrix")
  }

  if (is.list(data) && !is.data.frame(data)) {
    if (mode == "fix") {
      valid_types <- vapply(data, is_table_like, logical(1))
      if (!all(valid_types)) stop("Incorrect 'data' format.\nSupported formats: list of dataframes or matrices.")
    }
    return("list")
  }

  if (is_table_like(data)) return("single")
  if (is.character(data)) return("single")
  if (is.numeric(data) && mode == "flexible") return("single")

  if (mode == "fix") {
    stop("Incorrect 'data' format. \nSupported formats: data.frame, matrix, sparse matrix, or a list of these.")
  } else {
    stop("Incorrect 'data' format. \nSupported formats: table-like objects, numeric/chracter vectors, or lists.")
  }
}

FilterByMetadata <- function(data, metadata, id_col = NULL, verbose = TRUE) {

  conditions <- colnames(data)

  if (is.null(id_col)) {
    matching_counts <- vapply(metadata, function(col) sum(conditions %in% col), numeric(1))
    matching_col <- names(which.max(matching_counts))
    if (matching_counts[matching_col] == 0) stop("No matching Sample IDs found between 'data' and 'metadata.'")
    message(sprintf("- Auto-detected metadata ID column: '%s'", matching_col))
  } else {
    matching_col <- id_col
  }

  target_col <- metadata[[matching_col]]

  # data_dups <- conditions[duplicated(conditions) | duplicated(conditions, fromLast = TRUE)]
  # meta_dups <- target_col[duplicated(target_col) | duplicated(target_col, fromLast = TRUE)]

  data_dups <- unique(conditions[duplicated(conditions)])
  meta_dups <- unique(target_col[duplicated(target_col)])
  all_forbidden <- unique(c(data_dups, meta_dups))
  if (verbose && length(all_forbidden) > 0) message("- Found ", length(all_forbidden), " duplicated IDs. Removing all instances from analysis.")

  shared_unique <- setdiff(intersect(conditions, target_col), all_forbidden)
  if (length(shared_unique) == 0) stop("No unique, non-duplicated matching IDs found between 'data' and 'metadata.'")

  filtered_data <- data[, shared_unique, drop = FALSE]
  # temp_metadata <- metadata[metadata[[matching_col]] %in% shared_unique, , drop = FALSE]

  reorder_idx <- match(shared_unique, target_col)
  filtered_metadata <- metadata[reorder_idx, , drop = FALSE]

  rownames(filtered_metadata) <- shared_unique

  if (verbose) message("- Final dataset contains ", ncol(filtered_data), " unique samples.")

  return(list(data = filtered_data, metadata = filtered_metadata))
}
