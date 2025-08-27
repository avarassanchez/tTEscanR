SelectDefaultData <- function(species, filter = NULL){

  ###
  # CALL: Multiple internal functions
  # DESCRIPTION: This function enables to retrieve from the tTEscanR memory the required codon frequency-per-gene table.
  ###

  # Control steps
  if (!(species %in% c("hg38", "mm39"))) stop("Incorrect `species` input name.\n", "Supported formats: hg38 for human and mm39 for mouse.")
  if (!(filter %in% c("canonical", "length"))) stop("Incorrect `filter` input name.\n", "Supported formats: canonical, length.")

  message(paste("- The default dataset:", species, "will be used as `codon_freq`."))

  # HUMAN DEFAULT DATA
  if (species == "hg38" && filter == "canonical") codon_freq <- codon_freq_table_hg38_canonical
  if (species == "hg38" && filter == "length") codon_freq <- codon_freq_table_hg38_length

  # MOUSE DEFAULT DATA
  if (species == "mm39" && filter == "canonical") codon_freq <- codon_freq_table_mm39_canonical
  if (species == "mm39" && filter == "length") codon_freq <- codon_freq_table_mm39_length

  return(codon_freq) # Return the retrieved codon frequency-per-gene table
}

IdentifyInputFormat <- function(data,  mode = "fix"){

  ###
  # CALL: Multiple internal functions
  # DESCRIPTION: This function checks the format of the input data.
  # Returns a character string with the actual format of the input data.
  ###

  if(is.character(data)){ # DEALING WITH VECTORS
    format <- "vector"
  } else if (is.data.frame(data) || is.matrix(data) || inherits(data, "dgCMatrix")){ # DEALING WITH DATAFRAME OR MATRICES
    format <- "table"
    CheckDataFrame(data)
  } else if (inherits(data, "list")){ # DEALING WITH LISTS

    if (mode == "fix"){ # Check the actual content of the list
      for (i in 1:length(data)){
        if (!(is.data.frame(data[[i]]) || is.matrix(data[[i]]) || inherits(data[[i]], "dgCMatrix"))) stop("Incorrect data format.\n", "Supported formats: vector, dataframe or matrix.")
      }
    }

    format <- "list"

  } else { # WRONG INPUT FORMAT
    stop("Incorrect data format.\n", "Supported formats: vector, dataframe or matrix.")
  }
  return(format)
}

FilterByMetadata <- function(data, metadata, verbose = TRUE){

  ###
  # CALL: ExecuteDESeq2runner() and Compute_tTE()
  # DESCRIPTION: This function evaluates the entries in `data` and `metadata`.
  # If there is not a perfect match it filters them accordingly.
  # Returns the filtered matrix and metadata in a list.
  ###

  conditions <- colnames(data) # Extract the column names

  # Assess how many column names are present in each metadata column
  matching_counts <- sapply(metadata, function(col) sum(conditions %in% col))

  # Retrieves the name of the metadata with more matches
  matching_col <- names(which.max(matching_counts)) # Upon a tie retains the 1st appearing column

  # Retrieve the conditions that are present in `data` and `metadata`
  shared_conditions <- intersect(conditions, metadata[[matching_col]])

  # Filter the `data` and `metadata` accordingly
  filtered_metadata <- metadata[which(metadata[[matching_col]] %in% shared_conditions), ]
  filtered_data <- data[, shared_conditions]

  if (!identical(dim(data), dim(filtered_data)) && verbose) message("- Filtering steps applied to `data` to match the features in `metadata`.")
  if (!identical(dim(metadata), dim(filtered_metadata)) && verbose) message("- Filtering steps applied to `metadata` to match the features in `data`.")

  # Check that the data and metadata have a suitable format
  CheckDataFrame(filtered_data)
  CheckDataFrame(filtered_metadata)
  return(list(filtered_data, filtered_metadata))
}
