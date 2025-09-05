CheckCodonFreqTable <- function(data, species, filter = "canonical", verbose = TRUE){

  ###
  # CALL: ConsistencyWithCodonFreq()
  # DESCRIPTION: This function checks the codon frequency per gene table needed to compute the codon usage.
  # If no user-defined codon frequency table is given, it loads (if enabled) the codon frequency matrix from the tTEscanR memory.
  ###

  if (is.null(data)){ # No user-defined codon frequency table has been given

    if (is.null(species)) stop("No `codon_freq` has been given as input and no `species` is specified to use the default codon frequency per gene table.\n",
                               "Please, specify `hg38` for human or `mm39` for mouse default data.")
    if (verbose) message("- No `codon_freq` has been given as input.")

    # There are 2 filter options: canonical and length
    if(!(filter %in% c("canonical", "length"))) stop("An invalid `filter` has been given as input.\n", "Please, specify `canonical` or `length`.")

    # Extract from the tTEscanR memory the default codon frequency table
    data <- SelectDefaultData(species = species, filter = filter)

  } else{ # A user-defined codon-frequency table has been given
    CheckDataFrame(data = data, names = TRUE) # Check that the format of the table is suitable
  }

  # Ensure that all the codons in the codon frequency table are sense codons
  if (verbose) message("- Filtering the codon frequency per gene table by the sense codons.")
  data <- data[rownames(data) %in% sense_codons, ]

  CheckDataFrame(data = data, names = TRUE) # Check the final codon frequency table, ensure that it has proper dimensions
  return(data) # Returns the codon frequency per gene table filtered by sense codons
}

CheckDataFrame <- function(data, names = TRUE){

  ###
  # CALL: CheckCodonFreqTable()
  # DESCRIPTION: This function evaluates the format of the input data. There is no return value, the execution will stop upon an error.
  ###

  # Control - in case the data object is empty
  if (is.null(data)) stop("The dataset is empty or could not be found.")

  # We use at least 2 columns and 2 rows to define a table - evaluate that it is not empty
  if (ncol(data) < 2 || nrow(data) < 2) stop("Incorrect dimensions of the dataset.\n", "- Number of columns: ", ncol(data), ".\n", "- Number of rows: ", nrow(data), ".\n")

  # The names parameter is used to specify if the data needs to contain row and column names
  if (names == TRUE && (is.null(rownames(data)) || is.null(colnames(data)))) stop("The labels of the dataset could not be found.\n",
                                                                                  "Check that rownames() and columnames() do not give empty outcomes.")
}

CheckGeneAnnotation <- function(vector1, vector2, verbose = TRUE){

  ###
  # CALL: ConsistencyWithCodonFreq()
  # DESCRIPTION: This function compares the gene annotation of two vectors and, if required to translate one of them to ensure consistency.
  ###

  # Check in which format each vector of genes is
  format1 <- IsEnsemblID(vector1, action = "check")
  format2 <- IsEnsemblID(vector2, action = "check")

  if (format1 != format2) stop("The vectors are not in the same format: vector1 in ", format1, " and vector2 in ", format2, ".")
  if (verbose) message("Both vectors are in the same format: ", format1, ".")
}
