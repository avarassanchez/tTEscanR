GroupAA <- function(data, amino_acids_function){

  ###
  # CALL: ComputeAAUsage()
  # DESCRIPTION: This function takes a codon usage or anticodon usage and groups the features into amino acids.
  ###

  # Retrieve the AA from the codons or anticodons
  AAextracted <- as.vector(Biostrings::translate(amino_acids_function, no.init.codon = TRUE))
  if (is.null(AAextracted)) stop("No amino acids could be retrieved.")
  AAunique <- sort(unique(AAextracted))

  # Generate empty matrix (all values are 0) based on the number of codon/anticodon (rows) and the conditions in data (columns)
  AAmatrix <- matrix(data = 0, nrow = length(AAunique), ncol = length(colnames(data)))
  dimnames(AAmatrix) <- list(AAunique, colnames(data))

  for (i in 1:length(AAunique)) { # Fill the matrix by summing the counts of the codons/anticodons
    indices <- which(AAextracted == AAunique[i])
    AAmatrix[i, ] <- Matrix::colSums(data[indices, , drop = FALSE])
  }

  CheckDataFrame(data = AAmatrix) # Ensure that the generated matrix is valid
  return(AAmatrix) # Return amino acid by condition matrix
}

RetrieveAAUsageData <- function(object, data_section, assay_id, data_function){

  ###
  # CALL: ComputeAAUsage()
  # DESCRIPTION: This function runs the intermediate checks for each of the sections considered in ComputeAAUsage().
  ###

  IsIn_tTEscanR_Object(object = object, slot = "assays", section = data_section, verbose = FALSE) # Check if the required assays are present in the tTEscanR object
  CheckDataFrame(data = object@assays[[data_section]]) # Check that the data is in a suitable format

  data <- list(object@assays[[data_section]]) # Store the data and the assay_id into the previously defined lists

  if(data_function == "sense"){ # Extract AA from the codons - AADemand
    amino_acids_function <- list(Biostrings::DNAStringSet(sense_codons))
  } else { # Extract AA from the anticodons - AASupply
    amino_acids_function <- list(Biostrings::reverseComplement(Biostrings::DNAStringSet(rownames(object@assays$AnticodonUsage))))
  }

  return(list(data, assay_id, amino_acids_function)) # The output consists of list with the data, the id, and the function corresponding to the input section
}
