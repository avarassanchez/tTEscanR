ExtractGenes <- function(ensembl, filter, retain.mitochondrial, verbose = TRUE){

  ###
  # CALL: CallingEnsembl()
  # DESCRIPTION: This function uses biomaRt to get the genes of a target organism (based on the input ensembl) and returns a table.
  ###

  if (verbose) message("- Retriving protein coding genes & filtering out mitochondrial genes.")
  if (isFALSE(retain.mitochondrial)){ # Exclude mitochondrial genes
    prot_coding_ensembl <- biomaRt::getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "chromosome_name", "external_gene_name",
                                                         "transcript_start", "transcript_end", "transcript_is_canonical"),
                                          filters = c("biotype", "chromosome_name"), values = list("protein_coding", setdiff(c(1:22, "X", "Y"), "MT")), mart = ensembl)
  } else { # Retain all protein coding genes
    prot_coding_ensembl <- biomaRt::getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "chromosome_name", "external_gene_name",
                                                         "transcript_start", "transcript_end", "transcript_is_canonical"),
                                          filters = "biotype", values = "protein_coding", mart = ensembl)
  }
  filter_name <- if (filter == "length") "largest" else "canonical"
  if (verbose) message("- Selecting (if many) the ", filter_name, " transcript for each gene.")

  if (filter == "length"){ # If many transcripts are available uses the filter parameter to select a single case
    if (verbose) message("- Computing the length of each transcript.")
    prot_coding_ensembl$transcript_length <- as.numeric(abs((prot_coding_ensembl$transcript_start+1) - prot_coding_ensembl$transcript_end))
    prot_coding_ensembl <- prot_coding_ensembl[order(prot_coding_ensembl$ensembl_gene_id, -prot_coding_ensembl$transcript_length), ]
    prot_coding_ensembl <- prot_coding_ensembl[!duplicated(prot_coding_ensembl$ensembl_gene_id), ]
  } else {
    prot_coding_ensembl <- subset(prot_coding_ensembl, prot_coding_ensembl$transcript_is_canonical == 1)
  }
  return(prot_coding_ensembl) # Returns a table with the features retrieved by getBM()
}

ExtractSequences <- function(transcripts, ensembl){

  ###
  # CALL: CallingEnsembl()
  # DESCRIPTION: This function takes a list of transcript ids and retrieves the actual nucleotide sequence.
  ###

  transcript_sequences <- biomaRt::getSequence(id = transcripts, type = "ensembl_transcript_id", seqType = "coding", mart = ensembl)
  return(transcript_sequences) # Returns a table with 2 columns: (i) transcript id, and (ii) nucleotide sequence
}

CallingEnsembl <- function(dataset_name, transcripts, filter, retain.mitochondrial, verbose){

  ###
  # CALL: ObtainCodonFreqPerGene()
  # DESCRIPTION: This function retrieves nucleotide sequence of the protein-coding genes of a target organism, in an exploratory or targeted (user-defined list) way.
  ###

  results.list <- list() # Empty list to store the results

  message("2 . Retrieving Ensembl dataset.")
  ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = dataset_name)

  if (is.null(transcripts)){
    if (verbose) message("- Retrieving the protein-coding transcripts.")

    transcripts <- ExtractGenes(ensembl = ensembl, filter = filter, retain.mitochondrial = retain.mitochondrial) # Get table with the features of all protein-coding transcripts

    # Define a translator table with the different ids of a gene
    translator_table <- data.frame(transcripts$ensembl_transcript_id, transcripts$ensembl_gene_id, transcripts$external_gene_name)
    colnames(translator_table) <- c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name")
    results.list <- append(results.list, list(translator_table))

    transcripts <- transcripts$ensembl_transcript_id # Retrieve the transcript ids
  }

  message("  2 . COMPLETED\n", "3 . Extracting the genomic sequence of each transcript.")
  transcript_sequences <- ExtractSequences(transcripts = transcripts, ensembl = ensembl) # Based on a list of ids retrieve the actual sequence
  if (nrow(transcript_sequences) == 0) stop(paste("No transcripts were found.\n", "Please revise the `transcripts` if applicable."))

  unavailable <- which(transcript_sequences$coding == "Sequence unavailable")
  if (length(unavailable) != 0){
    if (verbose) message(paste("- There are", length(unavailable), "transcripts which sequence is unavailable.\n"), "- These transcripts will be removed.")
    if(length(unavailable) == length(transcript_sequences)) stop("No transcripts were found.\n", "Please revise the `transcripts` if applicable.")
    transcript_sequences <- transcript_sequences[-unavailable, ]
  }

  if (verbose) message(paste("- Number of protein coding transcripts:", nrow(transcript_sequences)))
  results.list <- append(results.list, list(transcript_sequences)) # Use the list to ensure that all the entries are interpreted as a single object
  message("  3 . COMPLETED")
  return(results.list) # Returns a list with the translator table between gene identifiers (if possible) and the sequences of the genes
}

CheckFASTAFormat <- function(file){

  ###
  # CALL: ObtainCodonFreqPerGene()
  # DESCRIPTION: This function checks if the input file corresponds to a FASTA file, otherwise reports an error.
  ###

  lines <- readLines(file)
  header <- grepl("^>", lines)
  if (!any(header)) stop("No headers found ('>'). Not a valid FASTA file.")
  if (length(lines) < 2) stop("File too short to contain valid FASTA records.")

  invalid_lines <- which(!header && !grepl("^[ATGCatgcnNUu]+$", lines))
  if (length(invalid_lines) > 0) stop("Invalid sequence characters found on lines:", paste(invalid_lines, collapse = ", "))
}

FromFASTAtoTable <- function(data, transcripts, retain.mitochondrial, verbose){

  ###
  # CALL: ObtainCodonFreqPerGene()
  # DESCRIPTION: This function transforms an input FASTA file into a table with 2 columns: (i) transcript id, and (ii) nucleotide sequence.
  ###

  message("2. Extracting the genomic sequence of each transcript.")
  if (verbose) message("- Retrieving the protein-coding transcripts.")
  protein_coding_data <- data[grep("transcript_biotype:protein_coding", names(data))] # Extract the protein-coding genes

  # Generate a table with the ids and the nucleotide sequences
  transcript_sequences <- data.frame(coding = unname(as.character(protein_coding_data)), ensembl_transcript_id = names(protein_coding_data), stringsAsFactors = FALSE)
  if (nrow(transcript_sequences) == 0) stop("There are no protein coding transcripts.")

  if (isFALSE(retain.mitochondrial)){ # Filter out mitochondrial genes
    if (verbose) message("- Mitochondrial genes will be removed.")
    mitochondrial_index <- grepl("chromosome:[^:]+:MT|gene_symbol:MT-|mitochondrial", transcript_sequences$ensembl_transcript_id, ignore.case = TRUE)
    if (nrow(transcript_sequences) == length(mitochondrial_index)) stop("All the protein coding transcripts correspond to mitochondrial genes.")
    transcript_sequences <- transcript_sequences[!mitochondrial_index, ]
  }

  # The entries of the FASTA file contain all the pieces of information - Here we retrieve each part independently
  ensembl_transcript_id <- stringr::str_extract(transcript_sequences$ensembl_transcript_id, "ENST[0-9]+\\.?[0-9]*")
  ensembl_gene_id <- stringr::str_extract(transcript_sequences$ensembl_transcript_id, "ENSG[0-9]+\\.?[0-9]*")
  external_gene_name <- stringr::str_extract(transcript_sequences$ensembl_transcript_id, "gene_symbol:[^ ]+")
  external_gene_name <- stringr::str_remove(external_gene_name, "gene_symbol:")

  transcript_sequences <- data.frame(coding = transcript_sequences$coding, ensembl_gene_id = ensembl_gene_id,
                                     ensembl_transcript_id = ensembl_transcript_id, external_gene_name = external_gene_name,
                                     transcript_length = nchar(transcript_sequences$coding), stringsAsFactors = FALSE)

  if(!is.null(transcripts)){ # Targeted approach based on the transcripts input parameter
    if (verbose) message("- Trimming by the ids included as `transcripts`.")

    targeted_indexes <- list(targeted_indexes_gene_id = which(transcript_sequences$ensembl_gene_id %in% transcripts),
                             targeted_indexes_transcript_id = which(transcript_sequences$ensembl_transcript_id %in% transcripts),
                             targeted_indexes_gene_name = which(transcript_sequences$external_gene_name %in% transcripts))

    if ((length(targeted_indexes$targeted_indexes_gene_id) == 0) && (length(targeted_indexes$targeted_indexes_transcript_id) == 0) && (length(targeted_indexes$targeted_indexes_gene_name) == 0)){
      stop("None of the ids given in `transcripts` have been found.\n", "Supported formats: Ensembl transcript id, Ensembl gene id and external gene name.")
    }

    # Select the column of the ids that have a higher match with the input list of ids
    targeted_indexes <- targeted_indexes[[which.max(lengths(targeted_indexes))]]
    transcript_sequences <- transcript_sequences[targeted_indexes, ]
  }

  if (nrow(transcript_sequences) == 0) stop("There are no protein coding transcripts.")
  if (verbose) message(paste("- Number of protein coding transcripts:", nrow(transcript_sequences)))

  # Generate a translator table with the above variables
  translator_table <- data.frame(transcript_sequences$ensembl_transcript_id, transcript_sequences$ensembl_gene_id, transcript_sequences$external_gene_name)
  translator_table <- translator_table[!is.na(translator_table$external_gene_name), ] # Remove the entries that do not have external gene names
  colnames(translator_table) <- c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name")

  return(list(transcript_sequences, translator_table)) # Returns the table with the ids and the sequences filtered (if applicable), and the translator table of the gene annotation formats
}

#' Extract Codon Composition of Sequences
#' @description
#' This function analyzes a given set of nucleotide sequences and computes the count of each codon present.
#'
#' @param sequences A \code{list} of nucleotide sequences (character strings) from which to extract the codon composition.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return Codon frequency per gene table of the \code{sequences}.
#' @export
#'
#' @examples
#' codon_composition <- ExtractCodonComposition(sequences = list("ATGCGTACG", "TTAAGGCCG"))

ExtractCodonComposition <- function(sequences, verbose = TRUE){

  ###
  # CALL: User or ObtainCodonFreqPerGene()
  # DESCRIPTION: This function takes nucleotide sequences and counts their codons to return a matrix with codons (rows) per genes/transcripts (columns).
  ###

  if (verbose) message("- Extracting the nucleotide sequences:")
  n <- if (!is.null(nrow(sequences))) nrow(sequences) else length(sequences) # Number of sequences to analyze
  pb <- utils::txtProgressBar(min = 0, max = n, style = 3) # Start a time counter to track the progress of the function execution

  for (i in 1:n) { # Iterates over each sequence

    if(!is.null(nrow(sequences))){ # Dealing with a table: (i) nucleotide sequence, (ii) sequence id
      transcript_id <- sequences[, 2][i] # Extract the sequence id
      sequence <- sequences[, 1][i] # Extract the actual nucleotide sequence
    } else { # Dealing with a list
      transcript_id <- paste("sequence", i, sep = "_") # Define a standard id for each sequence
      sequence <- sequences[i] # Retrieve the actual nucleotide sequence
    }

    # Ensure the validity of the sequence: (i) proper nucleotides, (ii) divided into triplets
    valid_nucleotides <- grepl("^[ATGCatgcnNUu]+$", sequence)
    if (isFALSE(valid_nucleotides)) stop(paste("Invalid sequence characters found in sequence:", transcript_id))

    # Generate the table with the counts - Count the number of appearances of each codon in each sequence
    if (nchar(sequence) %% 3 != 0)  stop("The sequence length is not a multiple of 3.")
    sequence <- substr(sequence, 1, (nchar(sequence) %/% 3) * 3) # Divide the nucleotide sequence into triplets
    codon_counts <- as.data.frame(table(substring(sequence, seq(1, nchar(sequence) - 2, 3), seq(3, nchar(sequence), 3))))
    colnames(codon_counts) <- c("codons", transcript_id) # Add the codons as a column not as rownames

    codon_freq_per_gene_matrix <- if (i == 1) codon_counts else merge(codon_freq_per_gene_matrix, codon_counts, all = TRUE)
    if (ncol(codon_freq_per_gene_matrix) == 0 || nrow(codon_freq_per_gene_matrix) == 0) stop("Incorrect dimensions.")

    utils::setTxtProgressBar(pb, i) # Increase the progress bar
  }

  close(pb) # Stop the progress bar
  if (verbose) message("- Extraction completed.")

  # Give proper format to the matrix
  rownames(codon_freq_per_gene_matrix) <- codon_freq_per_gene_matrix[, 1]
  codon_freq_per_gene_matrix[is.na(codon_freq_per_gene_matrix)] <- 0
  codon_freq_per_gene_matrix <- codon_freq_per_gene_matrix[, 2:ncol(codon_freq_per_gene_matrix)]

  codonsN <- grep(pattern = "N", x = rownames(codon_freq_per_gene_matrix)) # Removing undefined codons (if any)
  if (length(codonsN) != 0){
    if (verbose)  message("- There are codons with unknown bases (N nucleotides).\n", "- These codons will be removed.")
    codon_freq_per_gene_matrix <- codon_freq_per_gene_matrix[-codonsN, ]
  }
  return(codon_freq_per_gene_matrix) # Output the codon per gene matrix
}
