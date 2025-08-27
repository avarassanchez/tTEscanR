#' Compute Codon Frequency-per-Gene Table
#' @description
#' This function computes codon usage frequencies for each gene based on a provided set of gene sequences.
#' It can optionally subset the analysis to specific transcripts and apply filtering criteria when multiple transcripts are available per gene.
#' Consequently, no `filter`parameter will be considered if parameters `transcripts` and `genes_file` are given.
#'
#' @param dataset_name A character string specifying the Ensembl species dataset name (e.g. \code{"hsapiens_gene_ensembl"}).
#' @param genes_file Optional; a path to a FASTA file.
#' @param transcripts Optional; a character vector of transcripts or gene IDs to subset the analysis.
#' @param filter Either \code{"canonical"} (default) or \code{"length"} (longest transcript) to specify which transcript to choose if several are available for the same gene.
#' @param retain.mitochondrial Logical; if \code{FALSE} filters outs the mitochondrial genes. Defaults to \code{FALSE}.
#' @param out_format Either \code{"external_gene_name"} (default), \code{"ensembl_gene_id"} or \code{"ensembl_transcript_id"} to specify annotation to use in the output codon frequency-per-gene table.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return Codon frequency-per-gene table and a translator gene annotation table (if available).
#' @export

ObtainCodonFreqPerGene <- function(dataset_name = NULL, genes_file = NULL, transcripts = NULL, filter = "canonical", retain.mitochondrial = FALSE, out_format = "external_gene_name", verbose = TRUE){

  ###
  # CALL: User
  # DESCRIPTION: This function gets a Ensembl organism name or a FASTA file and counts the codons in the genes or transcripts included.
  ###

  results.list <- list() # Define an empty list to store the results

  message("1 . Checking the format of the input data.")
  if ((!is.null(dataset_name) && !is.null(genes_file)) || (is.null(dataset_name) && is.null(genes_file))) stop("Please specify either `dataset_name` or `genes_files` to proceed.\n",
                                                                                                               "Note: Both input parameters can not be specified at the same time.")
  if (!(out_format %in% c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"))) stop("Please specify a suitable `output_format` parameter.\n",
                                                                                                   "Supported formats: ensembl_transcript_id, ensembl_gene_id, external_gene_name")
  count <- 2 # Count variable to keep the steps properly labelled

  if (!is.null(dataset_name)){ # The input is an Ensembl organism name - obtain the sequence of the genes or transcripts

    message("  1 . COMPLETED")
    if (is.null(transcripts) && (!(filter %in% c("length", "canonical"))) ) stop("Please specify `length` or `canonical` as `filter` input parameter.")
    ensembl_results <- CallingEnsembl(dataset_name = dataset_name, transcripts = transcripts, filter = filter, retain.mitochondrial = retain.mitochondrial, verbose = verbose) # Performs steps 2 and 3

    if (length(ensembl_results) == 2) { # There is a translator_table to add as a result
      transcript_sequences <- ensembl_results[[2]] # list of transcripts ids
      translator_table <- ensembl_results[[1]] # translator table: transcript ids and gene names (ensembl and short format)
      results.list <- append(results.list, list(translator_table))

    } else { # There is NO translator_table to add as a result
      transcript_sequences <- ensembl_results[[1]] # list of transcripts ids
      if (out_format != "ensembl_transcript_id") message("The out_format is forced to be an ensembl_transcript_id.") # translator_table can not be extracted when the transcripts are given by the user
      out_format <- "ensembl_transcript_id"
    }
    count <- 4
  } else { # The input is a FASTA file with the transcript sequences

    CheckFASTAFormat(genes_file) # Assessing the format of the input file
    genes_data <- Biostrings::readDNAStringSet(genes_file) # Read the sequences in the FASTA file

    FASTA_transformation <- FromFASTAtoTable(data = genes_data, transcripts = transcripts, retain.mitochondrial = retain.mitochondrial, verbose = verbose)
    message("  1 . COMPLETED")

    results.list <- append(results.list, list(FASTA_transformation[[2]]))
  }

  message(paste(as.character(count), ". Retrieving the codon composition of each transcript."))
  codon_freq_per_gene_matrix <- ExtractCodonComposition(sequences = FASTA_transformation[[1]])
  CheckDataFrame(data = codon_freq_per_gene_matrix, names = TRUE)

  if (out_format != "ensembl_transcript_id"){
    if (verbose) message(paste("- Translating the transcript identifiers to their corresponding", out_format, "format."))
    tr <- stats::setNames(FASTA_transformation[[2]][[out_format]], nm = FASTA_transformation[[2]][["ensembl_transcript_id"]])
    translated_transcripts <- dplyr::recode(colnames(codon_freq_per_gene_matrix), !!!tr)
    colnames(codon_freq_per_gene_matrix) <- translated_transcripts
  }
  message(paste("  ", as.character(count), ". COMPLETED"))
  count <- count + 1

  message(paste(as.character(count), ". Validating and returning the codon frequency per gene matrix."))
  CheckDataFrame(data = codon_freq_per_gene_matrix, names = TRUE)
  results.list <- append(results.list, list(codon_freq_per_gene_matrix))
  message(paste("  ", as.character(count), ". COMPLETED"))

  return(results.list) # Returns a list with the codon per gene table and a gene annotation translator table (if applicable).
}
