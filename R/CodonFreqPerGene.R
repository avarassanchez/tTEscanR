#' Compute Codon Frequency-per-Gene Table
#' @description
#' This function computes codon usage frequencies for each gene based on a provided set of gene sequences.
#' It can optionally subset the analysis to specific transcripts and apply filtering criteria when multiple transcripts are available for the same gene.
#' Consequently, no \code{filter} parameter will be considered if parameters \code{transcripts} and \code{genes_file} are given.
#'
#' @param dataset_name A character string specifying the Ensembl species dataset name (e.g. \code{"hsapiens_gene_ensembl"}).
#' @param genes_file Optional; a path to a FASTA file.
#' @param transcripts Optional; a character vector of transcripts or gene IDs to subset the analysis.
#' @param filter Either \code{"canonical"} (default) or \code{"length"} (longest transcript) to specify which transcript to choose if several are available for the same gene.
#' @param retain_mitochondrial Logical; if \code{FALSE} filters out the mitochondrial genes. Defaults to \code{FALSE}.
#' @param retain_unannotated Logical; if \code{FALSE} filters out the gene names that do not have an \code{"external_gene_name"} identifier. Defaults to \code{FALSE}.
#' @param retain_geneversion Logical; if \code{FALSE} retains the gene versions from the \code{"ensembl_gene_id"} identifier. Defaults to \code{TRUE}.
#' @param out_format Either \code{"external_gene_name"} (default), \code{"ensembl_gene_id"} or \code{"ensembl_transcript_id"} to specify annotation to use in the output codon frequency-per-gene table.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return Codon frequency-per-gene table and a translator gene annotation table (if available).
#' @export

GetCodonFreq <- function(dataset_name = NULL, genes_file = NULL, transcripts = NULL, filter = c("canonical", "length"), retain_mitochondrial = FALSE,
                         out_format = c("external_gene_name", "ensembl_transcript_id", "ensembl_gene_id"), retain_unannotated = FALSE, retain_geneversion = TRUE, verbose = TRUE){

  ###
  # CALL: User
  # DESCRIPTION: This function gets a Ensembl organism name or a FASTA file and counts the codons in the genes or transcripts included.
  ###

  message("1 . Checking the format of the input data.")
  if ((!is.null(dataset_name) && !is.null(genes_file)) || (is.null(dataset_name) && is.null(genes_file))) stop("Please specify either 'dataset_name' or 'genes_files' to proceed.\n",
                                                                                                               "Note: Both input parameters can not be specified at the same time.")

  out_format <- match.arg(out_format)
  count <- 2 # Count variable to keep the steps properly labelled

  if (!is.null(dataset_name)){ # The input is an Ensembl organism name - obtain the sequence of the genes or transcripts

    # if (is.null(transcripts) && (!(filter %in% c("length", "canonical"))) ) stop("Please specify `length` or `canonical` as `filter` input parameter.")
    if (is.null(transcripts)) filter <- match.arg(filter)
    message("1 . COMPLETED")

    ensembl_results <- CallingEnsembl(dataset_name = dataset_name, transcripts = transcripts, filter = filter, retain_mitochondrial = retain_mitochondrial,
                                      retain_geneversion = retain_geneversion, verbose = verbose) # Performs steps 2 and 3

    transcript_sequences <- ensembl_results$transcript_sequences # list of transcripts ids
    translator_table <- ensembl_results$translator_table # translator table: transcript ids and gene names (ensembl and short format)
    if (is.null(translator_table)) {
      if (out_format != "ensembl_transcript_id") message("- The 'out_format' will be an 'ensembl_transcript_id'.") # translator_table can not be extracted when the transcripts are given by the user
      out_format <- "ensembl_transcript_id"
    }

    count <- 4
  } else { # The input is a FASTA file with the transcript sequences

    CheckFASTAFormat(genes_file) # Assessing the format of the input file
    genes_data <- Biostrings::readDNAStringSet(genes_file) # Read the sequences in the FASTA file

    FASTA_transformation <- FromFASTAtoTable(data = genes_data, transcripts = transcripts, retain_mitochondrial = retain_mitochondrial, verbose = verbose)
    transcript_sequences <- FASTA_transformation$transcript_sequences
    translator_table <- FASTA_transformation$translator_table
    message("1 . COMPLETED")
  }

  message(paste(as.character(count), ". Retrieving the codon composition of each transcript."))
  codon_freq_per_gene_matrix <- ExtractCodons(sequences = transcript_sequences)
  CheckDataFrame(data = codon_freq_per_gene_matrix, required_names = TRUE)
  if (verbose) message(paste("- Number of protein coding transcripts:", ncol(codon_freq_per_gene_matrix)))

  if (out_format != "ensembl_transcript_id"){
    if (verbose) message(paste("- Translating the transcript identifiers to their corresponding", out_format, "format."))
    tr <- stats::setNames(translator_table[[out_format]], nm = translator_table[["ensembl_transcript_id"]])
    translated_transcripts <- dplyr::recode(colnames(codon_freq_per_gene_matrix), !!!tr)
    colnames(codon_freq_per_gene_matrix) <- translated_transcripts
  }
  message(paste(as.character(count), ". COMPLETED"))
  count <- count + 1

  message(paste(as.character(count), ". Validating and returning the codon frequency per gene matrix."))
  empty_names <- (colnames(codon_freq_per_gene_matrix) == "")
  if (any(empty_names)) {
    if (isFALSE(retain_unannotated)) {
      if (verbose) message(sprintf("- Removing %d genes without official external names.", sum(empty_names)))
      codon_freq_per_gene_matrix <- codon_freq_per_gene_matrix[ , !empty_names]
      translator_table <- translator_table[which(translator_table$external_gene_name %in% colnames(codon_freq_per_gene_matrix)), ]
    }
  }
  CheckDataFrame(data = codon_freq_per_gene_matrix, required_names = TRUE)
  message(paste(as.character(count), ". COMPLETED"))

  return(list(codon_freq_per_gene_matrix = codon_freq_per_gene_matrix,
              translator_table = translator_table)) # Returns a list with the codon per gene table and a gene annotation translator table (if applicable).
}
