#' Compute Codon Frequency-per-Gene Table
#' @description
#' This function computes codon usage frequencies for each gene based on a
#' provided set of gene sequences. It can optionally subset the analysis to
#' specific transcripts and apply filtering criteria when multiple transcripts
#' are available for the same gene. Consequently, no \code{filter} parameter
#' will be considered if parameters \code{transcripts} and \code{genes_file}
#' are given.
#'
#' @param dataset_name A character string specifying the Ensembl species
#'     dataset name (e.g. \code{"hsapiens_gene_ensembl"}).
#' @param genes_file Optional; a path to a FASTA file.
#' @param transcripts Optional; a character vector of transcripts or gene IDs
#'     to subset the analysis.
#' @param filter Either \code{"canonical"} or \code{"length"}
#'     (longest transcript) to specify which transcript to choose if several
#'     are available for the same gene.
#' @param retain_mitochondrial Logical; if \code{FALSE} filters out the
#'     mitochondrial genes. Defaults to \code{FALSE}.
#' @param retain_unannotated Logical; if \code{FALSE} filters out the gene
#'     names that do not have an \code{"external_gene_name"} identifier.
#'     Defaults to \code{FALSE}.
#' @param retain_geneversion Logical; if \code{FALSE} retains the gene versions
#'     from the \code{"ensembl_gene_id"} identifier. Defaults to \code{TRUE}.
#' @param out_format Either \code{"external_gene_name"},
#'     \code{"ensembl_gene_id"} or \code{"ensembl_transcript_id"} to specify
#'     annotation to use in the output codon frequency-per-gene table.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#'     Defaults to \code{TRUE}.
#'
#' @return Codon frequency-per-gene table and a translator gene annotation
#'     table (if available).
#' @export
#'
#' @examples
#' \donttest{
#' codon_freq_results_canonical_hg38 <- getCodonFreq(
#'     dataset_name = "hsapiens_gene_ensembl", filter = "canonical",
#'     retain_geneversion = TRUE, retain_mitochondrial = FALSE,
#'     out_format = "external_gene_name"
#' )
#' }
getCodonFreq <- function(dataset_name = NULL, genes_file = NULL, verbose = TRUE,
    transcripts = NULL, retain_mitochondrial = FALSE,
    filter = c("canonical", "length"), retain_unannotated = FALSE,
    retain_geneversion = TRUE, out_format = c(
        "external_gene_name", "ensembl_transcript_id", "ensembl_gene_id")) {
    extract_data <- getTranscripts(
        dataset_name = dataset_name, transcripts = transcripts,
        filter = match.arg(filter), out_format = match.arg(out_format),
        verbose = verbose, genes_file = genes_file,
        retain_version = retain_geneversion, retain_mt = retain_mitochondrial
    )
    trans <- extract_data$translator_table
    out_format <- extract_data$out_format
    count <- extract_data$count # Variable to keep the steps properly labelled
    if (verbose) message(count, ". Analyzing the codon composition.")
    freq <- extractCodons(sequences = extract_data$transcript_sequences)
    checkDataFrame(data = freq, required_names = TRUE)
    if (verbose) message("- Protein-coding transcripts: ", ncol(freq))
    if (out_format != "ensembl_transcript_id" && !is.null(trans)) {
        if (verbose) message("- Changing transcript ids format to ", out_format)
        tr_map <- stats::setNames(
            as.character(trans[[out_format]]),
            nm = as.character(trans[["ensembl_transcript_id"]])
        )
        idx <- match(colnames(freq), names(tr_map))
        matched <- !is.na(idx)
        new_colnames <- colnames(freq)
        new_colnames[matched] <- tr_map[idx[matched]]
        colnames(freq) <- new_colnames
    }
    if (verbose) message(count, ". COMPLETED")
    count <- count + 1
    if (verbose) message(count, ". Validating the codon frequency matrix.")
    empty <- is.na(colnames(freq)) | colnames(freq) == ""
    if (any(empty) && !retain_unannotated) {
        if (verbose) {
            message("- Removing ", sum(empty), " genes without external names.")
        }
        freq <- freq[, !empty, drop = FALSE]
        if (!is.null(trans) && "external_gene_name" %in% colnames(trans)) {
            trans <- trans[
                which(trans$external_gene_name %in% colnames(freq)), ,
                drop = FALSE
            ]
        }
    }
    checkDataFrame(data = freq, required_names = TRUE)
    if (verbose) message(count, ". COMPLETED")
    return(list(codon_freq_per_gene_matrix = freq, translator_table = trans))
}

getTranscripts <- function(dataset_name, transcripts, retain_mt,
    filter, verbose, retain_version, genes_file, out_format) {
    if (verbose) message("1 . Checking the format of the input data.")
    if (is.null(dataset_name) == is.null(genes_file))  {
        stop(
            "Specify either 'dataset_name' or 'genes_files' to proceed.",
            "\nNote: Both input parameters can not be specified together."
        )
    }

    if (!is.null(dataset_name)) { # Obtain sequence of the genes/transcripts
        if (verbose) message("1 . COMPLETED")
        ensembl_results <- callingEnsembl(
            dataset_name = dataset_name, transcripts = transcripts,
            filter = filter, retain_mitochondrial = retain_mt,
            verbose = verbose, retain_geneversion = retain_version
        ) # Performs steps 2 and 3
        transcript_seq <- ensembl_results$transcript_seq # List transcripts ids
        translator_table <- ensembl_results$translator_table # Ids and names

        if (is.null(translator_table)) {
            if (out_format != "ensembl_transcript_id" && verbose) {
                message("- The 'out_format' will be 'ensembl_transcript_id'.")
            }
            out_format <- "ensembl_transcript_id"
        }
        count <- 4

    } else { # The input is a FASTA file with the transcript sequences
        checkFASTAFormat(genes_file) # Assessing the format of the input file
        genes_data <- Biostrings::readDNAStringSet(genes_file) # Get FASTA seqs
        FASTA_transformation <- fromFASTAtoTable(
            data = genes_data, transcripts = transcripts,
            retain_mitochondrial = retain_mt, verbose = verbose
        )
        transcript_seq <- FASTA_transformation$transcript_seq
        translator_table <- FASTA_transformation$translator_table
        if (verbose) message("1 . COMPLETED")
        count <- 2
    }
    return(list(
        transcript_sequences = transcript_seq, out_format = out_format,
        translator_table = translator_table, count = count
    ))
}

#' Extract Codon Composition of Sequences
#' @description
#' This function analyzes a given set of nucleotide sequences and computes the
#' count of each codon present.
#'
#' @param sequences A \code{list} of nucleotide sequences (character strings)
#'     from which to extract the codon composition.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#'     Defaults to \code{TRUE}.
#'
#' @return Codon frequency per gene table of the \code{sequences}.
#' @export
#'
#' @examples
#' codon_composition <- extractCodons(sequences = list(
#'     "ATGCGTACG",
#'     "TTAAGGCCG"
#' ))
extractCodons <- function(sequences, verbose = TRUE) {
    if (verbose) message("- Extracting the nucleotide sequences:")
    bases <- c("A", "C", "G", "T")
    grid <- expand.grid(bases, bases, bases, KEEP.OUT.ATTRS = FALSE)
    all_64_codons <- do.call(paste0, grid)

    codon_freq_per_gene_matrix <- helperExtractCodons(
        sequences = sequences, all_64_codons = all_64_codons
    )

    if (ncol(codon_freq_per_gene_matrix) == 0) {
        stop("No valid sequences remained after filtering.")
    }

    if (verbose) message("\n- Extraction completed.")

    return(codon_freq_per_gene_matrix)
}

helperExtractCodons <- function(sequences, all_64_codons) {
    is_tabular <- is.data.frame(sequences) || is.matrix(sequences)
    n <- if (is_tabular) nrow(sequences) else length(sequences)
    if (is_tabular) {
        seq_vec <- as.character(sequences[, 1])
        transcript_names <- if (ncol(sequences) >= 2) {
            as.character(sequences[, 2])
        } else {
            paste("sequence", seq_len(n), sep = "_")
        }
    } else {
        seq_vec <- as.character(sequences)
        transcript_names <- paste("sequence", seq_len(n), sep = "_")
    }

    ## Ensure the validity of the sequence
    seq_vec <- chartr("U", "T", toupper(seq_vec))
    invalid_chars <- !grepl("^[ATGCN]*$", seq_vec)
    if (any(invalid_chars)) {
        stop(
            "Invalid sequence characters found in sequence: ",
            transcript_names[invalid_chars][1]
        )
    }
    seq_lens <- nchar(seq_vec)
    invalid_len <- (seq_lens %% 3 != 0) | (seq_lens == 0)
    if (any(invalid_len)) {
        warning(sprintf(
                "Skipping transcripts with length not a multiple of 3: %s",
                paste0(transcript_names[invalid_len], collapse = ","))
        )
    }
    valid_mask <- !invalid_len
    if (!any(valid_mask)) {
        return(
            matrix(
                numeric(0), nrow = 64, ncol = 0,
                dimnames = list(all_64_codons, character(0))
            )
        )
    }
    dna_seqs <- Biostrings::DNAStringSet(seq_vec[valid_mask])
    counts_mat <- Biostrings::trinucleotideFrequency(dna_seqs, step = 3)
    counts_mat <- t(counts_mat[, all_64_codons, drop = FALSE])
    colnames(counts_mat) <- transcript_names[valid_mask]
    return(counts_mat)
}
