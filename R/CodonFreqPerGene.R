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
#' dataset name (e.g. \code{"hsapiens_gene_ensembl"}).
#' @param genes_file Optional; a path to a FASTA file.
#' @param transcripts Optional; a character vector of transcripts or gene IDs
#' to subset the analysis.
#' @param filter Either \code{"canonical"} or \code{"length"}
#' (longest transcript) to specify which transcript to choose if several are
#' available for the same gene.
#' @param retain_mitochondrial Logical; if \code{FALSE} filters out the
#' mitochondrial genes. Defaults to \code{FALSE}.
#' @param retain_unannotated Logical; if \code{FALSE} filters out the gene
#' names that do not have an \code{"external_gene_name"} identifier. Defaults
#' to \code{FALSE}.
#' @param retain_geneversion Logical; if \code{FALSE} retains the gene versions
#' from the \code{"ensembl_gene_id"} identifier. Defaults to \code{TRUE}.
#' @param out_format Either \code{"external_gene_name"},
#' \code{"ensembl_gene_id"} or \code{"ensembl_transcript_id"} to specify
#' annotation to use in the output codon frequency-per-gene table.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#' Defaults to \code{TRUE}.
#'
#' @return Codon frequency-per-gene table and a translator gene annotation
#' table (if available).
#' @export
#'
GetCodonFreq <- function(
    dataset_name = NULL, genes_file = NULL, verbose = TRUE, transcripts = NULL,
    retain_mitochondrial = FALSE, filter = c("canonical", "length"),
    retain_unannotated = FALSE, retain_geneversion = TRUE, out_format = c(
        "external_gene_name", "ensembl_transcript_id", "ensembl_gene_id"
    )
) {
    extract_data <- GetTranscripts(
        dataset_name = dataset_name, transcripts = transcripts, filter = filter,
        out_format = out_format, verbose = verbose, genes_file = genes_file,
        retain_version = retain_geneversion, retain_mt = retain_mitochondrial
    )
    translate <- extract_data$translator_table
    out_format <- extract_data$out_format
    count <- extract_data$count # Variable to keep the steps properly labelled
    if (verbose) {
        message(as.character(count), ". Analyzing the codon composition.")
    }
    codon_freq <- ExtractCodons(sequences = extract_data$transcript_sequences)
    CheckDataFrame(data = codon_freq, required_names = TRUE)
    if (verbose) message("- Protein-coding transcripts: ", ncol(codon_freq))
    if (out_format != "ensembl_transcript_id") {
        if (verbose) message("- Changing transcript ids format to ", out_format)
        tr <- stats::setNames(translate[[out_format]],
            nm = translate[["ensembl_transcript_id"]]
        )
        colnames(codon_freq) <- dplyr::recode(colnames(codon_freq), !!!tr)
    }
    if (verbose) message(as.character(count), ". COMPLETED")
    count <- count + 1
    if (verbose) {
        message(as.character(count), ". Validating the codon frequency matrix.")
    }
    empty <- (colnames(codon_freq) == "")
    if (any(empty) && isFALSE(retain_unannotated)) {
        if (verbose) {
            message("- Removing ", sum(empty), " genes without external names.")
        }
        codon_freq <- codon_freq[, !empty]
        translate <- translate[
            which(translate$external_gene_name %in% colnames(codon_freq)),
        ]
    }
    CheckDataFrame(data = codon_freq, required_names = TRUE)
    if (verbose) message(as.character(count), ". COMPLETED")
    return(list(
        codon_freq_per_gene_matrix = codon_freq, translator_table = translate
    ))
}

GetTranscripts <- function(
    dataset_name, transcripts, filter = c("canonical", "length"),
    retain_mt, verbose, retain_version, genes_file, out_format = c(
        "external_gene_name", "ensembl_transcript_id", "ensembl_gene_id"
    )
) {
    if (verbose) message("1 . Checking the format of the input data.")
    if ((!is.null(dataset_name) && !is.null(genes_file)) ||
        (is.null(dataset_name) && is.null(genes_file))) {
        stop(
            "Specify either 'dataset_name' or 'genes_files' to proceed.",
            "\nNote: Both input parameters can not be specified together."
        )
    }
    out_format <- match.arg(out_format)
    if (!is.null(dataset_name)) { # Obtain sequence of the genes/transcripts
        if (is.null(transcripts)) filter <- match.arg(filter)
        if (verbose) message("1 . COMPLETED")
        ensembl_results <- CallingEnsembl(
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
        CheckFASTAFormat(genes_file) # Assessing the format of the input file
        genes_data <- Biostrings::readDNAStringSet(genes_file) # Get FASTA seqs
        FASTA_transformation <- FromFASTAtoTable(
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
#' from which to extract the codon composition.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#' Defaults to \code{TRUE}.
#'
#' @return Codon frequency per gene table of the \code{sequences}.
#' @export
#'
#' @examples
#' codon_composition <- ExtractCodons(sequences = list(
#'     "ATGCGTACG",
#'     "TTAAGGCCG"
#' ))
ExtractCodons <- function(sequences, verbose = TRUE) {
    if (verbose) message("- Extracting the nucleotide sequences:")
    bases <- c("A", "C", "G", "T")
    all_64_codons <- apply(
        expand.grid(bases, bases, bases), 1, paste, collapse = ""
    )

    extract_values <- HelperExtractCodons(
        sequences = sequences, all_64_codons = all_64_codons
    )
    counts_list <- extract_values$counts_list
    transcript_names <- extract_values$transcript_names
    if (verbose) message("\n- Extraction completed.\n- Assembling matrix...")

    # Filter out the skipped sequences
    valid_indices <- !vapply(counts_list, is.null, FUN.VALUE = numeric(1))
    if (!any(valid_indices)) {
        stop("No valid sequences remained after filtering.")
    }

    # Give proper format to the matrix (do.call cbind is incredibly fast here!)
    codon_freq_per_gene_matrix <- do.call(cbind, counts_list[valid_indices])
    rownames(codon_freq_per_gene_matrix) <- all_64_codons
    colnames(codon_freq_per_gene_matrix) <- transcript_names[valid_indices]

    return(codon_freq_per_gene_matrix)
}

HelperExtractCodons <- function(sequences, all_64_codons) {
    is_tabular <- is.data.frame(sequences) || is.matrix(sequences)
    n <- if (is_tabular) nrow(sequences) else length(sequences)
    transcript_names <- character(n)
    counts_list <- vector("list", n)

    pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
    for (i in seq_len(n)) { # Iterates over each sequence
        if (is_tabular) {
            sequence <- as.character(sequences[i, 1])
            if (ncol(sequences) >= 2) {
                transcript_id <- as.character(sequences[i, 2])
            } else {
                transcript_id <- paste("sequence", i, sep = "_")
            }
        } else {
            sequence <- as.character(sequences[i])
            transcript_id <- paste("sequence", i, sep = "_")
        }
        transcript_names[i] <- transcript_id

        # Ensure the validity of the sequence
        sequence <- toupper(sequence)
        sequence <- gsub("U", "T", sequence)
        valid_nucleotides <- grepl("^[ATGCN]+$", sequence)
        if (isFALSE(valid_nucleotides)) {
            stop(
                "Invalid sequence characters found in sequence: ", transcript_id
            )
        }
        if (nchar(sequence) %% 3 != 0) { # Strict multiple of 3 check
            warning(sprintf(
                "Transcript %s length is not a multiple of 3. Skipping.",
                transcript_id
            ))
            utils::setTxtProgressBar(pb, i)
            next # Skip this sequence, but don't stop the whole function
        }
        # Generate the table with the counts
        triplets <- substring(
            sequence, seq(1, nchar(sequence) - 2, 3), seq(3, nchar(sequence), 3)
        )
        codon_counts <- table(factor(triplets, levels = all_64_codons))
        counts_list[[i]] <- as.numeric(codon_counts)
        utils::setTxtProgressBar(pb, i)
    }
    close(pb)
    return(list(counts_list = counts_list, transcript_names = transcript_names))
}
