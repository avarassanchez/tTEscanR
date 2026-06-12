ExtractGenes <- function(
    ensembl, filter, retain_mitochondrial, retain_geneversion, verbose = TRUE
) {
    if (verbose) {
        message(
            "- Retrieving protein coding genes & ",
            "filtering out mitochondrial genes (if applicable)."
        )
    }
    if (isTRUE(retain_geneversion)) {
        id_names <- c(
            "ensembl_gene_id_version", "ensembl_transcript_id_version"
        )
    } else {
        id_names <- c("ensembl_gene_id", "ensembl_transcript_id")
    }
    if (isFALSE(retain_mitochondrial)) { # Exclude mitochondrial genes
        filters_set <- c("biotype", "chromosome_name")
        values_set <- list("protein_coding", setdiff(
            c(seq_len(22), "X", "Y"), "MT"
        ))
    } else { # Retain all protein coding genes
        filters_set <- "biotype"
        values_set <- "protein_coding"
    }
    coding <- biomaRt::getBM(attributes = c(
            id_names, "chromosome_name", "external_gene_name",
            "transcript_start", "transcript_end", "transcript_is_canonical"
        ), filters = filters_set, values = values_set, mart = ensembl
    )
    name <- if (filter == "length") "largest" else "canonical"
    if (verbose) {
        message(
            "- Selecting (if many) the ", name, " transcript for each gene.")
    }
    if (filter == "length") { # If many transcripts uses 'filter' to select
        if (verbose) message("- Computing the length of each transcript.")
        coding$transcript_length <- as.numeric(
            abs((coding$transcript_start + 1) - coding$transcript_end)
        )
        coding <- coding[order(
            coding$ensembl_gene_id, -coding$transcript_length
        ), ]
        coding <- coding[!duplicated(coding$ensembl_gene_id), ]
    } else {
        coding <- subset(coding, coding$transcript_is_canonical == 1)
    }
    return(coding) # Returns a table with the features retrieved by getBM()
}

ExtractSequences <- function(transcripts, ensembl, retain_geneversion) {
    if (isTRUE(retain_geneversion)) {
        id_attribute <- "ensembl_transcript_id_version"
    } else {
        id_attribute <- "ensembl_transcript_id"
    }

    transcript_sequences <- biomaRt::getSequence(
        id = transcripts, type = id_attribute, seqType = "coding",
        mart = ensembl
    )

    # Returns a table with: (i) transcript id, and (ii) nucleotide sequence
    return(transcript_sequences)
}

CallingEnsembl <- function(
    dataset_name, transcripts, filter, retain_mitochondrial,
    retain_geneversion, verbose
) {
    trans <- NULL
    if (verbose) message("2 . Retrieving the Ensembl dataset.") # "ensembl"
    ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = dataset_name)
    if (is.null(transcripts)) {
        if (verbose) message("- Retrieving the protein-coding transcripts.")
        transcripts <- ExtractGenes( # Get table with features of transcripts
            ensembl = ensembl, retain_geneversion = retain_geneversion,
            retain_mitochondrial = retain_mitochondrial, filter = filter
        )
        trans <- data.frame( # Define table with the different ids
            transcripts$ensembl_transcript_id, transcripts$ensembl_gene_id,
            transcripts$external_gene_name
        )
        colnames(trans) <- c(
            "ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"
        )
        transcripts <- transcripts$ensembl_transcript_id # Retrieve ids
    }
    if (verbose) {
        message("2 . COMPLETED\n3 . Extracting the genomic sequences.")
    }
    seq <- ExtractSequences( # Based on id, retrieve sequence
        transcripts = transcripts, ensembl = ensembl,
        retain_geneversion = retain_geneversion
    )
    if (nrow(seq) == 0) {
        stop("No transcripts found.\nCheck 'transcripts' if applicable.")
    }
    unavailable <- which(seq$coding == "Sequence unavailable")
    len <- length(unavailable)
    if (len != 0) {
        if (verbose) {
            message("- Removed", len, " transcripts with unavailable sequence.")
        }
        if (len == length(seq)) {
            stop(
                "No transcripts found with available sequences.\n",
                "Check 'transcripts' if applicable."
            )
        }
        seq <- seq[-unavailable, ]
    }
    if (verbose) message("- Number protein coding transcripts: ", nrow(seq))
    if (verbose) message("3 . COMPLETED")
    return(list(translator_table = trans, transcript_sequences = seq))
}

CheckFASTAFormat <- function(file) {
    # Perform a check considering a subset of the file so that not all is loaded
    lines <- readLines(file, n = 100, warn = FALSE)

    header <- grepl("^>", lines)
    if (!grepl("^>", lines[1])) {
        stop("Not a valid FASTA file. The first line must be a header ('>').")
    }
    if (length(lines) < 2) {
        stop("File too short or empty. Not a valid FASTA file.")
    }

    invalid_lines <- which(!header && !grepl("^[ATGCatgcnNUu]+$", lines))
    if (length(invalid_lines) > 0) {
        stop(
            "Invalid sequence characters found on lines:",
            paste(invalid_lines, collapse = ", ")
        )
    }

    return(invisible(TRUE))
}

FromFASTAtoTable <- function(data, transcripts, retain_mitochondrial, verbose) {
    if (verbose) {
        message("2. Extracting the genomic sequence of each transcript.")
        message("- Retrieving the protein-coding transcripts.")
    }
    prot_data <- data[grep("transcript_biotype:protein_coding", names(data))]
    tran_seq <- data.frame( # Generate table: ids & nucleotide sequences
        coding = unname(as.character(prot_data)),
        ensembl_transcript_id = names(prot_data), stringsAsFactors = FALSE
    )
    if (nrow(tran_seq) == 0) stop("There are no protein-coding transcripts.")
    if (isFALSE(retain_mitochondrial)) { # Filter out mitochondrial genes
        if (verbose) message("- Mitochondrial genes will be removed.")
        mitochondrial_index <- grepl(
            "chromosome:[^:]+:MT|gene_symbol:MT-|mitochondrial",
            tran_seq$ensembl_transcript_id, ignore.case = TRUE
        )
        if (all(mitochondrial_index)) {
            stop("All protein-coding transcripts are mitochondrial.")
        }
        tran_seq <- tran_seq[!mitochondrial_index, ]
    }
    tran_seq <- ExtractFromFASTA(transcript_seq = tran_seq)
    if (!is.null(transcripts)) { # Targeted approach based on the transcripts
        if (verbose) message("- Trimming by the ids included as 'transcripts'.")
        targeted_indx <- list(
            gene_id = which(tran_seq$ensembl_gene_id %in% transcripts),
            transcript_id = which(
                tran_seq$ensembl_transcript_id %in% transcripts
            ),
            gene_name = which(tran_seq$external_gene_name %in% transcripts)
        )
        if ((lengths(targeted_indx) == 0)) {
            stop(
                "None of the ids given in 'transcripts' have been found.\n",
                "Supported formats: Ensembl transcript id, Ensembl gene id ",
                "and external gene name."
            )
        }
        # Select column of the ids with a higher match with the input ids list
        targeted_indx <- targeted_indx[[which.max(lengths(targeted_indx))]]
        tran_seq <- tran_seq[targeted_indx, ]
    }
    if (nrow(tran_seq) == 0) stop("There are no protein-coding transcripts.")
    if (verbose) {
        message("- Number of protein-coding transcripts: ", nrow(tran_seq))
    }
    trans <- DefineTranslatorTable(tran_seq)
    return(list(transcript_sequences = tran_seq, translator_table = trans))
}

DefineTranslatorTable <- function(transcript_sequences) {
    # Generate a translator table for the gene ids
    translator_table <- transcript_sequences[, c(
        "ensembl_transcript_id",
        "ensembl_gene_id",
        "external_gene_name"
    )]

    # Remove entries without external gene names
    translator_table <- translator_table[
        !is.na(translator_table$external_gene_name),
    ]

    return(translator_table)
}

ExtractFromFASTA <- function(transcript_seq) {
    # The FASTA file contains all the info - retrieve each part independently
    transcript_seq <- data.frame(
        coding = transcript_seq$coding,
        ensembl_gene_id = stringr::str_extract(
            transcript_seq$ensembl_transcript_id, "ENSG[0-9]+\\.?[0-9]*"
        ),
        ensembl_transcript_id = stringr::str_extract(
            transcript_seq$ensembl_transcript_id, "ENST[0-9]+\\.?[0-9]*"
        ),
        external_gene_name = stringr::str_extract(
            transcript_seq$ensembl_transcript_id, "(?<=gene_symbol:)[^ ]+"
        ),
        transcript_length = nchar(transcript_seq$coding),
        stringsAsFactors = FALSE
    )

    return(transcript_seq = transcript_seq)
}
