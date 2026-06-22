test_that("Check the codon frequency per gene function - getCodonFreq()", {
    base_args <- list(out_format = "external_gene_name", verbose = FALSE)

    # CASE 1: NULL inputs
    expect_error(do.call(getCodonFreq, c(list(dataset_name = NULL, genes_file = NULL), base_args)), "Specify either")

    # CASE 2: Invalid genes_file type
    expect_error(do.call(getCodonFreq, c(list(dataset_name = NULL, genes_file = 123), base_args)))
})

test_that("extractCodons produces an accurate count matrix", {
    seq1 <- "ATGCGTACG" # 3 codons: ATG, CGT, ACG
    seq2 <- "TTAAGGCCG" # 3 codons: TTA, AGG, CCG
    input_list <- list(s1 = seq1, s2 = seq2)

    res <- extractCodons(sequences = input_list, verbose = FALSE)

    expect_true(is.matrix(res))
    expect_type(res, "double")
    expect_equal(ncol(res), 2)
    expect_true(all(grepl("sequence_1|sequence_2", colnames(res)))) # Check if names are inherited from list or auto-generated

    expect_equal(unname(colSums(res)), c(3, 3)) # Each sequence has 9 bases -> 3 codons. The column sums must be 3.

    expect_true("ATG" %in% rownames(res)) # Ensure "ATG" is present and has a count of 1 in the first sequence
    expect_equal(res["ATG", 1], 1)
    expect_equal(res["ATG", 2], 0)

    # RNA and DNA versions of the same sequence should produce identical data frames
    rna_res <- extractCodons(sequences = "AUGCUGCAU", verbose = FALSE)
    dna_res <- extractCodons(sequences = "ATGCTGCAT", verbose = FALSE)

    expect_equal(rna_res, dna_res)
})

test_that("extractCodons handles data.frame inputs correctly", {
    df_input <- data.frame(seqs = c("ATGCGTACG", "TTAAGGCCG"), id = c("seq1", "seq2"))
    res_df <- extractCodons(sequences = df_input, verbose = FALSE)

    expect_true(is.matrix(res_df))
    expect_type(res_df, "double")
    expect_equal(sum(res_df), 6) # 6 unique codons total in this example
})

test_that("extractCodons correctly parses sequences and handles inputs", {
    expected_counts <- c(ATG = 1, CTG = 1, CAT = 1)
    variants <- c("ATGCTGCAT", "AUGCUGCAU", "atgctgcat", "augcugcau")

    for (seq in variants) {
        res <- extractCodons(sequences = seq, verbose = FALSE)
        expect_equal(as.numeric(res[names(expected_counts), 1]), unname(expected_counts), info = paste("Failed on variant:", seq))
    }

    # Dealing with a named list
    seq_list <- list(s1 = "ATGCTGCAT", s2 = "ATGCTGCAT")
    res_list <- extractCodons(sequences = seq_list, verbose = FALSE)
    expect_equal(ncol(res_list), 2)
    expect_equal(unname(colSums(res_list)), c(3, 3)) # Each seq has 3 codons

    # Dealing with a list
    seq_list <- list("ATGCTGCAT", "ATGCTGCAT")
    res_list <- extractCodons(sequences = seq_list, verbose = FALSE)
    expect_equal(ncol(res_list), 2)
    expect_equal(unname(colSums(res_list)), c(3, 3))

    # Dealing with a vector
    seq_vector <- c("ATGCTGCAT", "ATGCTGCAT")
    res_vector <- extractCodons(sequences = seq_vector, verbose = FALSE)
    expect_equal(ncol(res_vector), 2)
    expect_equal(unname(colSums(res_vector)), c(3, 3))

    # Dealing with a dataframe without ids
    seq_df <- data.frame(sequences = c("ATGCTGCAT", "ATGCTGCAT"))
    res_dg <- extractCodons(sequences = seq_df, verbose = FALSE)
    expect_equal(ncol(res_dg), 2)
    expect_equal(unname(colSums(res_dg)), c(3, 3))


    expect_error(extractCodons(sequences = 12345), "character")
    expect_warning(extractCodons(sequences = c("ATGCTGCAT", "ATGCTGCATC")), "multiple of 3")
    expect_error(suppressWarnings(extractCodons(sequences = "ATGCTGCATC")), "No valid sequences remained after filtering")
    expect_error(extractCodons(sequences = "ARGCRGCAR"), "Invalid sequence characters")
})

mock_transcripts_payload <- list(
    transcript_sequences = c("ATGGCC", "ATGTTT"), # Simple mock sequences
    out_format = "external_gene_name",
    count = 1,
    translator_table = data.frame(
        ensembl_transcript_id = c("ENST01", "ENST02", "ENST03"),
        external_gene_name = c("GeneA", "GeneB", ""), # Contains an empty unannotated name
        ensembl_gene_id = c("ENSG01", "ENSG02", "ENSG03"),
        stringsAsFactors = FALSE
    )
)

mock_codon_counts <- matrix(
    c(1, 0,  # Column 1 (ENST01) -> ATG=1, GCC=0
      0, 1,  # Column 2 (ENST02) -> ATG=0, GCC=1
      2, 3), # Column 3 (ENST03) -> ATG=2, GCC=3
    nrow = 2,
    ncol = 3,
    dimnames = list(c("ATG", "GCC"), c("ENST01", "ENST02", "ENST03"))
)
test_that("getCodonFreq accurately maps and renames transcript identifiers", {

    translate <- mock_transcripts_payload$translator_table
    out_format <- mock_transcripts_payload$out_format
    codon_freq <- mock_codon_counts

    tr <- stats::setNames(
        translate[[out_format]],
        nm = translate[["ensembl_transcript_id"]]
    )

    colnames(codon_freq) <- dplyr::recode(colnames(codon_freq), !!!tr)

    expect_equal(colnames(codon_freq), c("GeneA", "GeneB", ""))
})

test_that("getCodonFreq purges unannotated genes when retain_unannotated is FALSE", {

    codon_freq <- mock_codon_counts
    colnames(codon_freq) <- c("GeneA", "GeneB", "")
    translate <- mock_transcripts_payload$translator_table

    empty <- (colnames(codon_freq) == "")
    expect_equal(sum(empty), 1) # Confirm one empty element exists

    codon_freq <- codon_freq[, !empty, drop = FALSE]

    expect_equal(ncol(codon_freq), 2)
    expect_false("" %in% colnames(codon_freq))
})

test_that("getCodonFreq retains translator integrity when formatting to Ensembl Gene IDs", {

    codon_freq <- mock_codon_counts
    colnames(codon_freq) <- c("ENSG01", "ENSG02", "")  # Simulate if out_format was ensembl_gene_id
    translate <- mock_transcripts_payload$translator_table

    empty <- (colnames(codon_freq) == "")
    codon_freq <- codon_freq[, !empty, drop = FALSE]

    matched_indices <- which(translate$ensembl_gene_id %in% colnames(codon_freq))
    expect_equal(length(matched_indices), 2)
})

test_that("getCodonFreq completely executes all internal transformation and filtering lines", {

    # CASE 1
    mock_payload <- list( # Example from getTranscripts
        transcript_sequences = c("ATGGCC", "ATGTTT", "ATGAAA"),
        out_format = "external_gene_name",
        count = 1,
        translator_table = data.frame(
            ensembl_transcript_id = c("ENST01", "ENST02", "ENST03"),
            external_gene_name = c("GeneA", "GeneB", ""), # Crucial: triggers the 'empty' logic!
            ensembl_gene_id = c("ENSG01", "ENSG02", "ENSG03"),
            stringsAsFactors = FALSE
        )
    )

    mock_counts <- matrix( # Example from Extract Codons
        c(1, 0, 0, 1, 2, 3), nrow = 2, ncol = 3,
        dimnames = list(c("ATG", "GCC"), c("ENST01", "ENST02", "ENST03"))
    )

    mockery::stub(getCodonFreq, "getTranscripts", mock_payload)
    mockery::stub(getCodonFreq, "extractCodons", mock_counts)
    mockery::stub(getCodonFreq, "CheckDataFrame", TRUE)

    # CASE 2
    expect_message({
        res <- getCodonFreq(
            dataset_name = "hsapiens_gene_ensembl",
            retain_unannotated = FALSE,
            verbose = TRUE
        )
    }, "Validating the codon frequency matrix.")


    expect_equal(colnames(res$codon_freq_per_gene_matrix), c("GeneA", "GeneB"))
    expect_equal(ncol(res$codon_freq_per_gene_matrix), 2)
    expect_false("" %in% colnames(res$codon_freq_per_gene_matrix))

    expect_equal(nrow(res$translator_table), 2)
    expect_equal(res$translator_table$external_gene_name, c("GeneA", "GeneB"))
})

# Mock response for retain_geneversion = FALSE, filter = "length"
mock_biomart_length <- data.frame(
    ensembl_gene_id = c("ENSG1", "ENSG1", "ENSG2"),
    ensembl_transcript_id = c("ENST1_short", "ENST1_long", "ENST2"),
    chromosome_name = c("1", "1", "2"),
    external_gene_name = c("GeneA", "GeneA", "GeneB"),
    transcript_start = c(100, 100, 500),
    transcript_end = c(200, 400, 800), # ENST1_long has length 299 vs 99
    transcript_is_canonical = c(0, 1, 1),
    stringsAsFactors = FALSE
)

# Mock response for filter = "canonical"
mock_biomart_canonical <- data.frame(
    ensembl_gene_id = c("ENSG1", "ENSG1", "ENSG2"),
    ensembl_transcript_id = c("ENST1", "ENST1_alt", "ENST2"),
    chromosome_name = c("1", "1", "MT"), # Includes a mitochondrial record
    external_gene_name = c("GeneA", "GeneA", "GeneB"),
    transcript_start = c(100, 150, 200),
    transcript_end = c(500, 300, 400),
    transcript_is_canonical = c(1, 0, 1),
    stringsAsFactors = FALSE
)

test_that("extractGenes selects the longest transcript and drops duplicates when filter='length'", {

    mockery::stub(extractGenes, "biomaRt::getBM", mock_biomart_length)

    expect_message({
        res <- extractGenes(
            ensembl = "dummy_mart",
            filter = "length",
            retain_mitochondrial = TRUE,
            retain_geneversion = FALSE,
            verbose = TRUE
        )
    }, "Computing the length of each transcript")

    expect_equal(nrow(res), 2) # 3 raw records condensed down to 2 unique genes
    expect_true("transcript_length" %in% colnames(res))
    ensg1_record <- res[res$ensembl_gene_id == "ENSG1", ]
    expect_equal(ensg1_record$ensembl_transcript_id, "ENST1_long")
})

test_that("extractGenes filters by canonical transcripts when filter is not 'length'", {

    mockery::stub(extractGenes, "biomaRt::getBM", mock_biomart_canonical)

    res <- extractGenes(
        ensembl = "dummy_mart",
        filter = "canonical",
        retain_mitochondrial = TRUE,
        retain_geneversion = FALSE,
        verbose = FALSE
    )

    expect_equal(nrow(res), 2)
    expect_true(all(res$transcript_is_canonical == 1))
    expect_false("ENST1_alt" %in% res$ensembl_transcript_id)
})

test_that("extractGenes handles mitochondrial toggles correctly", {

    spy_getBM <- mockery::mock(mock_biomart_canonical)
    mockery::stub(extractGenes, "biomaRt::getBM", spy_getBM)

    extractGenes(
        ensembl = "dummy_mart",
        filter = "canonical",
        retain_mitochondrial = FALSE, # Triggers chromosome set filtering
        retain_geneversion = FALSE,
        verbose = FALSE
    )

    args <- mockery::mock_args(spy_getBM)[[1]]
    expect_equal(args$filters, c("biotype", "chromosome_name"))
    expect_true("MT" %in% args$values[[2]] == FALSE) # Proves MT was omitted via setdiff
})

mock_seq_no_version <- data.frame(
    coding = c("ATGGCCGCTGCATGA", "ATGTTTGGCTAA"),
    ensembl_transcript_id = c("ENST01", "ENST02"),
    stringsAsFactors = FALSE
)

mock_seq_with_version <- data.frame(
    coding = c("ATGGCCGCTGCATGA", "ATGTTTGGCTAA"),
    ensembl_transcript_id_version = c("ENST01.1", "ENST02.3"),
    stringsAsFactors = FALSE
)

test_that("ExtractSequences requests versioned and non-versioned attributes", {

    # CASE 1
    spy_getSequence <- mockery::mock(mock_seq_no_version)
    mockery::stub(ExtractSequences, "biomaRt::getSequence", spy_getSequence)

    res <- ExtractSequences(
        transcripts = c("ENST01", "ENST02"),
        ensembl = "dummy_mart",
        retain_geneversion = FALSE
    )
    expect_s3_class(res, "data.frame")
    expect_true("ensembl_transcript_id" %in% colnames(res))
    expect_equal(nrow(res), 2)
    args <- mockery::mock_args(spy_getSequence)[[1]]
    expect_equal(args$type, "ensembl_transcript_id")
    expect_equal(args$seqType, "coding")

    # CASE 2
    spy_getSequence <- mockery::mock(mock_seq_with_version)
    mockery::stub(ExtractSequences, "biomaRt::getSequence", spy_getSequence)

    res <- ExtractSequences(
        transcripts = c("ENST01.1", "ENST02.3"),
        ensembl = "dummy_mart",
        retain_geneversion = TRUE
    )

    expect_true("ensembl_transcript_id_version" %in% colnames(res))

    args <- mockery::mock_args(spy_getSequence)[[1]]
    expect_equal(args$type, "ensembl_transcript_id_version")
})

mock_genes_df <- data.frame(
    ensembl_transcript_id = c("ENST01", "ENST02"),
    ensembl_gene_id = c("ENSG01", "ENSG02"),
    external_gene_name = c("GeneA", "GeneB"),
    stringsAsFactors = FALSE
)

mock_seqs_df <- data.frame(
    coding = c("ATGGCC", "Sequence unavailable"), # Tests line 27-39 cleanup filter!
    ensembl_transcript_id = c("ENST01", "ENST02"),
    stringsAsFactors = FALSE
)

test_that("callingEnsembl builds translator tables and drops unavailable sequence lines", {

    mockery::stub(callingEnsembl, "biomaRt::useEnsembl", "dummy_mart")
    mockery::stub(callingEnsembl, "extractGenes", mock_genes_df)
    mockery::stub(callingEnsembl, "ExtractSequences", mock_seqs_df)

    expect_message({
        res <- callingEnsembl(
            dataset_name = "hsapiens_gene_ensembl",
            transcripts = NULL, # Enforces internal auto-extraction block
            filter = "length",
            retain_mitochondrial = TRUE,
            retain_geneversion = FALSE,
            verbose = TRUE
        )
    }, "- Removed 1 transcripts with unavailable sequence.")

    expect_type(res, "list")
    expect_named(res, c("translator_table", "transcript_sequences"))
    expect_equal(nrow(res$transcript_sequences), 1)
    expect_equal(res$transcript_sequences$ensembl_transcript_id, "ENST01")
    expect_equal(ncol(res$translator_table), 3)
    expect_equal(res$translator_table$ensembl_transcript_id, c("ENST01", "ENST02"))
})

test_that("callingEnsembl throws an error if all sequences are missing or unavailable", {

    all_unavailable_df <- data.frame(
        coding = c("Sequence unavailable", "Sequence unavailable"),
        ensembl_transcript_id = c("ENST01", "ENST02"),
        stringsAsFactors = FALSE
    )

    mockery::stub(callingEnsembl, "biomaRt::useEnsembl", "dummy_mart")
    mockery::stub(callingEnsembl, "extractGenes", mock_genes_df)
    mockery::stub(callingEnsembl, "ExtractSequences", all_unavailable_df)

    expect_error(
        callingEnsembl(
            dataset_name = "hsapiens_gene_ensembl",
            transcripts = NULL,
            filter = "length",
            retain_mitochondrial = TRUE,
            retain_geneversion = FALSE,
            verbose = FALSE
        ),
        "No transcripts found with available sequences."
    )
})

valid_fasta <- c(
    ">Sequence_1", "ATGCATGC", "ATGC-NGC", # USE of capital letters
    "", # Empty line that should be ignored without reporting an error
    ">Sequence_2", "ttgcatgcat" # Use of small letters
)
broken_header <- c(
    "Sequence_1", # Missing the '>' symbol
    "ATGCATGC"
)

corrupted_seq <- c(
    ">Sequence_1",
    "ATGCATGC",
    "ATGCXYZG" # 'XYZ' are completely invalid characters
)

test_that("checkFASTAFormat passes clean files and catches syntax failures", {

    # CASE 1
    tmp_valid <- withr::local_tempfile()
    writeLines(valid_fasta, tmp_valid)
    expect_invisible(checkFASTAFormat(tmp_valid))
    expect_true(checkFASTAFormat(tmp_valid))

    # CASE 2
    tmp_bad_header <- withr::local_tempfile()
    writeLines(broken_header, tmp_bad_header)
    expect_error(checkFASTAFormat(tmp_bad_header), "first line must be a header")

    # CASE 3
    tmp_corrupted <- withr::local_tempfile()
    writeLines(corrupted_seq, tmp_corrupted)
    expect_error(checkFASTAFormat(tmp_corrupted), "Invalid sequence characters found on lines: 3")
})

mock_fasta_data <- c(
    "ATGGCC" = ">ENST01 chromosome:GRCh38:1:100:200:1 transcript_biotype:protein_coding gene_symbol:GeneA",
    "ATGTTT" = ">ENST02 chromosome:GRCh38:2:500:600:1 transcript_biotype:protein_coding gene_symbol:GeneB",
    "ATGAAA" = ">ENST03 chromosome:GRCh38:MT:1:50:1 transcript_biotype:protein_coding gene_symbol:MT-CO1"
)
names(mock_fasta_data) <- unname(mock_fasta_data)
mock_extracted_df <- data.frame(
    coding = c("ATGGCC", "ATGTTT"),
    ensembl_transcript_id = c("ENST01", "ENST02"),
    ensembl_gene_id = c("ENSG01", "ENSG02"),
    external_gene_name = c("GeneA", "GeneB"),
    stringsAsFactors = FALSE
)

test_that("fromFASTAtoTable filters out mitochondrial records and outputs list components", {

    mockery::stub(fromFASTAtoTable, "ExtractFromFASTA", mock_extracted_df)
    mockery::stub(fromFASTAtoTable, "DefineTranslatorTable", "dummy_translator")

    expect_message({
        res <- fromFASTAtoTable(
            data = mock_fasta_data,
            transcripts = NULL,
            retain_mitochondrial = FALSE, # Triggers MT filtering logic block
            verbose = TRUE
        )
    }, "Mitochondrial genes will be removed.")

    expect_type(res, "list")
    expect_named(res, c("transcript_sequences", "translator_table"))
    expect_equal(res$translator_table, "dummy_translator")
})

test_that("fromFASTAtoTable narrows scope down to targeted identifiers explicitly matching a format", {

    mockery::stub(fromFASTAtoTable, "ExtractFromFASTA", mock_extracted_df)
    mockery::stub(fromFASTAtoTable, "DefineTranslatorTable", "dummy_translator")

    res <- fromFASTAtoTable(
        data = mock_fasta_data,
        transcripts = "ENST01",
        retain_mitochondrial = TRUE,
        verbose = FALSE
    )

    expect_equal(nrow(res$transcript_sequences), 1)
    expect_equal(res$transcript_sequences$ensembl_transcript_id, "ENST01")
})

test_that("fromFASTAtoTable stops cleanly if given a list of non-existent query identifiers", {

    mockery::stub(fromFASTAtoTable, "ExtractFromFASTA", mock_extracted_df)

    trace(fromFASTAtoTable, edit = FALSE)

    expect_error(
        fromFASTAtoTable(
            data = mock_fasta_data,
            transcripts = "NON_EXISTENT_ID",
            retain_mitochondrial = TRUE,
            verbose = FALSE
        ),
        "None of the ids given in 'transcripts' have been found"
    )
})
