test_that("featuresToAA stops when no overlapping features are found", {
    non_matching_codons <- c("XXX", "YYY", "ZZZ")
    expect_error(featuresToAA(data = non_matching_codons, notation_from = "codon", notation_to = "aa"), "No overlapping features identified")
})

test_that("featuresToAA correctly handles a data.frame with non-matching row names", {
    non_matching_df <- data.frame("value" = 1:5, row.names = c("XYZ", "ABC", "DEF", "GHI", "JKL"))
    expect_error(featuresToAA(data = non_matching_df, position = "row", notation_from = "codon", notation_to = "aa"), "No overlapping features identified")
})

test_that("featuresToAA correctly handles a data.frame with some matching row names", {
    partial_match_df <- data.frame("value" = 1:5, row.names = c("UUU", "GCA", "CGC", "XYZ", "ABC"))
    expect_no_error(translated_df <- featuresToAA(data = partial_match_df, position = "row", notation_from = "codon", notation_to = "aa"))

    expect_equal(nrow(translated_df), nrow(partial_match_df)) # Check that the number of rows is the same

    # Check if the row names have been translated for the matching features
    expect_true("F" %in% rownames(translated_df))
    expect_true("A" %in% rownames(translated_df))
    expect_true("R" %in% rownames(translated_df))
})

test_that("featuresToAA correctly handles a vector with some matching codons", {
    partial_match_vector <- c("TCT", "CGC", "GGC", "XYZ", "ABC")
    expect_no_error(translated_vector <- featuresToAA(data = partial_match_vector, notation_from = "codon", notation_to = "aa"))

    expect_equal(length(translated_vector), length(partial_match_vector)) # Check that the number of elements is the same

    # Check if the codons were translated correctly
    expect_true("S" %in% translated_vector)
    expect_true("R" %in% translated_vector)
    expect_true("G" %in% translated_vector)

    expect_true("XYZ" %in% translated_vector) # Check that the non-matching values were preserved
})

test_that("featuresToAA translates a data.frame column from codon to anticodon", {
    codon_df <- data.frame("gene_name" = paste0("gene_", 1:3), "codon" = c("GGC", "TCT", "CGC"), "value" = 1:3)
    expect_no_error(translated_df <- featuresToAA(data = codon_df, position = "codon", notation_from = "codon", notation_to = "anticodon"))

    expect_equal(nrow(translated_df), nrow(codon_df))
    expect_equal(translated_df$codon, c("GCC", "AGA", "GCG")) # Check that the 'codon' column has been correctly translated to anticodons
    expect_equal(translated_df$gene_name, codon_df$gene_name) # Ensure other columns remain unchanged
})

test_that("featuresToAA translates a data.frame from anticodon to codon", {
    anticodon_df <- data.frame("value" = 1:2, row.names = c("UUC", "CCG"))
    expect_no_error(translated_df <- featuresToAA(data = anticodon_df, position = "row", notation_from = "anticodon", notation_to = "codon"))

    expect_true("GAA" %in% rownames(translated_df))
    expect_true("CGG" %in% rownames(translated_df))
})

test_that("The features translation works - featuresToAA()", {
    data(default_tTEscanR_mRNA_data)
    tTEobject <- createObject(counts = default_tTEscanR_mRNA_data, assay = "mRNA")
    tTEobject <- computeCodonUsage(object = tTEobject, species = "hg38", additional_metrics = FALSE)
    codon_usage <- tTEobject@assays$CodonUsage

    # CASE 1: no error - base case using a vector of genes
    expect_no_error(featuresToAA(data = rownames(codon_usage), notation_from = "codon", notation_to = "aa"))
    expect_no_error(featuresToAA(data = rownames(codon_usage), notation_from = "anticodon", notation_to = "aa"))
    expect_no_error(featuresToAA(data = rownames(codon_usage), notation_from = "codon", notation_to = "anticodon"))
    expect_no_error(featuresToAA(data = rownames(codon_usage), notation_from = "anticodon", notation_to = "codon"))

    expect_no_error(featuresToAA(data = colnames(t(codon_usage)), notation_from = "codon", notation_to = "aa"))
    expect_no_error(featuresToAA(data = colnames(t(codon_usage)), notation_from = "anticodon", notation_to = "aa"))
    expect_no_error(featuresToAA(data = colnames(t(codon_usage)), notation_from = "codon", notation_to = "anticodon"))
    expect_no_error(featuresToAA(data = colnames(t(codon_usage)), notation_from = "anticodon", notation_to = "codon"))

    # CASE 1.1: no error - base case using a data table
    expect_no_error(featuresToAA(data = codon_usage, position = "row", notation_from = "codon", notation_to = "aa"))
    expect_no_error(featuresToAA(data = codon_usage, position = "row", notation_from = "anticodon", notation_to = "aa"))
    expect_no_error(featuresToAA(data = codon_usage, position = "row", notation_from = "codon", notation_to = "anticodon"))
    expect_no_error(featuresToAA(data = codon_usage, position = "row", notation_from = "anticodon", notation_to = "codon"))

    expect_no_error(featuresToAA(data = t(codon_usage), position = "column", notation_from = "codon", notation_to = "aa"))
    expect_no_error(featuresToAA(data = t(codon_usage), position = "column", notation_from = "anticodon", notation_to = "aa"))
    expect_no_error(featuresToAA(data = t(codon_usage), position = "column", notation_from = "codon", notation_to = "anticodon"))
    expect_no_error(featuresToAA(data = t(codon_usage), position = "column", notation_from = "anticodon", notation_to = "codon"))

    # CASE 1.2: error - base case using a data table (missing position parameter)
    expect_error(featuresToAA(data = codon_usage, notation_from = "codon", notation_to = "aa"))
    expect_error(featuresToAA(data = codon_usage, notation_from = "anticodon", notation_to = "aa"))
    expect_error(featuresToAA(data = codon_usage, notation_from = "codon", notation_to = "anticodon"))
    expect_error(featuresToAA(data = codon_usage, notation_from = "anticodon", notation_to = "codon"))

    expect_error(featuresToAA(data = t(codon_usage), notation_from = "codon", notation_to = "aa"))
    expect_error(featuresToAA(data = t(codon_usage), notation_from = "anticodon", notation_to = "aa"))
    expect_error(featuresToAA(data = t(codon_usage), notation_from = "codon", notation_to = "anticodon"))
    expect_error(featuresToAA(data = t(codon_usage), notation_from = "anticodon", notation_to = "codon"))

    # CASE 2: error - same input and output format using a vector of genes
    expect_error(featuresToAA(data = rownames(codon_usage), notation_from = "codon", notation_to = "codon"))
    expect_error(featuresToAA(data = rownames(codon_usage), notation_from = "anticodon", notation_to = "anticodon"))
    expect_error(featuresToAA(data = colnames(t(codon_usage)), notation_from = "codon", notation_to = "codon"))
    expect_error(featuresToAA(data = colnames(t(codon_usage)), notation_from = "anticodon", notation_to = "anticodon"))

    # CASE 2.1: error - same input and output format using a data table
    expect_error(featuresToAA(data = codon_usage, position = "row", notation_from = "codon", notation_to = "codon"))
    expect_error(featuresToAA(data = codon_usage, position = "row", notation_from = "anticodon", notation_to = "anticodon"))
    expect_error(featuresToAA(data = t(codon_usage), position = "column", notation_from = "codon", notation_to = "codon"))
    expect_error(featuresToAA(data = t(codon_usage), position = "column", notation_from = "anticodon", notation_to = "anticodon"))

    # CASE 3: error - wrong usage of the parameters using a vector of genes
    expect_error(featuresToAA(data = codon_usage, position = "row", notation_from = "aa", notation_to = "codon")) # wrong notation_to parameter
})
