test_that("featuresToAA fails when no overlapping features exist", {
    ## CASE 1: Non-matching vector input
    expect_error(featuresToAA(data = c("XXX", "YYY", "ZZZ"), notation_from = "codon", notation_to = "aa"), "No overlapping features identified")

    ## CASE 2: Non-matching data frame row names
    non_matching_df <- data.frame(value = 1:5, row.names = c("XYZ", "ABC", "DEF", "GHI", "JKL"))
    expect_error(featuresToAA(data = non_matching_df, position = "row", notation_from = "codon", notation_to = "aa"), "No overlapping features identified")
})

test_that("featuresToAA translates partial matching vectors and data frames", {
    ## CASE 1: Partial match data frame
    partial_df <- data.frame(value = 1:5, row.names = c("UUU", "GCA", "CGC", "XYZ", "ABC"))
    expect_no_error(translated_df <- featuresToAA(data = partial_df, position = "row", notation_from = "codon", notation_to = "aa"))
    expect_equal(nrow(translated_df), nrow(partial_df))
    expect_true(all(c("F", "A", "R") %in% rownames(translated_df)))

    ## CASE 2: Partial match vector
    partial_vec <- c("TCT", "CGC", "GGC", "XYZ", "ABC")
    expect_no_error(translated_vec <- featuresToAA(data = partial_vec, notation_from = "codon", notation_to = "aa"))
    expect_equal(length(translated_vec), length(partial_vec))
    expect_true(all(c("S", "R", "G", "XYZ") %in% translated_vec))
})

test_that("featuresToAA handles specific data frame column and anticodon inputs", {
    ## CASE 1: Column-based translation (codon to anticodon)
    codon_df <- data.frame(gene_name = paste0("gene_", 1:3), codon = c("GGC", "TCT", "CGC"), value = 1:3)
    expect_no_error(translated_col <- featuresToAA(data = codon_df, position = "codon", notation_from = "codon", notation_to = "anticodon"))
    expect_equal(translated_col$codon, c("GCC", "AGA", "GCG"))
    expect_equal(translated_col$gene_name, codon_df$gene_name)

    ## CASE 2: Anticodon to codon translation on row names
    anticodon_df <- data.frame(value = 1:2, row.names = c("UUC", "CCG"))
    expect_no_error(translated_ac <- featuresToAA(data = anticodon_df, position = "row", notation_from = "anticodon", notation_to = "codon"))
    expect_true(all(c("GAA", "CGG") %in% rownames(translated_ac)))
})

test_that("featuresToAA handles valid notation conversions across matrix formats", {
    mock_codons <- c("UUU", "GCA", "CGC", "AUG")
    mock_matrix <- matrix(1, nrow = 4, ncol = 2, dimnames = list(mock_codons, c("S1", "S2")))

    valid_pairs <- list(
        c("codon", "aa"),
        c("anticodon", "aa"),
        c("codon", "anticodon"),
        c("anticodon", "codon")
    )

    for (pair in valid_pairs) {
        ## CASE 1: Vector input
        expect_no_error(featuresToAA(data = mock_codons, notation_from = pair[1], notation_to = pair[2]))

        ## CASE 2: Matrix row input
        expect_no_error(featuresToAA(data = mock_matrix, position = "row", notation_from = pair[1], notation_to = pair[2]))

        ## CASE 3: Matrix column input (transposed)
        expect_no_error(featuresToAA(data = t(mock_matrix), position = "column", notation_from = pair[1], notation_to = pair[2]))
    }
})

test_that("featuresToAA validates position parameters and notation constraints", {
    mock_codons <- c("UUU", "GCA", "CGC")
    mock_matrix <- matrix(1, nrow = 3, ncol = 2, dimnames = list(mock_codons, c("S1", "S2")))

    ## CASE 1: Missing required position parameter for tabular data
    expect_error(featuresToAA(data = mock_matrix, notation_from = "codon", notation_to = "aa"))

    ## CASE 2: Identical notation_from and notation_to arguments
    expect_error(featuresToAA(data = mock_codons, notation_from = "codon", notation_to = "codon"))
    expect_error(featuresToAA(data = mock_matrix, position = "row", notation_from = "anticodon", notation_to = "anticodon"))

    ## CASE 3: Unsupported notation argument
    expect_error(featuresToAA(data = mock_matrix, position = "row", notation_from = "aa", notation_to = "codon"))
})
