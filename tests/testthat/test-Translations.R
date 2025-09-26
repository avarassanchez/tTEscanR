test_that("FeaturesToAA stops when no overlapping features are found", {
  non_matching_codons <- c("XXX", "YYY", "ZZZ")
  expect_error(FeaturesToAA(data_to_translate = non_matching_codons, notation.from = "codon", notation.to = "aa"), "No overlapping features have been identified belonging to codon")
})

test_that("FeaturesToAA correctly handles a data.frame with non-matching row names", {
  non_matching_df <- data.frame("value" = 1:5, row.names = c("XYZ", "ABC", "DEF", "GHI", "JKL"))
  expect_error(FeaturesToAA(data_to_translate = non_matching_df, position = "row", notation.from = "codon", notation.to = "aa"), "No overlapping features have been identified belonging to codon")
})

test_that("FeaturesToAA correctly handles a data.frame with some matching row names", {
  partial_match_df <- data.frame("value" = 1:5, row.names = c("UUU", "GCA", "CGC", "XYZ", "ABC"))
  expect_no_error(translated_df <- FeaturesToAA(data_to_translate = partial_match_df, position = "row", notation.from = "codon", notation.to = "aa"))

  expect_equal(nrow(translated_df), nrow(partial_match_df)) # Check that the number of rows is the same

  # Check if the row names have been translated for the matching features
  expect_true("Phe" %in% rownames(translated_df))
  expect_true("Ala" %in% rownames(translated_df))
  expect_true("Arg" %in% rownames(translated_df))
})

test_that("FeaturesToAA correctly handles a vector with some matching codons", {
  partial_match_vector <- c("TCT", "CGC", "GGC", "XYZ", "ABC")
  expect_no_error(translated_vector <- FeaturesToAA( data_to_translate = partial_match_vector, notation.from = "codon", notation.to = "aa"))

  expect_equal(length(translated_vector), length(partial_match_vector)) # Check that the number of elements is the same

  # Check if the codons were translated correctly
  expect_true("Ser" %in% translated_vector)
  expect_true("Arg" %in% translated_vector)
  expect_true("Gly" %in% translated_vector)

  expect_true("XYZ" %in% translated_vector) # Check that the non-matching values were preserved
})

test_that("FeaturesToAA translates a data.frame column from codon to anticodon", {
  codon_df <- data.frame("gene_name" = paste0("gene_", 1:3), "codon" = c("GGC", "TCT", "CGC"), "value" = 1:3)
  expect_no_error(translated_df <- FeaturesToAA(data_to_translate = codon_df, position = "codon", notation.from = "codon", notation.to = "anticodon"))

  expect_equal(nrow(translated_df), nrow(codon_df))
  expect_equal(translated_df$codon, c("GCC", "AGA", "GCG")) # Check that the 'codon' column has been correctly translated to anticodons
  expect_equal(translated_df$gene_name, codon_df$gene_name) # Ensure other columns remain unchanged
})

test_that("FeaturesToAA translates a data.frame from anticodon to codon", {
  anticodon_df <- data.frame("value" = 1:2, row.names = c("UUC", "CCG"))
  expect_no_error(translated_df <- FeaturesToAA(data_to_translate = anticodon_df, position = "row", notation.from = "anticodon", notation.to = "codon"))

  expect_true("GAA" %in% rownames(translated_df))
  expect_true("CGG" %in% rownames(translated_df))
})

test_that("The features translation works - FeaturesToAA()", {

  data(subset_mRNA_data)
  tTEobject <- Create_tTEscanR_Object(counts = subset_mRNA_data, assay = "mRNA")
  tTEobject <- ComputeCodonUsage(object = tTEobject, species = "hg38", filter = "canonical", additional.metrics = FALSE)
  codon_usage <- tTEobject@assays$CodonUsage

  # CASE 1: no error - base case using a vector of genes
  expect_no_error(FeaturesToAA(data_to_translate = rownames(codon_usage), notation.from = "codon", notation.to = "aa"))
  expect_no_error(FeaturesToAA(data_to_translate = rownames(codon_usage), notation.from = "anticodon", notation.to = "aa"))
  expect_no_error(FeaturesToAA(data_to_translate = rownames(codon_usage), notation.from = "codon", notation.to = "anticodon"))
  expect_no_error(FeaturesToAA(data_to_translate = rownames(codon_usage), notation.from = "anticodon", notation.to = "codon"))

  expect_no_error(FeaturesToAA(data_to_translate = colnames(t(codon_usage)), notation.from = "codon", notation.to = "aa"))
  expect_no_error(FeaturesToAA(data_to_translate = colnames(t(codon_usage)), notation.from = "anticodon", notation.to = "aa"))
  expect_no_error(FeaturesToAA(data_to_translate = colnames(t(codon_usage)), notation.from = "codon", notation.to = "anticodon"))
  expect_no_error(FeaturesToAA(data_to_translate = colnames(t(codon_usage)), notation.from = "anticodon", notation.to = "codon"))

  # CASE 1.1: no error - base case using a data table
  expect_no_error(FeaturesToAA(data_to_translate = codon_usage, position = "row", notation.from = "codon", notation.to = "aa"))
  expect_no_error(FeaturesToAA(data_to_translate = codon_usage, position = "row", notation.from = "anticodon", notation.to = "aa"))
  expect_no_error(FeaturesToAA(data_to_translate = codon_usage, position = "row", notation.from = "codon", notation.to = "anticodon"))
  expect_no_error(FeaturesToAA(data_to_translate = codon_usage, position = "row", notation.from = "anticodon", notation.to = "codon"))

  expect_no_error(FeaturesToAA(data_to_translate = t(codon_usage), position = "column", notation.from = "codon", notation.to = "aa"))
  expect_no_error(FeaturesToAA(data_to_translate = t(codon_usage), position = "column", notation.from = "anticodon", notation.to = "aa"))
  expect_no_error(FeaturesToAA(data_to_translate = t(codon_usage), position = "column", notation.from = "codon", notation.to = "anticodon"))
  expect_no_error(FeaturesToAA(data_to_translate = t(codon_usage), position = "column", notation.from = "anticodon", notation.to = "codon"))

  # CASE 1.2: error - base case using a data table (missing position parameter)
  expect_error(FeaturesToAA(data_to_translate = codon_usage, notation.from = "codon", notation.to = "aa"))
  expect_error(FeaturesToAA(data_to_translate = codon_usage, notation.from = "anticodon", notation.to = "aa"))
  expect_error(FeaturesToAA(data_to_translate = codon_usage, notation.from = "codon", notation.to = "anticodon"))
  expect_error(FeaturesToAA(data_to_translate = codon_usage, notation.from = "anticodon", notation.to = "codon"))

  expect_error(FeaturesToAA(data_to_translate = t(codon_usage), notation.from = "codon", notation.to = "aa"))
  expect_error(FeaturesToAA(data_to_translate = t(codon_usage), notation.from = "anticodon", notation.to = "aa"))
  expect_error(FeaturesToAA(data_to_translate = t(codon_usage), notation.from = "codon", notation.to = "anticodon"))
  expect_error(FeaturesToAA(data_to_translate = t(codon_usage), notation.from = "anticodon", notation.to = "codon"))

  # CASE 2: error - same input and output format using a vector of genes
  expect_error(FeaturesToAA(data_to_translate = rownames(codon_usage), notation.from = "codon", notation.to = "codon"))
  expect_error(FeaturesToAA(data_to_translate = rownames(codon_usage), notation.from = "anticodon", notation.to = "anticodon"))
  expect_error(FeaturesToAA(data_to_translate = colnames(t(codon_usage)), notation.from = "codon", notation.to = "codon"))
  expect_error(FeaturesToAA(data_to_translate = colnames(t(codon_usage)), notation.from = "anticodon", notation.to = "anticodon"))

  # CASE 2.1: error - same input and output format using a data table
  expect_error(FeaturesToAA(data_to_translate = codon_usage, position = "row", notation.from = "codon", notation.to = "codon"))
  expect_error(FeaturesToAA(data_to_translate = codon_usage, position = "row", notation.from = "anticodon", notation.to = "anticodon"))
  expect_error(FeaturesToAA(data_to_translate = t(codon_usage), position = "column", notation.from = "codon", notation.to = "codon"))
  expect_error(FeaturesToAA(data_to_translate = t(codon_usage), position = "column", notation.from = "anticodon", notation.to = "anticodon"))

  # CASE 3: error - wrong usage of the parameters using a vector of genes
  expect_error(FeaturesToAA(data_to_translate = codon_usage, position = "row", notation.from = "aa", notation.to = "codon")) # wrong notation.to parameter

  # expect_error(FeaturesToAA(data_to_translate = codon_usage, position = "rows", notation.from = "anticodon", notation.to = "codon")) # wrong position parameter
})
