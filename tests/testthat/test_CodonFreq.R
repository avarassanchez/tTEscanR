test_that("Check the codon frequency per gene function - GetCodonFreq()", {

  base_args <- list(out_format = "external_gene_name", verbose = FALSE)

  # CASE 1: NULL inputs
  expect_error(do.call(GetCodonFreq, c(list(dataset_name = NULL, genes_file = NULL), base_args)), "specify either")

  # CASE 2: Invalid genes_file type
  expect_error(do.call(GetCodonFreq, c(list(dataset_name = NULL, genes_file = 123), base_args)))
})

test_that("ExtractCodons produces an accurate count matrix", {

  seq1 <- "ATGCGTACG" # 3 codons: ATG, CGT, ACG
  seq2 <- "TTAAGGCCG" # 3 codons: TTA, AGG, CCG
  input_list <- list(s1 = seq1, s2 = seq2)

  res <- ExtractCodons(sequences = input_list, verbose = FALSE)

  expect_true(is.matrix(res))
  expect_type(res, "double")
  expect_equal(ncol(res), 2)
  expect_true(all(grepl("sequence_1|sequence_2", colnames(res)))) # Check if names are inherited from list or auto-generated

  expect_equal(unname(colSums(res)), c(3, 3)) # Each sequence has 9 bases -> 3 codons. The column sums must be 3.

  expect_true("ATG" %in% rownames(res)) # Ensure "ATG" is present and has a count of 1 in the first sequence
  expect_equal(res["ATG", 1], 1)
  expect_equal(res["ATG", 2], 0)

  # RNA and DNA versions of the same sequence should produce identical data frames
  rna_res <- ExtractCodons(sequences = "AUGCUGCAU", verbose = FALSE)
  dna_res <- ExtractCodons(sequences = "ATGCTGCAT", verbose = FALSE)

  expect_equal(rna_res, dna_res)
})

test_that("ExtractCodons handles data.frame inputs correctly", {
  df_input <- data.frame(seqs = c("ATGCGTACG", "TTAAGGCCG"), id = c("seq1", "seq2"))
  res_df <- ExtractCodons(sequences = df_input, verbose = FALSE)

  expect_true(is.matrix(res_df))
  expect_type(res_df, "double")
  expect_equal(sum(res_df), 6) # 6 unique codons total in this example
})

test_that("ExtractCodons correctly parses sequences and handles inputs", {

  expected_counts <- c(ATG = 1, CTG = 1, CAT = 1)
  variants <- c("ATGCTGCAT", "AUGCUGCAU", "atgctgcat", "augcugcau")

  for (seq in variants) {
    res <- ExtractCodons(sequences = seq, verbose = FALSE)
    expect_equal(as.numeric(res[names(expected_counts), 1]), unname(expected_counts), info = paste("Failed on variant:", seq))
  }

  # Dealing with a named list
  seq_list <- list(s1 = "ATGCTGCAT", s2 = "ATGCTGCAT")
  res_list <- ExtractCodons(sequences = seq_list, verbose = FALSE)
  expect_equal(ncol(res_list), 2)
  expect_equal(unname(colSums(res_list)), c(3, 3)) # Each seq has 3 codons

  # Dealing with a list
  seq_list <- list("ATGCTGCAT", "ATGCTGCAT")
  res_list <- ExtractCodons(sequences = seq_list, verbose = FALSE)
  expect_equal(ncol(res_list), 2)
  expect_equal(unname(colSums(res_list)), c(3, 3))

  # Dealing with a vector
  seq_vector <- c("ATGCTGCAT", "ATGCTGCAT")
  res_vector <- ExtractCodons(sequences = seq_vector, verbose = FALSE)
  expect_equal(ncol(res_vector), 2)
  expect_equal(unname(colSums(res_vector)), c(3, 3))

  # Dealing with a dataframe without ids
  seq_df <- data.frame(sequences = c("ATGCTGCAT", "ATGCTGCAT"))
  res_dg <- ExtractCodons(sequences = seq_df, verbose = FALSE)
  expect_equal(ncol(res_dg), 2)
  expect_equal(unname(colSums(res_dg)), c(3, 3))


  expect_error(ExtractCodons(sequences = 12345), "character")
  expect_warning(ExtractCodons(sequences = c("ATGCTGCAT", "ATGCTGCATC")), "multiple of 3")
  expect_error(suppressWarnings(ExtractCodons(sequences = "ATGCTGCATC")), "No valid sequences remained after filtering")
  expect_error(ExtractCodons(sequences = "ARGCRGCAR"), "Invalid sequence characters")
})
