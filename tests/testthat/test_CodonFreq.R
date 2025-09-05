test_that("Check the codon frequency per gene function - ObtainCodonFreqPerGene()", {

  # CASE 1: error - no input specified
  expect_error(ObtainCodonFreqPerGene(dataset_name = NULL, genes_file = NULL, transcripts = NULL, filter = "canonical", out_format = "external_gene_name", verbose = TRUE))
  # CASE 2: error - wrong genes_file
  expect_error(ObtainCodonFreqPerGene(dataset_name = NULL, genes_file = genes, transcripts = NULL, filter = "canonical", out_format = "external_gene_name", verbose = TRUE))
})

test_that("Check the codon composition - ExtractCodonComposition()", {

  # CASE 1: no error - base cases
  expect_no_error(ExtractCodonComposition(sequences = "ATGCTGCAT", verbose = FALSE)) # using T bases
  expect_no_error(ExtractCodonComposition(sequences = "AUGCUGCAU", verbose = FALSE)) # using U bases
  expect_no_error(ExtractCodonComposition(sequences = "atgctgcat", verbose = FALSE)) # using T bases
  expect_no_error(ExtractCodonComposition(sequences = "augcugcau", verbose = FALSE)) # using U bases

  expect_no_error(ExtractCodonComposition(sequences = list("ATGCTGCAT"), verbose = FALSE)) # list of one sequence
  expect_no_error(ExtractCodonComposition(sequences = list("ATGCTGCAT", "ATGCTGCAT"), verbose = FALSE)) # list of multiple sequences

  sequences_example <- data.frame(sequences = c("ATGCTGCAT", "ATGCTGCAT"), ids = c("sequence1", "sequence2"))
  expect_no_error(ExtractCodonComposition(sequences = sequences_example, verbose = FALSE)) # dataframe as input

  # CASE 2: error - wrong sequences input
  expect_error(ExtractCodonComposition(sequences = ATGCTGCAT, verbose = FALSE)) # not a string
  expect_error(ExtractCodonComposition(sequences = "ATGCTGCATC", verbose = FALSE)) # not a multiple of 3 - ERROR
  expect_error(ExtractCodonComposition(sequences = "ARGCRGCAR", verbose = FALSE)) # wrong nucleotide bases
})
