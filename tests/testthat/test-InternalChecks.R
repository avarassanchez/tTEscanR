test_that("The codon frequency per gene matrix is properly assessed", {

  # CASE 1: no error - loading the default codon_freq tables
  expect_no_error(CheckCodonFreqTable(data = NULL, species = "hg38"))
  expect_no_error(CheckCodonFreqTable(data = NULL, species = "mm39"))

  # CASE 2: error - incorrect species parameter
  expect_error(CheckCodonFreqTable(data = NULL, species = "human"))
  expect_error(CheckCodonFreqTable(data = NULL, species = NULL))

  # CASE 4: no error - loading user-defined codon_freq table
  expect_no_error(CheckCodonFreqTable(data = codon_freq_table_canonical_hg38, species = NULL))
  expect_no_error(CheckCodonFreqTable(data = codon_freq_table_canonical_mm39, species = NULL))
})

test_that("The data frames are properly assessed", {

  # CASE 1: no error -  checking well defined matrices
  expect_no_error(CheckDataFrame(data = mRNA_data_test))
  expect_no_error(CheckDataFrame(data = mRNA_data_test, required_names = TRUE))
  expect_no_error(CheckDataFrame(data = mRNA_data_test, required_names = FALSE))

  # CASE 2: error - the input data is a vector
  expect_error(CheckDataFrame(data = mRNA_data_test[1, ]))
  expect_error(CheckDataFrame(data = mRNA_data_test[, 1]))
  expect_error(CheckDataFrame(data = ENSG_gene_names_mRNA_data))
  expect_error(CheckDataFrame(data = list("A", "B", "C")))
  expect_error(CheckDataFrame(data = "mRNA_data"))

  # CASE 3: error - incorrect colnames and/or rownames
  example_data <- data(mRNA_data_test)
  rownames(example_data) <- NULL
  expect_error(CheckDataFrame(data = example_data, required_names = TRUE))

  colnames(example_data) <- NULL
  expect_error(CheckDataFrame(data = example_data, required_names = TRUE))
})
