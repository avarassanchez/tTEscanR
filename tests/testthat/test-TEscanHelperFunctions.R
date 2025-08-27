test_that("The default data is correctly loaded", {

  # CASE 1
  # no error - base case
  expect_equal(SelectDefaultData(species = "hg38", filter = "length"), codon_freq_table_hg38_length)
  expect_equal(SelectDefaultData(species = "hg38", filter = "canonical"), codon_freq_table_hg38_canonical)
  expect_equal(SelectDefaultData(species = "mm39", filter = "length"), codon_freq_table_mm39_length)
  expect_equal(SelectDefaultData(species = "mm39", filter = "canonical"), codon_freq_table_mm39_canonical)

  # CASE 2
  # error - wrong or missing parameters
  expect_error(SelectDefaultData(species = "mm39"))
  expect_error(SelectDefaultData(species = "hg38"))
  expect_error(SelectDefaultData(species = "mm39", filter = "larger"))

})

test_that("The data loaded has a correct format", {

  expect_error(IdentifyInputFormat(data = list(4, 5, 6, 7, 8)))
  expect_error(IdentifyInputFormat(data = c()))
  expect_no_error(IdentifyInputFormat(data = mRNA_data_test))
  expect_no_error(IdentifyInputFormat(data = ENSG_gene_names_mRNA_data))

})

