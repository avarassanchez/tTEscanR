test_that("The features translation works - FeaturesToAA()", {

  data(subset_mRNA_data)
  tTEobject <- Create_tTEscanR_Object(counts = subset_mRNA_data, assay = "mRNA")
  tTEobject <- ComputeCodonUsage(object = tTEobject, species = "hg38", filter = "canonical", additional.metrics = FALSE)
  codon_usage <- tTEobject@assays$CodonUsage

  # CASE 1
  # no error - base case using a vector of genes
  expect_no_error(FeaturesToAA(data_to_translate = rownames(codon_usage), notation.from = "codon", notation.to = "aa"))
  expect_no_error(FeaturesToAA(data_to_translate = rownames(codon_usage), notation.from = "anticodon", notation.to = "aa"))
  expect_no_error(FeaturesToAA(data_to_translate = rownames(codon_usage), notation.from = "codon", notation.to = "anticodon"))
  expect_no_error(FeaturesToAA(data_to_translate = rownames(codon_usage), notation.from = "anticodon", notation.to = "codon"))

  expect_no_error(FeaturesToAA(data_to_translate = colnames(t(codon_usage)), notation.from = "codon", notation.to = "aa"))
  expect_no_error(FeaturesToAA(data_to_translate = colnames(t(codon_usage)), notation.from = "anticodon", notation.to = "aa"))
  expect_no_error(FeaturesToAA(data_to_translate = colnames(t(codon_usage)), notation.from = "codon", notation.to = "anticodon"))
  expect_no_error(FeaturesToAA(data_to_translate = colnames(t(codon_usage)), notation.from = "anticodon", notation.to = "codon"))

  # CASE 1.1
  # no error - base case using a data table
  expect_no_error(FeaturesToAA(data_to_translate = codon_usage, position = "row", notation.from = "codon", notation.to = "aa"))
  expect_no_error(FeaturesToAA(data_to_translate = codon_usage, position = "row", notation.from = "anticodon", notation.to = "aa"))
  expect_no_error(FeaturesToAA(data_to_translate = codon_usage, position = "row", notation.from = "codon", notation.to = "anticodon"))
  expect_no_error(FeaturesToAA(data_to_translate = codon_usage, position = "row", notation.from = "anticodon", notation.to = "codon"))

  expect_no_error(FeaturesToAA(data_to_translate = t(codon_usage), position = "column", notation.from = "codon", notation.to = "aa"))
  expect_no_error(FeaturesToAA(data_to_translate = t(codon_usage), position = "column", notation.from = "anticodon", notation.to = "aa"))
  expect_no_error(FeaturesToAA(data_to_translate = t(codon_usage), position = "column", notation.from = "codon", notation.to = "anticodon"))
  expect_no_error(FeaturesToAA(data_to_translate = t(codon_usage), position = "column", notation.from = "anticodon", notation.to = "codon"))

  # CASE 1.2
  # error - base case using a data table (missing position parameter)
  expect_error(FeaturesToAA(data_to_translate = codon_usage, notation.from = "codon", notation.to = "aa"))
  expect_error(FeaturesToAA(data_to_translate = codon_usage, notation.from = "anticodon", notation.to = "aa"))
  expect_error(FeaturesToAA(data_to_translate = codon_usage, notation.from = "codon", notation.to = "anticodon"))
  expect_error(FeaturesToAA(data_to_translate = codon_usage, notation.from = "anticodon", notation.to = "codon"))

  expect_error(FeaturesToAA(data_to_translate = t(codon_usage), notation.from = "codon", notation.to = "aa"))
  expect_error(FeaturesToAA(data_to_translate = t(codon_usage), notation.from = "anticodon", notation.to = "aa"))
  expect_error(FeaturesToAA(data_to_translate = t(codon_usage), notation.from = "codon", notation.to = "anticodon"))
  expect_error(FeaturesToAA(data_to_translate = t(codon_usage), notation.from = "anticodon", notation.to = "codon"))

  # CASE 2
  # error - same input and output format using a vector of genes
  expect_error(FeaturesToAA(data_to_translate = rownames(codon_usage), notation.from = "codon", notation.to = "codon"))
  expect_error(FeaturesToAA(data_to_translate = rownames(codon_usage), notation.from = "anticodon", notation.to = "anticodon"))
  expect_error(FeaturesToAA(data_to_translate = colnames(t(codon_usage)), notation.from = "codon", notation.to = "codon"))
  expect_error(FeaturesToAA(data_to_translate = colnames(t(codon_usage)), notation.from = "anticodon", notation.to = "anticodon"))

  # CASE 2.1
  # error - same input and output format using a data table
  expect_error(FeaturesToAA(data_to_translate = codon_usage, position = "row", notation.from = "codon", notation.to = "codon"))
  expect_error(FeaturesToAA(data_to_translate = codon_usage, position = "row", notation.from = "anticodon", notation.to = "anticodon"))
  expect_error(FeaturesToAA(data_to_translate = t(codon_usage), position = "column", notation.from = "codon", notation.to = "codon"))
  expect_error(FeaturesToAA(data_to_translate = t(codon_usage), position = "column", notation.from = "anticodon", notation.to = "anticodon"))

  # CASE 3
  # error - wrong usage of the parameters using a vector of genes
  expect_error(FeaturesToAA(data_to_translate = codon_usage, position = "row", notation.from = "aa", notation.to = "codon")) # wrong notation.to parameter
  expect_error(FeaturesToAA(data_to_translate = codon_usage, position = "rows", notation.from = "anticodon", notation.to = "codon")) # wrong position parameter

})
