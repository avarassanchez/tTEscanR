test_that("The different modes to compute the theoretical translation efficiency scores", {

  data(subset_mRNA_data, subset_tRNA_data, metadata)

  tTEscanR_obj <- Create_tTEscanR_Object(counts = list(mRNA = subset_mRNA_data, tRNA = subset_tRNA_data), meta.data = list(ConditionsLabels = metadata, CorrectionFactor = "tissue"))

  # CASE 1: assays not present in the object
  expect_error(Compute_tTE(object = tTEscanR_obj, level = "codon"))
  expect_error(Compute_tTE(object = tTEscanR_obj, level = "aa"))
  expect_error(Compute_tTE(object = tTEscanR_obj, level = "both"))

  tTEscanR_obj <- ComputeCodonUsage(object = tTEscanR_obj, species = "hg38", additional.metrics = FALSE, reduce = 10000)
  tTEscanR_obj <- ComputeAnticodonUsage(object = tTEscanR_obj)

  expect_no_error(Compute_tTE(object = tTEscanR_obj, level = "codon"))
  expect_error(Compute_tTE(object = tTEscanR_obj, level = "aa"))
  expect_error(Compute_tTE(object = tTEscanR_obj, level = "both"))

  tTEscanR_obj <- ComputeAAUsage(object = tTEscanR_obj, level = "both")

  expect_no_error(Compute_tTE(object = tTEscanR_obj, level = "codon"))
  expect_no_error(Compute_tTE(object = tTEscanR_obj, level = "aa"))
  expect_no_error(Compute_tTE(object = tTEscanR_obj, level = "both"))

  # CASE 2: wrong level parameter
  expect_error(Compute_tTE(object = tTEscanR_obj, level = codon)) # level needs to be a string
  expect_error(Compute_tTE(object = tTEscanR_obj, level = "AA")) # case sensitive
  expect_error(Compute_tTE(object = tTEscanR_obj, level = "anticodon")) # level can take codon, aa or both
})
