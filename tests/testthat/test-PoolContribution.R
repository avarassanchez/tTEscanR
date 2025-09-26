test_that("Checking the codon pool contribution", {

  data(subset_mRNA_data, metadata)

  tTEscanR_obj <- Create_tTEscanR_Object(counts = subset_mRNA_data, assay = "mRNA")
  expect_error(ExaminePoolContribution(object = tTEscanR_obj, species = "hg38", verbose = FALSE)) # codon usage not present in object

  tTEscanR_obj <- ComputeCodonUsage(object = tTEscanR_obj, species = "hg38", filter = "canonical", additional.metrics = FALSE)
  expect_error(ExaminePoolContribution(object = tTEscanR_obj, species = "hg38", verbose = FALSE)) # metadata and corr_factor not present in object

  tTEscanR_obj <- Update_tTEscanR_Object(object = tTEscanR_obj, meta.data = list(metadata, "tissue"), meta.data.ids = list("ConditionsLabels", "CorrectionFactor"))
  expect_no_error(ExaminePoolContribution(object = tTEscanR_obj, species = "hg38", verbose = FALSE))
})
