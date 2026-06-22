test_that("The different modes to compute the theoretical translation efficiency scores", {
    data(default_tTEscanR_mRNA_data, default_tTEscanR_tRNA_data, default_tTEscanR_metadata)

    tTEscanR_obj <- createObject(counts = list(mRNA = default_tTEscanR_mRNA_data, tRNA = default_tTEscanR_tRNA_data), meta.data = list(ConditionsLabels = default_tTEscanR_metadata, CorrectionFactor = "tissue"))

    # CASE 1: assays not present in the object
    expect_error(computeTheoreticalTE(object = tTEscanR_obj, level = "codon"), "not found")
    expect_error(computeTheoreticalTE(object = tTEscanR_obj, level = "aa"), "not found")
    expect_error(computeTheoreticalTE(object = tTEscanR_obj, level = "both"), "not found")

    tTEscanR_obj <- computeCodonUsage(object = tTEscanR_obj, species = "hg38", additional_metrics = FALSE, reduce = 10000)
    tTEscanR_obj <- computeAnticodonUsage(object = tTEscanR_obj)

    expect_no_error(computeTheoreticalTE(object = tTEscanR_obj, level = "codon"))
    expect_error(computeTheoreticalTE(object = tTEscanR_obj, level = "aa"))
    expect_error(computeTheoreticalTE(object = tTEscanR_obj, level = "both"))

    tTEscanR_obj <- computeAAUsage(object = tTEscanR_obj, level = "both")

    expect_no_error(computeTheoreticalTE(object = tTEscanR_obj, level = "codon"))
    expect_no_error(computeTheoreticalTE(object = tTEscanR_obj, level = "aa"))

    expect_no_error(computeTheoreticalTE(object = tTEscanR_obj, level = "codon", compute_significance = FALSE))
    expect_no_error(computeTheoreticalTE(object = tTEscanR_obj, level = "aa", compute_significance = FALSE))
    expect_no_error(computeTheoreticalTE(object = tTEscanR_obj, level = "both", compute_significance = FALSE))
    tTEscanR_obj <- computeTheoreticalTE(object = tTEscanR_obj, level = "both")
    expect_true(all(c("tTEresults_codon", "tTEresults_AA") %in% names(tTEscanR_obj@meta.data)))

    # CASE 2: wrong level parameter
    expect_error(computeTheoreticalTE(object = tTEscanR_obj, level = codon)) # level needs to be a string
    expect_error(computeTheoreticalTE(object = tTEscanR_obj, level = "AA")) # case sensitive
    expect_error(computeTheoreticalTE(object = tTEscanR_obj, level = "anticodon")) # level can take codon, aa or both
})
