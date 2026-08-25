data("default_tTEscanR_mRNA_data", package = "tTEscanR")
data("default_tTEscanR_tRNA_data", package = "tTEscanR")
data("default_tTEscanR_metadata", package = "tTEscanR")

test_that("computeTheoreticalTE executes correctly across different modes", {
    tTEscanR_obj <- suppressWarnings(createObject(counts = list(mRNA = default_tTEscanR_mRNA_data, tRNA = default_tTEscanR_tRNA_data),
                                                  meta_data = default_tTEscanR_metadata, sample_id = "conditions",
                                                  params = list(CorrectionFactor = "tissue")))

    ## CASE 1: Missing all usage assays
    expect_error(computeTheoreticalTE(object = tTEscanR_obj, level = "both"), "not found")

    ## CASE 2: Codon and Anticodon usage computed (AA usage missing)
    tTEscanR_obj <- computeCodonUsage(object = tTEscanR_obj, species = "hg38", additional_metrics = FALSE, reduce = 10000)
    tTEscanR_obj <- computeAnticodonUsage(object = tTEscanR_obj)

    expect_no_error(computeTheoreticalTE(object = tTEscanR_obj, level = "codon", compute_significance = FALSE))
    expect_error(computeTheoreticalTE(object = tTEscanR_obj, level = "aa"), "not found")

    ## CASE 3: All usage assays present
    tTEscanR_obj <- computeAAUsage(object = tTEscanR_obj, level = "both")

    for (lvl in c("codon", "aa", "both")) {
        expect_no_error(computeTheoreticalTE(object = tTEscanR_obj, level = lvl, compute_significance = FALSE))
    }

    ## CASE 4: Verify full metadata slot update
    tTEscanR_obj <- computeTheoreticalTE(object = tTEscanR_obj, level = "both")
    expect_true("tTE_results" %in% names(S4Vectors::metadata(tTEscanR_obj)))
})

test_that("computeTheoreticalTE rejects invalid level parameter inputs", {
    tTEscanR_obj <- suppressWarnings(createObject(
        counts = list(mRNA = default_tTEscanR_mRNA_data, tRNA = default_tTEscanR_tRNA_data),
        meta_data = default_tTEscanR_metadata,
        sample_id = "conditions",
        params = list(CorrectionFactor = "tissue")
    ))

    ## CASE 1: Invalid level arguments: unquoted symbol, incorrect casing, unsupported level string
    expect_error(computeTheoreticalTE(object = tTEscanR_obj, level = codon))
    expect_error(computeTheoreticalTE(object = tTEscanR_obj, level = "AA"))
    expect_error(computeTheoreticalTE(object = tTEscanR_obj, level = "anticodon"))
})
