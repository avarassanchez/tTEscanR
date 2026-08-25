data("default_tTEscanR_mRNA_data", package = "tTEscanR")
data("default_tTEscanR_metadata", package = "tTEscanR")

base_obj <- createObject(counts = default_tTEscanR_mRNA_data, assay = "mRNA")
obj_with_cu <- computeCodonUsage(object = base_obj, species = "hg38", additional_metrics = FALSE)
valid_obj <- updateObject(object = obj_with_cu, meta_data = default_tTEscanR_metadata, params = list("CorrectionFactor" = "tissue"), sample_id = "conditions")

test_that("showPoolContribution validates object requirements and input arguments", {
    ## CASE 1: Missing required codon usage slot
    expect_error(showPoolContribution(object = base_obj, species = "hg38", verbose = FALSE))

    ## CASE 2: Missing metadata and correction factor parameter
    expect_error(showPoolContribution(object = obj_with_cu, species = "hg38", verbose = FALSE))

    ## CASE 3: Invalid correction factor column name
    invalid_meta_obj <- updateObject(object = obj_with_cu, meta_data = default_tTEscanR_metadata, params = list("CorrectionFactor" = "FakeColumn"), sample_id = "conditions")
    expect_error(showPoolContribution(object = invalid_meta_obj, species = "hg38", verbose = FALSE), "correction factor was not found")

    ## CASE 4: Non-intersecting genes in codon frequency matrix
    bad_codon_freq <- matrix(1, nrow = 4, ncol = 2, dimnames = list(c("AAA", "AAC", "AAG", "AAT"), c("FakeGene1", "FakeGene2")))
    expect_error(showPoolContribution(object = valid_obj, codon_freq = bad_codon_freq, verbose = FALSE), "No mRNAs in common found")
})

test_that("showPoolContribution executes successfully and populates object slots", {
    ## CASE 1: Initial execution with top 10 genes (N = 10)
    res_obj <- showPoolContribution(object = valid_obj, species = "hg38", N = 10, verbose = FALSE)

    expect_true(all(c("SizeCorrectedCodonUsage", "SizeCorrected_mRNA") %in% names(res_obj)))

    meta_results <- S4Vectors::metadata(res_obj)[["CodonPoolContribution_Results"]]
    expect_true(all(c("CodonPoolContribution", "PoolContributorTop10Genes", "CorrelationTop10Genes") %in% names(meta_results)))

    ## CASE 2: Re-run with N = 5 to verify dynamic slot update
    res_obj_5 <- showPoolContribution(object = res_obj, species = "hg38", N = 5, verbose = FALSE, overwrite = TRUE)
    expect_true("PoolContributorTop5Genes" %in% names(S4Vectors::metadata(res_obj_5)[["CodonPoolContribution_Results"]]))
})
