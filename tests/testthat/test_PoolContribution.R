data(default_tTEscanR_metadata)

test_that("Checking the codon pool contribution", {
    tTEscanR_obj <- createObject(counts = mRNA_data_test, assay = "mRNA")
    expect_error(showPoolContribution(object = tTEscanR_obj, species = "hg38", verbose = FALSE)) # codon usage not present in object

    tTEscanR_obj <- computeCodonUsage(object = tTEscanR_obj, species = "hg38", additional_metrics = FALSE)
    expect_error(showPoolContribution(object = tTEscanR_obj, species = "hg38", verbose = FALSE)) # metadata and corr_factor not present in object

    tTEscanR_obj <- updateObject(object = tTEscanR_obj, meta.data = list(default_tTEscanR_metadata, "tissue"), meta.data.ids = list("ConditionsLabels", "CorrectionFactor"))
    expect_no_error(showPoolContribution(object = tTEscanR_obj, species = "hg38", verbose = FALSE))
})

test_that("Input validation for showPoolContribution catches bad arguments", {
    # CASE 1: wrong corr_factor
    tTEscanR_obj <- createObject(counts = mRNA_data_test, assay = "mRNA")
    tTEscanR_obj <- computeCodonUsage(object = tTEscanR_obj, species = "hg38", additional_metrics = FALSE)
    tTEscanR_obj <- updateObject(object = tTEscanR_obj, meta.data = list(default_tTEscanR_metadata, "FakeColumn"), meta.data.ids = list("ConditionsLabels", "CorrectionFactor"))
    expect_error(showPoolContribution(object = tTEscanR_obj, species = "hg38", verbose = FALSE), "The correction factor was not found in the metadata")

    # CASE 2: no-intersecting genes
    tTEscanR_obj <- updateObject(object = tTEscanR_obj, meta.data = list(default_tTEscanR_metadata, "tissue"), meta.data.ids = list("ConditionsLabels", "CorrectionFactor"), overwrite = TRUE)
    fake_genes <- c("FakeGene1", "FakeGene2")
    fake_codons <- apply(expand.grid(c("A", "C", "G", "T"), c("A", "C", "G", "T"), c("A", "C", "G", "T")), 1, paste, collapse = "")
    bad_codon_freq <- matrix(1, nrow = 64, ncol = 2, dimnames = list(fake_codons, fake_genes))
    expect_error(showPoolContribution(object = tTEscanR_obj, codon_freq = bad_codon_freq, verbose = FALSE), "No mRNAs in common found")
})

test_that("showPoolContribution correctly updates the tTEscanR object slots", {
    tTEscanR_obj <- createObject(counts = mRNA_data_test, assay = "mRNA")
    tTEscanR_obj <- computeCodonUsage(object = tTEscanR_obj, species = "hg38", additional_metrics = FALSE)
    tTEscanR_obj <- updateObject(object = tTEscanR_obj, meta.data = list(default_tTEscanR_metadata, "tissue"), meta.data.ids = list("ConditionsLabels", "CorrectionFactor"))

    tTEscanR_obj <- showPoolContribution(object = tTEscanR_obj, species = "hg38", N = 10, verbose = FALSE)

    expect_true("SizeCorrectedCodonUsage" %in% names(tTEscanR_obj@assays))
    expect_true("SizeCorrected_mRNA" %in% names(tTEscanR_obj@assays))

    expect_true("CodonPoolContribution_Results" %in% names(tTEscanR_obj@meta.data))

    meta_results <- tTEscanR_obj@meta.data$CodonPoolContribution_Results
    expect_true("CodonPoolContribution" %in% names(meta_results))
    expect_true("PoolContributorTop10Genes" %in% names(meta_results))
    expect_true("CorrelationTop10Genes" %in% names(meta_results))

    # Run it again with a different N to ensure the dynamic naming works
    tTEscanR_obj <- showPoolContribution(object = tTEscanR_obj, species = "hg38", N = 5, verbose = FALSE, overwrite = TRUE)
    expect_true("PoolContributorTop5Genes" %in% names(tTEscanR_obj@meta.data$CodonPoolContribution_Results))
})
