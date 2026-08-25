data("default_tTEscanR_mRNA_data", package = "tTEscanR")
data("default_tTEscanR_tRNA_data", package = "tTEscanR")
data("default_tTEscanR_metadata", package = "tTEscanR")

test_that("The different modes to compute the amino acid usage", {
    tTEobject <- suppressWarnings(createObject(
        counts = list(default_tTEscanR_mRNA_data, default_tTEscanR_tRNA_data), assay = list("mRNA", "tRNA"),
        meta_data = default_tTEscanR_metadata, sample_id = "conditions", params = list("CorrectionFactor" = "tissue")
    ))

    tTEobject <- computeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", reduce = 1000, additional_metrics = FALSE, verbose = FALSE)
    tTEobject <- computeAnticodonUsage(object = tTEobject, verbose = FALSE)

    ## CASE 1: base cases
    expect_no_error(suppressWarnings(computeAAUsage(object = tTEobject, level = "demand", verbose = FALSE)))
    expect_no_error(suppressWarnings(computeAAUsage(object = tTEobject, level = "supply", verbose = FALSE)))
    expect_no_error(suppressWarnings(computeAAUsage(object = tTEobject, level = "both", verbose = FALSE)))

    ## CASE 2: invalid level parameter
    expect_error(suppressWarnings(computeAAUsage(object = tTEobject, level = "invalid_level", verbose = FALSE)))

    ## CASE 3: missing dataset in the input object
    tTEobject <- createObject(
        counts = default_tTEscanR_mRNA_data, assay = "mRNA", meta_data = default_tTEscanR_metadata,
        sample_id = "conditions", params = list("CorrectionFactor" = "tissue")
    )

    ## Misisng CodonUsage
    expect_error(suppressWarnings(computeAAUsage(object = tTEobject, level = "demand", verbose = FALSE)))
    tTEobject <- computeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", reduce = 1000, additional_metrics = FALSE, verbose = FALSE)
    expect_no_error(suppressWarnings(computeAAUsage(object = tTEobject, level = "demand", verbose = FALSE)))

    ## Missing AnticodonUSage
    expect_error(suppressWarnings(computeAAUsage(object = tTEobject, level = "supply", verbose = FALSE)))

    # CASE 4: invalid input object class
    codon_usage <- SummarizedExperiment::assay(tTEobject, "CodonUsage")
    expect_error(suppressWarnings(computeAAUsage(object = codon_usage, level = "demand", verbose = FALSE)), "'object' must be a MultiAssayExperiment object.")
})
