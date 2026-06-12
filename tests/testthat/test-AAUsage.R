test_that("The different modes to compute the amino acid usage", {
    data(default_tTEscanR_mRNA_data, default_tTEscanR_tRNA_data, default_tTEscanR_metadata)

    tTEobject <- CreateObject(
        counts = list(default_tTEscanR_mRNA_data, default_tTEscanR_tRNA_data), assay = list("mRNA", "tRNA"),
        meta.data = list(default_tTEscanR_metadata, "tissue"), meta.data.ids = list("ConditionsLabels", "CorrectionFactor")
    )

    tTEobject <- ComputeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", reduce = 1000, additional_metrics = FALSE, verbose = FALSE)
    tTEobject <- ComputeAnticodonUsage(object = tTEobject, verbose = FALSE)

    # CASE 1: o error - base cases
    expect_no_error(suppressWarnings(ComputeAAUsage(object = tTEobject, level = "demand", verbose = FALSE)))
    expect_no_error(suppressWarnings(ComputeAAUsage(object = tTEobject, level = "supply", verbose = FALSE)))
    expect_no_error(suppressWarnings(ComputeAAUsage(object = tTEobject, level = "both", verbose = FALSE)))

    # CASE 2: error - wrong level parameter
    expect_error(suppressWarnings(ComputeAAUsage(object = tTEobject, verbose = FALSE)))
    expect_error(suppressWarnings(ComputeAAUsage(object = tTEobject, level = "AAdemand", verbose = FALSE)))
    expect_error(suppressWarnings(ComputeAAUsage(object = tTEobject, level = "AAsupply", verbose = FALSE)))
    expect_error(suppressWarnings(ComputeAAUsage(object = tTEobject, level = "demand and supply", verbose = FALSE)))

    # CASE 3: error - missing dataset in the input object
    tTEobject <- CreateObject(
        counts = default_tTEscanR_mRNA_data, assay = "mRNA", meta.data = list(default_tTEscanR_metadata, "tissue"),
        meta.data.ids = list("ConditionsLabels", "CorrectionFactor")
    )

    expect_error(suppressWarnings(ComputeAAUsage(object = tTEobject, level = "demand", verbose = FALSE)))
    expect_error(suppressWarnings(ComputeAAUsage(object = tTEobject, level = "supply", verbose = FALSE)))
    expect_error(suppressWarnings(ComputeAAUsage(object = tTEobject, level = "both", verbose = FALSE)))

    tTEobject <- ComputeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", reduce = 1000, additional_metrics = FALSE, verbose = FALSE)

    expect_no_error(suppressWarnings(ComputeAAUsage(object = tTEobject, level = "demand", verbose = FALSE)))
    expect_error(suppressWarnings(ComputeAAUsage(object = tTEobject, level = "supply", verbose = FALSE)))
    expect_error(suppressWarnings(ComputeAAUsage(object = tTEobject, level = "both", verbose = FALSE)))

    # CASE 4: error - using a dataset not a tTEscanR object
    expect_error(suppressWarnings(ComputeAAUsage(object = tTEobject@assays$CodonUsage, level = "demand", verbose = FALSE)))
    expect_error(suppressWarnings(ComputeAAUsage(object = tTEobject@assays$AnticodonUsage, level = "supply", verbose = FALSE)))
    expect_error(suppressWarnings(ComputeAAUsage(object = list(tTEobject@assays$CodonUsage, tTEobject@assays$AnticodonUsage), level = "both", verbose = FALSE)))
})
