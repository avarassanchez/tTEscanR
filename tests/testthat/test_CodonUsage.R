data("default_tTEscanR_mRNA_data", package = "tTEscanR")
data("default_tTEscanR_tRNA_data", package = "tTEscanR")
data("default_tTEscanR_metadata", package = "tTEscanR")

test_that("The different corr_method possibilities - inside computeCodonUsage()", {
    tTEobject <- createObject(counts = default_tTEscanR_mRNA_data, assay = "mRNA", meta_data = default_tTEscanR_metadata, sample_id = "conditions", params = list("CorrectionFactor" = "tissue"))

    ## CASE 1: valid parameter combinations
    expect_no_error(suppressWarnings(computeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", additional_metrics = FALSE, verbose = FALSE)))
    expect_no_error(suppressWarnings(computeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", additional_metrics = TRUE, corr_method = "spearman", verbose = FALSE)))
    expect_no_error(suppressWarnings(computeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", additional_metrics = TRUE, corr_method = "pearson", verbose = FALSE)))

    ## CASE 2:Invalid corr_method parameter (1 assertion tests argument validation)
    expect_error(suppressWarnings(computeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", additional_metrics = TRUE, corr_method = "invalid_method", verbose = FALSE)))

    expect_no_error(suppressWarnings(computeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", additional_metrics = FALSE, verbose = FALSE)))
    expect_no_error(suppressWarnings(computeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", additional_metrics = FALSE, corr_method = "spearman", verbose = FALSE))) # corr_method will be ignored as it is only required if additional_metrics = TRUE
})

test_that("The different corr_method possibilities - inside computeMetricsCodonUsage()", {
    tTEobject <- createObject(counts = default_tTEscanR_mRNA_data, assay = "mRNA", meta_data = default_tTEscanR_metadata, sample_id = "conditions")
    tTEobject <- computeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", additional_metrics = FALSE)

    meta <- SummarizedExperiment::colData(tTEobject)
    codon_usage <- SummarizedExperiment::assay(tTEobject, "CodonUsage")

    ## CASE 1: Valid methods
    expect_no_error(suppressWarnings(computeMetricsCodonUsage(codon_usage = codon_usage, codon_freq = codon_freq_table_canonical_hg38, metadata = meta, corr_method = "spearman", batch = "tissue", verbose = FALSE)))
    expect_no_error(suppressWarnings(computeMetricsCodonUsage(codon_usage = codon_usage, codon_freq = codon_freq_table_canonical_hg38, metadata = meta, corr_method = "kendall", batch = "tissue", verbose = FALSE)))

    ## CASE 2: Invalid method and missing required batch parameter
    expect_error(suppressWarnings(computeMetricsCodonUsage(codon_usage = codon_usage, codon_freq = codon_freq_table_canonical_hg38, metadata = meta, corr_method = "invalid_method", verbose = FALSE)))
    expect_error(suppressWarnings(computeMetricsCodonUsage(codon_usage = codon_usage, codon_freq = codon_freq_table_canonical_hg38, metadata = meta, corr_method = "spearman", verbose = FALSE)))
})

test_that("The different corr_method possibilities - inside computeCorrelationBackground()", {
    tTEobject <- createObject(counts = default_tTEscanR_mRNA_data, assay = "mRNA")
    tTEobject <- computeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", additional_metrics = FALSE)

    codon_usage <- SummarizedExperiment::assay(tTEobject, "CodonUsage")
    mean_codon_usage <- suppressWarnings(computeMeanUsage(data = codon_usage, metadata = default_tTEscanR_metadata, batch = "tissue", mode = "raw"))
    exonic_background <- suppressWarnings(computeExonicBackground(data = codon_usage))

    ## CASE 1: Valid methods
    expect_no_error(suppressWarnings(computeCorrelationBackground(mean = mean_codon_usage, background = exonic_background, corr_method = "spearman")))
    expect_no_error(suppressWarnings(computeCorrelationBackground(mean = mean_codon_usage, background = exonic_background, corr_method = "pearson")))

    ## CASE 2: Invalid method
    expect_error(suppressWarnings(computeCorrelationBackground(mean = mean_codon_usage, background = exonic_background, corr_method = "invalid_method")))
})

test_that("The different parameters to perform the mean codon usage - computeMeanUsage()", {
    tTEobject <- createObject(counts = default_tTEscanR_mRNA_data, assay = "mRNA", meta_data = default_tTEscanR_metadata, sample_id = "conditions", params = list("CorrectionFactor" = "tissue"))
    tTEobject <- computeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", additional_metrics = FALSE)

    meta <- SummarizedExperiment::colData(tTEobject)
    codon_usage <- SummarizedExperiment::assay(tTEobject, "CodonUsage")

    # CASE 1: Valid object and inputs
    expect_no_error(suppressWarnings(computeMeanUsage(data = tTEobject, assay = "CodonUsage", mode = "raw")))
    expect_no_error(suppressWarnings(computeMeanUsage(data = tTEobject, assay = "CodonUsage", mode = "raw", metadata = meta, batch = "tissue"))) # metadata and batch will be ignored

    # CASE 2: Invalid object/assay specifications
    expect_error(computeMeanUsage(data = tTEobject, assay = CodonUsage, mode = "raw", metadata = meta, batch = "tissue"))
    expect_error(computeMeanUsage(data = tTEobject, assay = "NonExistentAssay", mode = "raw", metadata = meta))
    expect_error(computeMeanUsage(data = tTEobject, mode = "raw"))

    # CASE 3: missing specifying the batch or metadata
    expect_error(suppressWarnings(computeMeanUsage(data = codon_usage, mode = "raw", metadata = meta)))
    expect_error(suppressWarnings(computeMeanUsage(data = codon_usage, mode = "raw", batch = "tissue")))
})

test_that("The mean usage across different features - computeMeanUsage()", {
    tTEobject <- suppressWarnings(createObject(counts = list(default_tTEscanR_mRNA_data, default_tTEscanR_tRNA_data),
                                               assay = list("mRNA", "tRNA"), meta_data = default_tTEscanR_metadata,
                                               sample_id = "conditions", params = list("CorrectionFactor" = "tissue")))
    tTEobject <- computeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", additional_metrics = FALSE)
    tTEobject <- computeAnticodonUsage(object = tTEobject)
    tTEobject <- computeAAUsage(object = tTEobject, level = "both")

    meta <- SummarizedExperiment::colData(tTEobject)
    features <- c("AnticodonUsage", "AADemand", "AASupply")

    for (feature in features) {
        expect_no_error(suppressWarnings(computeMeanUsage(data = tTEobject, assay = feature, mode = "raw")))
        expect_no_error(suppressWarnings(computeMeanUsage(data = SummarizedExperiment::assay(tTEobject, feature), mode = "raw", metadata = meta, batch = "tissue")))
    }
})
