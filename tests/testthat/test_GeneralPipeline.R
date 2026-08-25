data("default_tTEscanR_mRNA_data", package = "tTEscanR")
data("default_tTEscanR_tRNA_data", package = "tTEscanR")
data("default_tTEscanR_metadata", package = "tTEscanR")

test_that("runPipeline validates required arguments before execution", {
    ## Baseline pipeline arguments
    base_args <- list(
        mRNA_data = default_tTEscanR_mRNA_data,
        tRNA_data = default_tTEscanR_tRNA_data,
        meta_data = default_tTEscanR_metadata,
        sample_id = "conditions",
        verbose = FALSE
    )

    ## CASE 1: Missing required batch parameter
    expect_error(do.call(runPipeline, base_args), "argument \"batch\" is missing")

    ## CASE 2: Missing species specification
    expect_error(do.call(runPipeline, c(base_args, list(batch = "tissue"))), "No 'codon_freq' provided and no 'species' specified")

    ## CASE 3: Invalid species value
    expect_error(do.call(runPipeline, c(base_args, list(batch = "tissue", species = "hg19"))), "Incorrect 'species'")

    ## CASE 4: Invalid batch column name (not in metadata)
    expect_error(do.call(runPipeline, c(base_args, list(batch = "tissues", species = "hg38"))), "correction factor was not found")
 })

test_that("runPipeline executes successfully across full and lightweight configurations", {
    ## CASE 1: Full execution with default DESeq and additional metrics
    res_full <- suppressWarnings(runPipeline(
        mRNA_data = default_tTEscanR_mRNA_data,
        tRNA_data = default_tTEscanR_tRNA_data,
        meta_data = default_tTEscanR_metadata,
        batch = "tissue",
        species = "hg38",
        sample_id = "conditions",
        compute_pairwise = FALSE,
        verbose = TRUE
    ))

    expect_s4_class(res_full, "MultiAssayExperiment")

    ## CASE 2: Fast execution with optional modules toggled off (runDESeq = FALSE, additional_metrics = FALSE)
    res_fast <- suppressWarnings(runPipeline(
        mRNA_data = default_tTEscanR_mRNA_data,
        tRNA_data = default_tTEscanR_tRNA_data,
        meta_data = default_tTEscanR_metadata,
        batch = "tissue",
        species = "hg38",
        sample_id = "conditions",
        runDESeq = FALSE,
        additional_metrics = FALSE,
        compute_pairwise = FALSE,
        verbose = FALSE
    ))

    expect_s4_class(res_fast, "MultiAssayExperiment")
})
