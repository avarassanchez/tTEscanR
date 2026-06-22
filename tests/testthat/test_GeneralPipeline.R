test_that("runPipeline executes successfully with default parameters", {
    data(default_tTEscanR_mRNA_data, default_tTEscanR_tRNA_data, default_tTEscanR_metadata)

    expect_error(
        res <- runPipeline(
            mRNA_data = default_tTEscanR_mRNA_data,
            tRNA_data = default_tTEscanR_tRNA_data,
            metadata = default_tTEscanR_metadata,
            verbose = TRUE
        ),
        "argument \"batch\" is missing, with no default"
    )

    expect_error(
        res <- runPipeline(
            mRNA_data = default_tTEscanR_mRNA_data,
            tRNA_data = default_tTEscanR_tRNA_data,
            metadata = default_tTEscanR_metadata,
            batch = "tissue",
            verbose = TRUE
        ),
        "No 'codon_freq' provided and no 'species' specified"
    )

    expect_error(
        res <- runPipeline(
            mRNA_data = default_tTEscanR_mRNA_data,
            tRNA_data = default_tTEscanR_tRNA_data,
            metadata = default_tTEscanR_metadata,
            batch = "tissue",
            species = "hg19",
            verbose = TRUE
        ),
        "Incorrect 'species'"
    )

    expect_error(
        res <- runPipeline(
            mRNA_data = default_tTEscanR_mRNA_data,
            tRNA_data = default_tTEscanR_tRNA_data,
            metadata = default_tTEscanR_metadata,
            batch = "tissues",
            species = "hg38",
            verbose = TRUE
        ),
        "The correction factor was not found in the metadata"
    )

    expect_no_error(
        res <- runPipeline(
            mRNA_data = default_tTEscanR_mRNA_data,
            tRNA_data = default_tTEscanR_tRNA_data,
            metadata = default_tTEscanR_metadata,
            batch = "tissue",
            species = "hg38",
            verbose = TRUE
        )
    )
    expect_false(is.null(res))
    expect_s4_class(res, "tTEscanR_Object")

    expect_no_error(
        res <- runPipeline(
            mRNA_data = default_tTEscanR_mRNA_data,
            tRNA_data = default_tTEscanR_tRNA_data,
            metadata = default_tTEscanR_metadata,
            batch = "tissue",
            species = "hg38",
            verbose = FALSE
        )
    )
    expect_false(is.null(res))
    expect_s4_class(res, "tTEscanR_Object")

    expect_no_error(
        res <- runPipeline(
            mRNA_data = default_tTEscanR_mRNA_data,
            tRNA_data = default_tTEscanR_tRNA_data,
            metadata = default_tTEscanR_metadata,
            batch = "tissue",
            species = "hg38",
            runDESeq = FALSE,
            verbose = FALSE
        )
    )
    expect_false(is.null(res))
    expect_s4_class(res, "tTEscanR_Object")

    expect_no_error(
        res <- runPipeline(
            mRNA_data = default_tTEscanR_mRNA_data,
            tRNA_data = default_tTEscanR_tRNA_data,
            metadata = default_tTEscanR_metadata,
            batch = "tissue",
            species = "hg38",
            additional_metrics = FALSE,
            verbose = FALSE
        )
    )
    expect_false(is.null(res))
    expect_s4_class(res, "tTEscanR_Object")

    expect_error(
        res <- runPipeline(
            mRNA_data = default_tTEscanR_mRNA_data,
            tRNA_data = default_tTEscanR_tRNA_data,
            metadata = default_tTEscanR_metadata,
            species = "hg38",
            additional_metrics = FALSE,
            runDESeq = FALSE,
            verbose = FALSE
        )
    ) # batch parameter needed to compute the tTE score
})
