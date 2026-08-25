data("default_tTEscanR_mRNA_data", package = "tTEscanR")

test_that("checkCodonFreqTable handles species lookup and custom tables", {
    ## CASE 1: Valid inputs: default species lookup and user-defined matrix
    expect_no_error(checkCodonFreqTable(data = NULL, species = "hg38"))
    expect_no_error(checkCodonFreqTable(data = codon_freq_table_canonical_hg38, species = NULL))

    ## CASE 2: Invalid species parameter (tests error handling for bad string and NULL)
    expect_error(checkCodonFreqTable(data = NULL, species = "human"))
    expect_error(checkCodonFreqTable(data = NULL, species = NULL))
})

test_that("checkDataFrame validates tabular structures and dimnames", {
    ## CASE 1: Valid matrix/data.frame checks
    expect_no_error(checkDataFrame(data = default_tTEscanR_mRNA_data, required_names = TRUE))
    expect_no_error(checkDataFrame(data = default_tTEscanR_mRNA_data, required_names = FALSE))

    ## CASE 2: Invalid non-tabular input structures
    expect_error(checkDataFrame(data = default_tTEscanR_mRNA_data[, 1]))
    expect_error(checkDataFrame(data = list("A", "B", "C")))

    ## CASE 3: Missing dimnames
    unnamed_df <- default_tTEscanR_mRNA_data
    colnames(unnamed_df) <- NULL
    expect_error(checkDataFrame(data = unnamed_df, required_names = TRUE))
})
