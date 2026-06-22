test_that("The codon frequency per gene matrix is properly assessed", {
    # CASE 1: no error - loading the default codon_freq tables
    expect_no_error(checkCodonFreqTable(data = NULL, species = "hg38"))
    expect_no_error(checkCodonFreqTable(data = NULL, species = "mm39"))

    # CASE 2: error - incorrect species parameter
    expect_error(checkCodonFreqTable(data = NULL, species = "human"))
    expect_error(checkCodonFreqTable(data = NULL, species = NULL))

    # CASE 4: no error - loading user-defined codon_freq table
    expect_no_error(checkCodonFreqTable(data = codon_freq_table_canonical_hg38, species = NULL))
    expect_no_error(checkCodonFreqTable(data = codon_freq_table_canonical_mm39, species = NULL))
})

test_that("The data frames are properly assessed", {
    data(mRNA_data_test, ENSG_gene_names_mRNA_data)

    # CASE 1: no error -  checking well defined matrices
    expect_no_error(checkDataFrame(data = mRNA_data_test))
    expect_no_error(checkDataFrame(data = mRNA_data_test, required_names = TRUE))
    expect_no_error(checkDataFrame(data = mRNA_data_test, required_names = FALSE))

    # CASE 2: error - the input data is a vector
    expect_error(checkDataFrame(data = mRNA_data_test[1, ]))
    expect_error(checkDataFrame(data = mRNA_data_test[, 1]))
    expect_error(checkDataFrame(data = ENSG_gene_names_mRNA_data))
    expect_error(checkDataFrame(data = list("A", "B", "C")))
    expect_error(checkDataFrame(data = "mRNA_data"))

    # CASE 3: error - incorrect colnames and/or rownames
    rownames(mRNA_data_test) <- NULL
    expect_error(checkDataFrame(data = mRNA_data_test, required_names = TRUE))

    colnames(mRNA_data_test) <- NULL
    expect_error(checkDataFrame(data = mRNA_data_test, required_names = TRUE))
})
