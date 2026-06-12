test_that("The default data is correctly loaded", {
    # CASE 1: no error - base case
    params <- list(
        list(s = "hg38", ref = codon_freq_table_canonical_hg38[rownames(codon_freq_table_canonical_hg38), ]),
        list(s = "mm39", ref = codon_freq_table_canonical_mm39[rownames(codon_freq_table_canonical_mm39), ])
    )

    for (p in params) {
        res <- SelectDefaultData(p$s)
        expect_identical(res, p$ref)
    }

    expect_no_error(SelectDefaultData(species = "mm39"))
    expect_no_error(SelectDefaultData(species = "hg38"))

    # CASE 2: error - wrong or missing parameters
    expect_error(SelectDefaultData(species = "mm10"))
    expect_error(SelectDefaultData(species = "HG38"))
})

test_that("The data loaded has a correct format", {
    # CONTROLS
    expect_no_error(IdentifyInputFormat(data = mRNA_data_test, mode = "flexible"))
    expect_no_error(IdentifyInputFormat(data = mRNA_data_test, mode = "fix"))
    expect_no_error(IdentifyInputFormat(data = ENSG_gene_names_mRNA_data, mode = "flexible"))
    expect_no_error(IdentifyInputFormat(data = ENSG_gene_names_mRNA_data, mode = "fix"))

    # CASE 1: error - wrong input formats
    expect_no_error(IdentifyInputFormat(data = list(4, 5, 6, 7, 8), mode = "flexible"))
    expect_error(IdentifyInputFormat(data = list(4, 5, 6, 7, 8), mode = "fix")) # the fix parameters requires that the elements are matrices

    # CASE 2: error - empty input data
    expect_error(IdentifyInputFormat(data = c(), mode = "flexible"))
    expect_error(IdentifyInputFormat(data = c(), mode = "fix"))

    # CASE 3: missing mode parameter
    expect_no_error(IdentifyInputFormat(data = mRNA_data_test)) # mode fix is selected by default
    expect_no_error(IdentifyInputFormat(data = ENSG_gene_names_mRNA_data)) # mode fix is selected by default
    expect_error(IdentifyInputFormat(data = list(4, 5, 6, 7, 8)))

    # CASE 4: error - wrong mode parameter
    expect_error(IdentifyInputFormat(data = mRNA_data_test, mode = flexible)) # mode needs to be a string
    expect_error(IdentifyInputFormat(data = mRNA_data_test, mode = "strict")) # mode has to be flexible or fix
})

test_that("The individual columns are properly groupped", {
    # CASE 1 - Basic dense matrix
    data <- data.frame(sample_1 = c(1, 2), sample_2 = c(3, 4), sample_3 = c(5, 6))
    rownames(data) <- c("gene1", "gene2")
    groups <- c("A", "A", "B")

    result <- GroupConditions(data, groups)

    expect_equal(colnames(result), c("A", "B"))
    expect_equal(result$A, c(4, 6)) # 1+3, 2+4
    expect_equal(result$B, c(5, 6))

    # CASE 2 - Sparse matrix
    data_sparse <- Matrix::Matrix(c(1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0), nrow = 2, ncol = 6, byrow = TRUE)
    rownames(data_sparse) <- c("gene1", "gene2")
    groups <- c("A", "A", "A", "B", "B", "B")

    result <- GroupConditions(data_sparse, groups)

    expect_equal(colnames(result), c("A", "B"))
    expect_equal(result$A, c(3, 9)) # 1+0+2, 4+0+5
    expect_equal(result$B, c(3, 6)) # 0+3+0, 0+6+0

    # Create a sparse matrix (column-major by default)
    data_sparse <- Matrix::Matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 6)
    rownames(data_sparse) <- c("gene1", "gene2")
    groups <- c("A", "A", "A", "B", "B", "B")

    result <- GroupConditions(data_sparse, groups)

    expected <- sapply(unique(groups), function(g) Matrix::rowSums(data_sparse[, groups == g, drop = FALSE]))
    expected <- as.data.frame(expected)
    colnames(expected) <- unique(groups)
    rownames(expected) <- rownames(data_sparse)
    expect_equal(result, expected)

    # CASE 3 - Single group
    data <- data.frame(sample_1 = c(1, 2), sample_2 = c(3, 4))
    rownames(data) <- c("gene1", "gene2")
    groups <- c("A", "A")

    result <- GroupConditions(data, groups)

    expect_equal(ncol(result), 1)
    expect_equal(result$A, c(4, 6))

    # CASE 4 - Unique groups (no repeats)
    data <- data.frame(sample_1 = c(1, 2), sample_2 = c(3, 4))
    rownames(data) <- c("gene1", "gene2")
    groups <- c("A", "B")

    result <- GroupConditions(data, groups)

    expect_equal(colnames(result), c("A", "B"))
    expect_equal(result$A, c(1, 2))
    expect_equal(result$B, c(3, 4))

    # CASE 5 - Non-numeric data fails gracefully
    data <- data.frame(sample_1 = c("a", "b"), sample_2 = c("c", "d"))
    rownames(data) <- c("gene1", "gene2")
    groups <- c("A", "B")

    result <- tryCatch(GroupConditions(data, groups), error = function(e) e)

    expect_true(inherits(result, "error"))

    # CASE 6 - Dealing with missing values
    data <- data.frame(sample_1 = c(1, NA), sample_2 = c(3, 4))
    rownames(data) <- c("gene1", "gene2")
    groups <- c("A", "A")

    result <- GroupConditions(data, groups)

    expect_equal(result$A, c(4, NA)) # NA propagates

    # CASE 6 - Empty matrix
    data <- data.frame()
    groups <- character(0)

    result <- tryCatch(GroupConditions(data, groups), error = function(e) e)

    expect_true(inherits(result, "error"))

    # CASE 7 - "Large random matrix
    set.seed(42)
    data <- matrix(runif(1000 * 20), nrow = 1000, byrow = TRUE)
    groups <- sample(c("A", "B", "C", "D"), 20, replace = TRUE)

    result <- GroupConditions(data, groups)

    expect_equal(nrow(result), 1000)
    expect_equal(ncol(result), length(unique(groups)))
})
