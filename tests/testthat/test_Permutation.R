mock_codon_table <- matrix(
    rpois(40, lambda = 50), nrow = 4, ncol = 10,
    dimnames = list(c("AUU", "AUC", "AUA", "AUG"), paste0("Gene", 1:10))
)

mock_target_data <- data.frame(
    Expression = c(12.5, 4.2),
    row.names = c("Gene1", "Gene2")
)

set.seed(42)
mock_permut_data <- data.frame(
    codon = rep(c("AUG", "UAA"), each = 100),
    freq = c(rnorm(100, mean = 0.25, sd = 0.03), rnorm(100, mean = 0.10, sd = 0.01)),
    stringsAsFactors = FALSE
)

mock_sig_data <- data.frame(
    codon = c("AUG", "UAA"),
    freq = c(0.34, 0.11),
    group = c("Treatment_A", "Treatment_A"),
    sig_adj = c(0.002, 0.045),
    p_val_adj = c(0.002, 0.045), # Added to prevent the mutate evaluation crash
    stringsAsFactors = FALSE
)

test_value_input <- mock_sig_data[, c("group", "codon", "freq")]
test_dist_input <- mock_permut_data[, c("codon", "freq")]

test_that("getPermutationDist computes empirical distributions under valid setups", {
    ## CASE 1: Control pathway (NULL target_data)
    res_ctrl <- getPermutationDist(
        n_permut = 5, n_features = 3, codon_freq = mock_codon_table,
        target_data = NULL, verbose = TRUE
    )
    expect_s3_class(res_ctrl, "data.frame")
    expect_equal(colnames(res_ctrl), c("codon", "freq"))
    expect_equal(nrow(res_ctrl), 20)
    expect_true(all(res_ctrl$freq[res_ctrl$codon == "AUU"] >= 0 & res_ctrl$freq[res_ctrl$codon == "AUU"] <= 1))

    ## CASE 2: Target-isolated feature sampling
    res_target <- getPermutationDist(
        n_permut = 2, target_data = mock_target_data, codon_freq = mock_codon_table,
        verbose = FALSE
    )
    expect_s3_class(res_target, "data.frame")
    expect_false(any(is.na(res_target$freq)))
})

test_that("getPermutationDist validates permutation counts and baseline gene matching", {
    ## CASE1 1: Invalid permutation count
    expect_error(
        getPermutationDist(n_permut = 1, codon_freq = mock_codon_table, verbose = FALSE),
        "must be at least 2"
    )

    ## CASE 2: Missing target genes in baseline profile
    corrupt_target <- data.frame(Expression = c(10, 20), row.names = c("Inexistent_Gene_X", "Inexistent_Gene_Y"))
    expect_error(
        getPermutationDist(n_permut = 10, target_data = corrupt_target, codon_freq = mock_codon_table, verbose = FALSE),
        "None of the genes in 'target_data' are present"
    )
})

test_that("obtainSignificance computes empirical p-values and outputs progress messages", {
    ## CASE 1: Compute significance and capture progress output in one call
    logs <- capture.output(
        res <- obtainSignificance(dist = test_dist_input, value = test_value_input, padj_threshold = 0.05, verbose = TRUE),
        type = "message"
    )

    expect_s3_class(res, "data.frame")
    expect_true(all(c("p_val", "tail", "p_val_adj", "sig_adj") %in% colnames(res)))
    expect_equal(res$tail[res$codon == "AUG"], "right")

    ## CASE 2: Sequential progress logs
    expect_true(any(grepl("--- Computing the statistical significance ---", logs, fixed = TRUE)))
    expect_true(any(grepl("1 . Checking the input data.", logs)))
    expect_true(any(grepl("2 . Calculating empirical p-values.", logs)))
    expect_true(any(grepl("3 . Performing FDR correction.", logs)))

    ## CASE 3: Unmatched observed codon assigns NA
    missing_codon_value <- rbind(test_value_input, data.frame(group = "Treatment_A", codon = "XYZ", freq = 0.50))
    res_missing <- obtainSignificance(dist = test_dist_input, value = missing_codon_value, verbose = FALSE)
    expect_true(is.na(res_missing$p_val[res_missing$codon == "XYZ"]))
    expect_true(is.na(res_missing$tail[res_missing$codon == "XYZ"]))
})

test_that("obtainSignificance validates schema structure and data types", {
    ## CASE 1: Non-conforming schema (4 or 5 columns instead of expected 3)
    corrupt_schema <- test_value_input
    corrupt_schema$ExtraColumn <- "Error"
    expect_error(obtainSignificance(dist = test_dist_input, value = corrupt_schema, verbose = FALSE), "Incorrect number of columns")
    expect_error(obtainSignificance(dist = mock_permut_data, value = mock_sig_data, verbose = FALSE), "Incorrect number of columns")

    ## CASE 2: Non-numeric tracking data
    corrupt_numeric <- test_value_input
    corrupt_numeric$freq <- c("High", "Low")
    expect_error(obtainSignificance(dist = test_dist_input, value = corrupt_numeric, verbose = FALSE), "must be numeric")
})
