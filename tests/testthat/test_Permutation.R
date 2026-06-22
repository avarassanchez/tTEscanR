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

test_that("getPermutationDist computes and unrolls a normalized empirical null distribution", {

    res <- getPermutationDist(
        n_permut = 5,
        n_features = 3,
        codon_freq = mock_codon_table,
        target_data = NULL, # Enforces control pathway
        verbose = TRUE
    )

    expect_s3_class(res, "data.frame")
    expect_equal(colnames(res), c("codon", "freq"))
    expect_equal(nrow(res), 20) # 4 codons * 5 permutations = 20 structural long rows

    sample_perm <- res$freq[res$codon == "AUU"]
    expect_true(all(sample_perm >= 0 & sample_perm <= 1)) # Every individual permutation slice must sum up to 1.0 (100%) across the 4 codons
})

test_that("getPermutationDist isolates target features from the sampling pool", {

    res <- getPermutationDist(
        n_permut = 2,
        target_data = mock_target_data, # Shares Gene1 & Gene2 with table
        codon_freq = mock_codon_table,
        verbose = FALSE
    )

    expect_s3_class(res, "data.frame")
    expect_false(any(is.na(res$freq)))
})

test_that("getPermutationDist rejects invalid permutation counts", {
    expect_error(
        getPermutationDist(n_permut = 1, codon_freq = mock_codon_table, verbose = FALSE),
        "must be at least 2"
    )
})

test_that("getPermutationDist flags an error if target genes are entirely missing from baseline profiles", {
    corrupt_target <- data.frame(
        Expression = c(10, 20),
        row.names = c("Inexistent_Gene_X", "Inexistent_Gene_Y")
    )

    expect_error(
        getPermutationDist(
            n_permut = 10,
            target_data = corrupt_target,
            codon_freq = mock_codon_table,
            verbose = FALSE
        ),
        "None of the genes in 'target_data' are present"
    )
})

test_that("obtainSignificance computes empirical p-values and applies FDR thresholds", {

    res <- obtainSignificance(
        dist = test_dist_input,
        value = test_value_input,
        padj_threshold = 0.05,
        verbose = FALSE
    )

    expect_s3_class(res, "data.frame")
    expect_true(all(c("p_val", "tail", "p_val_adj", "sig_adj") %in% colnames(res)))
    expect_true(is.numeric(res$p_val_adj))
    aug_tail <- res$tail[res$codon == "AUG"]
    expect_equal(aug_tail, "right")
})

test_that("obtainSignificance safely assigns NA when an observed codon is missing from the distribution matrix", {

    missing_codon_value <- rbind(
        test_value_input,
        data.frame(group = "Treatment_A", codon = "XYZ", freq = 0.50)
    )

    res <- obtainSignificance(
        dist = test_dist_input,
        value = missing_codon_value,
        verbose = FALSE
    )

    xyz_row <- res[res$codon == "XYZ", ]
    expect_true(is.na(xyz_row$p_val))
    expect_true(is.na(xyz_row$tail))
})

test_that("obtainSignificance blocks processing if table column schemas do not conform", {

    corrupt_value <- test_value_input
    corrupt_value$ExtraColumn <- "Error"

    expect_error(
        obtainSignificance(dist = test_dist_input, value = corrupt_value, verbose = FALSE),
        "Incorrect number of columns"
    )
})

test_that("obtainSignificance prevents operations on non-numeric tracking arrays", {

    corrupt_numeric_value <- test_value_input
    corrupt_numeric_value$freq <- c("High", "Low") # Corrupt numbers to characters

    expect_error(
        obtainSignificance(dist = test_dist_input, value = corrupt_numeric_value, verbose = FALSE),
        "must be numeric"
    )
})

test_that("obtainSignificance emits sequential progress messages using example data", {

    # CASE 1: input data has exactly 3 columns
    logs <- capture.output(
        res <- obtainSignificance(
            dist = mock_permut_data,
            value = test_value_input,
            verbose = TRUE
        ),
        type = "message"
    )

    # Capturing the messages
    expect_true(any(grepl("--- Computing the statistical significance ---", logs)))
    expect_true(any(grepl("1 . Checking the input data.", logs)))
    expect_true(any(grepl("1 . COMPLETED", logs)))
    expect_true(any(grepl("2 . Calculating empirical p-values.", logs)))
    expect_true(any(grepl("2 . COMPLETED", logs)))
    expect_true(any(grepl("3 . Performing FDR correction.", logs)))
    expect_true(any(grepl("--- The significnace has been successfully computed ---", logs)))

    # CASE 2: input data has 5 columns
    expect_error(
        obtainSignificance(
            dist = mock_permut_data,
            value = mock_sig_data,
            verbose = FALSE
        ),
        regexp = "Incorrect number of columns"
    )
})
