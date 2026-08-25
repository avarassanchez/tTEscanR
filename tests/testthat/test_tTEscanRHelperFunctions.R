data("default_tTEscanR_mRNA_data", package = "tTEscanR")
ENSG_gene_names <- c("ENSG00000162104", "ENSG00000058799", "ENSG00000184988", "ENSG00000131459")

test_that("selectDefaultData correctly retrieves canonical species tables", {
    ## CASE 1: Valid species checks
    expect_identical(selectDefaultData("hg38"), codon_freq_table_canonical_hg38)
    expect_identical(selectDefaultData("mm39"), codon_freq_table_canonical_mm39)

    ## CASE 2: Unsupported or misformatted species identifiers
    expect_error(selectDefaultData(species = "mm10"))
    expect_error(selectDefaultData(species = "HG38"))
})

test_that("identifyInputFormat validates input structures and mode arguments", {
    ## CASE 1: Valid inputs under flexible and default (fix) modes
    expect_no_error(identifyInputFormat(data = default_tTEscanR_mRNA_data, mode = "flexible"))
    expect_no_error(identifyInputFormat(data = default_tTEscanR_mRNA_data, mode = "fix"))
    expect_no_error(identifyInputFormat(data = ENSG_gene_names, mode = "flexible"))
    expect_no_error(identifyInputFormat(data = default_tTEscanR_mRNA_data))

    ## CASE 2: Invalid input formats
    expect_no_error(identifyInputFormat(data = list(4, 5, 6, 7, 8), mode = "flexible"))
    expect_error(identifyInputFormat(data = list(4, 5, 6, 7, 8), mode = "fix"))
    expect_error(identifyInputFormat(data = c(), mode = "flexible"))
    expect_error(identifyInputFormat(data = c(), mode = "fix"))

    ## CASE 3: Invalid mode parameters
    expect_error(identifyInputFormat(data = default_tTEscanR_mRNA_data, mode = flexible))
    expect_error(identifyInputFormat(data = default_tTEscanR_mRNA_data, mode = "strict"))
})

test_that("groupConditions properly aggregates dense, sparse, and edge-case matrices", {
    ## CASE 1: Dense matrix aggregation
    dense_df <- data.frame(sample_1 = c(1, 2), sample_2 = c(3, 4), sample_3 = c(5, 6), row.names = c("gene1", "gene2"))
    res_dense <- groupConditions(dense_df, c("A", "A", "B"))
    expect_equal(colnames(res_dense), c("A", "B"))
    expect_equal(res_dense$A, c(4, 6))
    expect_equal(res_dense$B, c(5, 6))

    ## CASE 2: Sparse matrix aggregation
    sparse_mat <- Matrix::Matrix(c(1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0), nrow = 2, ncol = 6, byrow = TRUE, dimnames = list(c("gene1", "gene2"), NULL))
    res_sparse <- groupConditions(sparse_mat, c("A", "A", "A", "B", "B", "B"))
    expect_equal(res_sparse$A, c(3, 9))
    expect_equal(res_sparse$B, c(3, 6))

    ## CASE 3: Column-major sparse matrix rowSums match
    sparse_cm <- Matrix::Matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 6, dimnames = list(c("gene1", "gene2"), NULL))
    res_cm <- groupConditions(sparse_cm, c("A", "A", "A", "B", "B", "B"))
    expected_cm <- as.data.frame(sapply(c("A", "B"), function(g) Matrix::rowSums(sparse_cm[, c("A", "A", "A", "B", "B", "B") == g, drop = FALSE])))
    expect_equal(res_cm, expected_cm)

    ## CASE 4: Unique & Single Group cases
    res_single <- groupConditions(dense_df[, 1:2], c("A", "A"))
    expect_equal(res_single$A, c(4, 6))

    res_unique <- groupConditions(dense_df[, 1:2], c("A", "B"))
    expect_equal(colnames(res_unique), c("A", "B"))

    ## CASE 5: NA handling and invalid data errors
    df_na <- data.frame(sample_1 = c(1, NA), sample_2 = c(3, 4), row.names = c("gene1", "gene2"))
    expect_equal(groupConditions(df_na, c("A", "A"))$A, c(4, NA))

    expect_error(groupConditions(data.frame(sample_1 = c("a", "b"), sample_2 = c("c", "d")), c("A", "B")))
    expect_error(groupConditions(data.frame(), character(0)))

    ## CASE 6: Large random matrix test
    set.seed(42)
    large_mat <- matrix(runif(1000 * 20), nrow = 1000, byrow = TRUE)
    groups <- sample(c("A", "B", "C", "D"), 20, replace = TRUE)
    res_large <- groupConditions(large_mat, groups)
    expect_equal(dim(res_large), c(1000, length(unique(groups))))
})

test_that("transformFormat converts wide matrices to long data frames accurately", {
    mock_matrix <- matrix(
        c(10, 20, 30, 40, 50, 60),
        nrow = 3, ncol = 2,
        dimnames = list(c("GeneA", "GeneB", "GeneC"), c("Sample1", "Sample2"))
    )

    ## CASE 1: Standard unnormalized conversion
    res_unnorm <- transformFormat(mock_matrix, normalize = FALSE, rownames_to_column = "FeatureID", names_to = "SampleID", values_to = "Counts")
    expect_s3_class(res_unnorm, "data.frame")
    expect_equal(nrow(res_unnorm), 6)
    expect_equal(colnames(res_unnorm), c("FeatureID", "SampleID", "Counts"))
    expect_equal(res_unnorm$Counts, c(10, 20, 30, 40, 50, 60))

    ## CASE 2: Normalized proportions
    res_norm <- transformFormat(mock_matrix, normalize = TRUE, rownames_to_column = "FeatureID", names_to = "SampleID", values_to = "Proportion")
    expect_equal(sum(res_norm$Proportion[res_norm$SampleID == "Sample1"]), 1.0)
    expect_equal(res_norm$Proportion[res_norm$SampleID == "Sample1"][1], 10 / 60)

    ## CASE 3: Zero counts normalization handling
    zero_mat <- mock_matrix
    zero_mat[, "Sample1"] <- 0
    expect_silent(res_zero <- transformFormat(zero_mat, normalize = TRUE, rownames_to_column = "FeatureID", names_to = "SampleID", values_to = "Proportion"))
    expect_equal(res_zero$Proportion[res_zero$SampleID == "Sample1"], c(0, 0, 0))

    ## CASE 4: Data frame input compatibility
    expect_silent(res_df <- transformFormat(as.data.frame(mock_matrix), normalize = FALSE, rownames_to_column = "FeatureID", names_to = "SampleID", values_to = "Counts"))
    expect_equal(nrow(res_df), 6)
})

test_that("mergeMatrices handles disjoint domains, zero counts, and cell aggregation", {
    matrix_A <- matrix(c(5, 0, 0, 12), nrow = 2, ncol = 2, dimnames = list(c("Gene_A", "Gene_B"), c("Sample_1", "Sample_2")))
    matrix_B <- matrix(c(7, 0, 0, 15), nrow = 2, ncol = 2, dimnames = list(c("Gene_B", "Gene_C"), c("Sample_2", "Sample_3")))
    matrix_C <- matrix(c(99), nrow = 1, ncol = 1, dimnames = list("Gene_C", "Sample_3"))

    ## CASE 1: Join non-overlapping domains into sparse matrix
    res_disjoint <- mergeMatrices(matrix_A, matrix_C)
    expect_s4_class(res_disjoint, "dgCMatrix")
    expect_equal(dim(res_disjoint), c(3, 3))
    expect_equal(as.numeric(res_disjoint["Gene_A", "Sample_1"]), 5)
    expect_equal(as.numeric(res_disjoint["Gene_C", "Sample_3"]), 99)
    expect_equal(as.numeric(res_disjoint["Gene_A", "Sample_3"]), 0)

    ## CASE 2: Convert character matrices safely
    mixed_matrix <- matrix(c("10", "0", "0", "20"), nrow = 2, dimnames = list(c("Gene_A", "Gene_B"), c("Sample_1", "Sample_2")))
    expect_silent(res_mixed <- mergeMatrices(mixed_matrix))
    expect_equal(as.numeric(res_mixed["Gene_A", "Sample_1"]), 10)

    ## CASE 3: All-zero matrix merging
    zero_matrix <- matrix(c(0, 0), nrow = 2, ncol = 1, dimnames = list(c("Gene_A", "Gene_B"), c("Sample_Null")))
    expect_silent(res_zero <- mergeMatrices(matrix_A, zero_matrix))
    expect_true("Sample_Null" %in% colnames(res_zero))
    expect_equal(as.numeric(sum(res_zero[, "Sample_Null"])), 0)

    ## CASE 4: Sum intersecting coordinates
    res_agg <- mergeMatrices(matrix_A, matrix_B)
    expect_equal(as.numeric(res_agg["Gene_B", "Sample_2"]), 19)
    expect_equal(as.numeric(res_agg["Gene_A", "Sample_1"]), 5)
    expect_equal(as.numeric(res_agg["Gene_C", "Sample_3"]), 15)
})

test_that("tRNAFilterCuts filters sample count cutoffs correctly", {
    tRNA_test_matrix <- matrix(c(100, 200, 4000, 6000), nrow = 2, ncol = 2, dimnames = list(c("tRNA-Ala", "tRNA-Gly"), c("Sample1", "Sample2")))
    tRNA_test_matrix_2 <- matrix(c(100, 200, 50, 150), nrow = 2, ncol = 2, dimnames = list(c("tRNA-Ala", "tRNA-Gly"), c("Sample1", "Sample2")))

    ## CASE 1: Successful sample retention
    expect_message(filtered_data <- tRNAFilterCuts(tRNA_test_matrix, cutoff = 5000, verbose = TRUE), "Initial number of samples")
    expect_equal(colnames(filtered_data), "Sample2")

    ## CASE 2: Exception when zero samples pass cutoff
    expect_error(tRNAFilterCuts(tRNA_test_matrix_2, cutoff = 5000, verbose = FALSE), "No samples pass the cutoff")

    ## CASE 3: Single column matrix retention
    single_keep_matrix <- matrix(c(6000, 4000), nrow = 2, ncol = 1, dimnames = list(c("tRNA-Ala", "tRNA-Gly"), c("Sample_High")))
    res_single <- tRNAFilterCuts(single_keep_matrix, cutoff = 5000, verbose = FALSE)
    expect_true(is.matrix(res_single))
    expect_equal(dim(res_single), c(2, 1))
})

test_that("tRNASetGenes handles customizable separators and boundary clamping", {
    tmp_bed <- withr::local_tempfile()
    writeLines("chr1\t500\t600\ttRNA-Ala", tmp_bed)

    ## CASE 1: Standard dash separator
    mock_dash <- matrix(1, nrow = 1, ncol = 1, dimnames = list("chr1-401-700", "S1"))
    expect_equal(rownames(tRNASetGenes(mock_dash, tmp_bed, flanking_region = 100)), "tRNA-Ala")

    ## CASE 2: Custom separator pair
    mock_custom <- matrix(1, nrow = 1, ncol = 1, dimnames = list("chr1:401-700", "S1"))
    expect_equal(rownames(tRNASetGenes(mock_custom, tmp_bed, flanking_region = 100, name_sep = c(":", "-"))), "tRNA-Ala")
})

test_that("tRNASetCutoff executes pipeline and formats output lists", {
    ref_out <- list(anticodon = c("Ala-AGC" = 120, "Gly-GCC" = 95), supply = c("Ala" = 120, "Gly" = 95))
    iter_out <- list(data.frame(optimal_cutoff = c(100, 100), type = c("anticodon", "supply"), stringsAsFactors = FALSE))
    cor_out <- data.frame(cutoff = c(50, 51), anticodon_spearman = c(0.96, 0.97), supply_spearman = c(0.95, 0.96), total_anticodon = c(255, 255), total_supply = c(255, 255))

    mockery::stub(tRNASetCutoff, "referenceObject", ref_out)
    mockery::stub(tRNASetCutoff, "transformCounts", data.frame(features = "dummy"))
    mockery::stub(tRNASetCutoff, "runIterations", iter_out)
    mockery::stub(tRNASetCutoff, "computeCorrelations", cor_out)
    mockery::stub(tRNASetCutoff, "correlationCutoffPlot", "mocked_correlation_plot")
    mockery::stub(tRNASetCutoff, "selectionCutoffPlot", "mocked_histogram_plot")

    mock_counts <- matrix(c(100, 150, 200, 50), nrow = 2, dimnames = list(c("tRNA-Ala-AGC", "tRNA-Gly-GCC"), c("S1", "S2")))

    ## CASE 1: Verbose execution with plots
    expect_message(res <- tRNASetCutoff(data = mock_counts, num_iter = 2, cutoffs_limits = c(100, 105), generate_plot = TRUE, compute_aa = TRUE, verbose = TRUE), "1 . Computing reference tTEscanR object.")
    expect_named(res, c("optimal_cutoff", "correlation_plot", "histogram_plot"))

    ## CASE 2: Execution without plots
    res_no_plots <- tRNASetCutoff(data = mock_counts, num_iter = 1, cutoffs_limits = c(50, 52), generate_plot = FALSE, compute_aa = FALSE, verbose = FALSE)
    expect_named(res_no_plots, "optimal_cutoff")
})

test_that("transformCounts converts wide counts into uncounted long formats", {
    mock_matrix <- matrix(c(2, 0, 0, 1), nrow = 2, dimnames = list(c("tRNA-1", "tRNA-2"), c("Sample_A", "Sample_B")))
    res <- transformCounts(mock_matrix)

    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), 3)
    expect_true(all(c("features", "conditions") %in% colnames(res)))
    expect_false("counts" %in% colnames(res))
    expect_equal(sum(res$features == "tRNA-1" & res$conditions == "Sample_A"), 2)
})

test_that("computeCorrelations calculates spearman profiles against references", {
    ref_anticodon <- c("Ala-AGC" = 10, "Gly-GCC" = 20, "Val-CAC" = 30)
    ref_supply    <- c("Ala" = 10, "Gly" = 20, "Val" = 30)

    mockery::stub(computeCorrelations, "filteringCutoffs", list(anticodon = ref_anticodon, supply = ref_supply))

    res <- computeCorrelations(
        data = "dummy",
        ref_anticodon = ref_anticodon,
        ref_supply = ref_supply,
        cutoffs = c(100, 200)
    )

    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), 2)
    expect_named(res, c("cutoff", "anticodon_spearman", "supply_spearman", "total_anticodon", "total_supply"))
    expect_equal(res$anticodon_spearman[1], 1.0)
    expect_equal(res$supply_spearman[1], 1.0)
})

test_that("referenceObject structures reference lists conditionally", {
    mock_codon_table <- matrix(
        1:40, nrow = 4, ncol = 10,
        dimnames = list(c("AUU", "AUC", "AUA", "AUG"), paste0("Gene", 1:10))
    )
    mock_ac <- c("AUU" = 50, "AUC" = 60)
    mock_sp <- c("Ile" = 110)

    # CASE 1: Without AA computation
    mockery::stub(referenceObject, "CreateObject", "dummy")
    mockery::stub(referenceObject, "computeAnticodonUsage", "dummy")
    mockery::stub(referenceObject, "SummarizedExperiment::assay", "dummy")
    mockery::stub(referenceObject, "getReference", mock_ac)

    res_no_aa <- referenceObject(data = mock_codon_table, compute_aa = FALSE)
    expect_named(res_no_aa, c("anticodon", "supply"))
    expect_null(res_no_aa$supply)

    # CASE 2: With AA computation
    mockery::stub(referenceObject, "CreateObject", "dummy")
    mockery::stub(referenceObject, "computeAnticodonUsage", "dummy")
    mockery::stub(referenceObject, "computeAAUsage", "dummy")
    mockery::stub(referenceObject, "SummarizedExperiment::assay", "dummy")
    mockery::stub(referenceObject, "getReference", mockery::mock(mock_ac, mock_sp))

    res_aa <- referenceObject(data = mock_codon_table, compute_aa = TRUE)
    expect_equal(res_aa$supply, mock_sp)
})
test_that("runIterations iterates across cutoffs and manages progress bar", {
    mock_iter <- data.frame(optimal_cutoff = 100, type = "anticodon", stringsAsFactors = FALSE)

    mockery::stub(runIterations, "iterateCutofftRNA", mock_iter)
    mockery::stub(runIterations, "utils::txtProgressBar", "mock_pb")
    mockery::stub(runIterations, "utils::setTxtProgressBar", NULL)
    mockery::stub(runIterations, "close", NULL)

    res <- runIterations(num_iter = 3, data = "dummy", anticodon = c("AUU" = 50), supply = NULL, cuts = c(50, 100), slope = 0.001, rho = 0.95)
    expect_type(res, "list")
    expect_length(res, 3)
    expect_equal(res[[1]], mock_iter)
})

test_that("Cutoff visualization functions construct expected ggplot objects", {
    cor_data <- data.frame(cutoff = c(50, 100), anticodon_spearman = c(0.95, 0.98), supply_spearman = c(0.91, 0.94), total_anticodon = c(500, 500), total_supply = c(500, 500))
    boot_data <- data.frame(optimal_cutoff = c(100, 100, 150), type = c("anticodon", "supply", "anticodon"))

    ## CASE 1: Correlation cutoff plot
    p_cor <- correlationCutoffPlot(data = cor_data, add_titles = TRUE)
    expect_s3_class(p_cor, "ggplot")
    expect_equal(p_cor$labels$title, "Correlation vs Reference Across Cutoff Thresholds")

    ## CASE 2: Correlation plot without titles
    p_cor_no_title <- correlationCutoffPlot(data = cor_data, add_titles = FALSE)
    expect_null(p_cor_no_title$labels$title)

    ## CASE 3: Selection cutoff plot
    p_sel <- selectionCutoffPlot(data = boot_data, add_titles = TRUE)
    expect_s3_class(p_sel, "ggplot")
    expect_s3_class(p_sel$facet, "FacetWrap")
    expect_equal(p_sel$labels$title, "Distribution of Cutoff Values")
})
