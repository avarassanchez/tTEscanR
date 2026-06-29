test_that("The default data is correctly loaded", {
    # CASE 1: no error - base case
    params <- list(
        list(s = "hg38", ref = codon_freq_table_canonical_hg38[rownames(codon_freq_table_canonical_hg38), ]),
        list(s = "mm39", ref = codon_freq_table_canonical_mm39[rownames(codon_freq_table_canonical_mm39), ])
    )

    for (p in params) {
        res <- selectDefaultData(p$s)
        expect_identical(res, p$ref)
    }

    expect_no_error(selectDefaultData(species = "mm39"))
    expect_no_error(selectDefaultData(species = "hg38"))

    # CASE 2: error - wrong or missing parameters
    expect_error(selectDefaultData(species = "mm10"))
    expect_error(selectDefaultData(species = "HG38"))
})

test_that("The data loaded has a correct format", {
    # CONTROLS
    expect_no_error(identifyInputFormat(data = mRNA_data_test, mode = "flexible"))
    expect_no_error(identifyInputFormat(data = mRNA_data_test, mode = "fix"))
    expect_no_error(identifyInputFormat(data = ENSG_gene_names_mRNA_data, mode = "flexible"))
    expect_no_error(identifyInputFormat(data = ENSG_gene_names_mRNA_data, mode = "fix"))

    # CASE 1: error - wrong input formats
    expect_no_error(identifyInputFormat(data = list(4, 5, 6, 7, 8), mode = "flexible"))
    expect_error(identifyInputFormat(data = list(4, 5, 6, 7, 8), mode = "fix")) # the fix parameters requires that the elements are matrices

    # CASE 2: error - empty input data
    expect_error(identifyInputFormat(data = c(), mode = "flexible"))
    expect_error(identifyInputFormat(data = c(), mode = "fix"))

    # CASE 3: missing mode parameter
    expect_no_error(identifyInputFormat(data = mRNA_data_test)) # mode fix is selected by default
    expect_no_error(identifyInputFormat(data = ENSG_gene_names_mRNA_data)) # mode fix is selected by default
    expect_error(identifyInputFormat(data = list(4, 5, 6, 7, 8)))

    # CASE 4: error - wrong mode parameter
    expect_error(identifyInputFormat(data = mRNA_data_test, mode = flexible)) # mode needs to be a string
    expect_error(identifyInputFormat(data = mRNA_data_test, mode = "strict")) # mode has to be flexible or fix
})

test_that("The individual columns are properly groupped", {
    # CASE 1 - Basic dense matrix
    data <- data.frame(sample_1 = c(1, 2), sample_2 = c(3, 4), sample_3 = c(5, 6))
    rownames(data) <- c("gene1", "gene2")
    groups <- c("A", "A", "B")

    result <- groupConditions(data, groups)

    expect_equal(colnames(result), c("A", "B"))
    expect_equal(result$A, c(4, 6)) # 1+3, 2+4
    expect_equal(result$B, c(5, 6))

    # CASE 2 - Sparse matrix
    data_sparse <- Matrix::Matrix(c(1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0), nrow = 2, ncol = 6, byrow = TRUE)
    rownames(data_sparse) <- c("gene1", "gene2")
    groups <- c("A", "A", "A", "B", "B", "B")

    result <- groupConditions(data_sparse, groups)

    expect_equal(colnames(result), c("A", "B"))
    expect_equal(result$A, c(3, 9)) # 1+0+2, 4+0+5
    expect_equal(result$B, c(3, 6)) # 0+3+0, 0+6+0

    # Create a sparse matrix (column-major by default)
    data_sparse <- Matrix::Matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 6)
    rownames(data_sparse) <- c("gene1", "gene2")
    groups <- c("A", "A", "A", "B", "B", "B")

    result <- groupConditions(data_sparse, groups)

    expected <- sapply(unique(groups), function(g) Matrix::rowSums(data_sparse[, groups == g, drop = FALSE]))
    expected <- as.data.frame(expected)
    colnames(expected) <- unique(groups)
    rownames(expected) <- rownames(data_sparse)
    expect_equal(result, expected)

    # CASE 3 - Single group
    data <- data.frame(sample_1 = c(1, 2), sample_2 = c(3, 4))
    rownames(data) <- c("gene1", "gene2")
    groups <- c("A", "A")

    result <- groupConditions(data, groups)

    expect_equal(ncol(result), 1)
    expect_equal(result$A, c(4, 6))

    # CASE 4 - Unique groups (no repeats)
    data <- data.frame(sample_1 = c(1, 2), sample_2 = c(3, 4))
    rownames(data) <- c("gene1", "gene2")
    groups <- c("A", "B")

    result <- groupConditions(data, groups)

    expect_equal(colnames(result), c("A", "B"))
    expect_equal(result$A, c(1, 2))
    expect_equal(result$B, c(3, 4))

    # CASE 5 - Non-numeric data fails gracefully
    data <- data.frame(sample_1 = c("a", "b"), sample_2 = c("c", "d"))
    rownames(data) <- c("gene1", "gene2")
    groups <- c("A", "B")

    result <- tryCatch(groupConditions(data, groups), error = function(e) e)

    expect_true(inherits(result, "error"))

    # CASE 6 - Dealing with missing values
    data <- data.frame(sample_1 = c(1, NA), sample_2 = c(3, 4))
    rownames(data) <- c("gene1", "gene2")
    groups <- c("A", "A")

    result <- groupConditions(data, groups)

    expect_equal(result$A, c(4, NA)) # NA propagates

    # CASE 6 - Empty matrix
    data <- data.frame()
    groups <- character(0)

    result <- tryCatch(groupConditions(data, groups), error = function(e) e)

    expect_true(inherits(result, "error"))

    # CASE 7 - "Large random matrix
    set.seed(42)
    data <- matrix(runif(1000 * 20), nrow = 1000, byrow = TRUE)
    groups <- sample(c("A", "B", "C", "D"), 20, replace = TRUE)

    result <- groupConditions(data, groups)

    expect_equal(nrow(result), 1000)
    expect_equal(ncol(result), length(unique(groups)))
})

mock_matrix <- matrix(
    c(10, 20, 30,  # Sample1
      40, 50, 60), # Sample2
    nrow = 3,
    ncol = 2,
    dimnames = list(c("GeneA", "GeneB", "GeneC"), c("Sample1", "Sample2"))
)

test_that("transformFormat converts a wide matrix into a long data frame cleanly", {
    # CASE 1
    res <- transformFormat(
        data = mock_matrix,
        normalize = FALSE,
        rownames_to_column = "FeatureID",
        names_to = "SampleID",
        values_to = "Counts"
    )
    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), 6) # 3 rows * 2 columns = 6 long records
    expect_equal(colnames(res), c("FeatureID", "SampleID", "Counts"))
    expect_equal(res$FeatureID, rep(c("GeneA", "GeneB", "GeneC"), 2))
    expect_equal(res$SampleID, c(rep("Sample1", 3), rep("Sample2", 3)))
    expect_equal(res$Counts, c(10, 20, 30, 40, 50, 60))

    # CASE 2
    res <- transformFormat(
        data = mock_matrix,
        normalize = TRUE,
        rownames_to_column = "FeatureID",
        names_to = "SampleID",
        values_to = "Proportion"
    )
    sample1_props <- res$Proportion[res$SampleID == "Sample1"]
    expect_equal(sum(sample1_props), 1.0)
    expect_equal(sample1_props[1], 10 / 60)

    sample2_props <- res$Proportion[res$SampleID == "Sample2"]
    expect_equal(sum(sample2_props), 1.0)
    expect_equal(sample2_props[1], 40 / 150)

    # CASE 3
    matrix_with_zero <- mock_matrix
    matrix_with_zero[, "Sample1"] <- 0 # Wipe out all counts in sample 1

    expect_silent({
        res <- transformFormat(
            data = matrix_with_zero,
            normalize = TRUE,
            rownames_to_column = "FeatureID",
            names_to = "SampleID",
            values_to = "Proportion"
        )
    })
    sample1_props <- res$Proportion[res$SampleID == "Sample1"]
    expect_false(any(is.nan(sample1_props)))
    expect_equal(sample1_props, c(0, 0, 0))

    # CASE 4
    df_input <- as.data.frame(mock_matrix)

    expect_silent({
        res <- transformFormat(
            data = df_input,
            normalize = FALSE,
            rownames_to_column = "FeatureID",
            names_to = "SampleID",
            values_to = "Counts"
        )
    })
    expect_equal(nrow(res), 6)
})

matrix_A <- matrix(
    c(5, 0, 0, 12),
    nrow = 2, ncol = 2,
    dimnames = list(c("Gene_A", "Gene_B"), c("Sample_1", "Sample_2"))
)

matrix_B <- matrix(
    c(7,  0, 0, 15),
    nrow = 2, ncol = 2,
    dimnames = list(c("Gene_B", "Gene_C"), c("Sample_2", "Sample_3"))
)

test_that("mergeMatrices cleanly joins non-overlapping matrix domains into a unified sparse matrix", {
    matrix_C <- matrix( # Does not share any row or columns
        c(99), nrow = 1, ncol = 1,
        dimnames = list("Gene_C", "Sample_3")
    )

    res <- mergeMatrices(matrix_A, matrix_C)

    expect_s4_class(res, "dgCMatrix") # Confirm structural sparse type
    expect_equal(nrow(res), 3) # 3 genes
    expect_equal(ncol(res), 3) # 3 samples
    expect_equal(as.numeric(res["Gene_A", "Sample_1"]), 5) # Proper coords
    expect_equal(as.numeric(res["Gene_C", "Sample_3"]), 99) # Proper coords
    expect_equal(as.numeric(res["Gene_A", "Sample_3"]), 0) # Fill with 0s
})

test_that("mergeMatrices successfully handles and forces alternative storage structures", {
    mixed_matrix <- matrix(
        c("10", "0", "0", "20"), nrow = 2,
        dimnames = list(c("Gene_A", "Gene_B"), c("Sample_1", "Sample_2"))
    )

    expect_silent({
        res <- mergeMatrices(mixed_matrix)
    })

    expect_equal(as.numeric(res["Gene_A", "Sample_1"]), 10)
    expect_equal(as.numeric(res["Gene_B", "Sample_2"]), 20)
})

test_that("mergeMatrices handles all-zero matrices without crashing coordinate counters", {
    zero_matrix <- matrix(
        c(0, 0), nrow = 2, ncol = 1,
        dimnames = list(c("Gene_A", "Gene_B"), c("Sample_Null"))
    )

    expect_silent({
        res <- mergeMatrices(matrix_A, zero_matrix)
    })

    expect_true("Sample_Null" %in% colnames(res))
    expect_equal(as.numeric(sum(res[, "Sample_Null"])), 0)
})

test_that("mergeMatrices aggregates intersecting coordinate cells by summation", {
    res <- mergeMatrices(matrix_A, matrix_B)

    expect_equal(as.numeric(res["Gene_B", "Sample_2"]), 19)
    expect_equal(as.numeric(res["Gene_A", "Sample_1"]), 5)
    expect_equal(as.numeric(res["Gene_C", "Sample_3"]), 15)
})

tRNA_test_matrix <- matrix(
    c(100, 200, 4000, 6000),
    nrow = 2, ncol = 2,
    dimnames = list(c("tRNA-Ala", "tRNA-Gly"), c("Sample1", "Sample2"))
)

tRNA_test_matrix_2 <- matrix(
    c(100, 200, 50, 150),
    nrow = 2, ncol = 2,
    dimnames = list(c("tRNA-Ala", "tRNA-Gly"), c("Sample1", "Sample2"))
)

single_keep_matrix <- matrix(
    c(6000, 4000),
    nrow = 2, ncol = 1,
    dimnames = list(c("tRNA-Ala", "tRNA-Gly"), c("Sample_High"))
)

test_that("tRNAFilterCuts filters out accurately", {
    # CASE 1
    expect_message(
        filtered_data <- tRNAFilterCuts(tRNA_test_matrix, cutoff = 5000, verbose = TRUE),
        "Initial number of samples"
    )
    expect_equal(ncol(filtered_data), 1)
    expect_equal(colnames(filtered_data), "Sample2")
    expect_equal(as.numeric(filtered_data["tRNA-Ala", "Sample2"]), 4000)

    # CASE 2
    expect_error(
        tRNAFilterCuts(tRNA_test_matrix_2, cutoff = 5000, verbose = FALSE),
        "No samples pass the cutoff"
    )

    # CASE 3
    res <- tRNAFilterCuts(single_keep_matrix, cutoff = 5000, verbose = FALSE)

    expect_true(is.matrix(res))
    expect_equal(dim(res), c(2, 1))
})

test_that("tRNASetGenes handles customizable separators and boundary clamping", {

    # Create mock BED file
    mock_bed_contents <- c("chr1\t500\t600\ttRNA-Ala")
    tmp_bed <- withr::local_tempfile()
    writeLines(mock_bed_contents, tmp_bed)

    # CASE 1: Using the default c("-", "-") format -> "chr1-401-700"
    mock_matrix_dash <- matrix(1, nrow = 1, ncol = 1)
    rownames(mock_matrix_dash) <- "chr1-401-700"

    res_dash <- tRNASetGenes(mock_matrix_dash, tmp_bed, flanking_region = 100)
    expect_equal(rownames(res_dash), "tRNA-Ala")

    # CASE 2: Custom mismatched separators c(":", "-") -> "chr1:401-700"
    mock_matrix_custom <- matrix(1, nrow = 1, ncol = 1)
    rownames(mock_matrix_custom) <- "chr1:401-700"

    res_custom <- tRNASetGenes(
        data = mock_matrix_custom,
        tRNA_bed = tmp_bed,
        flanking_region = 100,
        name_sep = c(":", "-")
    )
    expect_equal(rownames(res_custom), "tRNA-Ala")
})

mock_ref_output <- list(
    anticodon = c("Ala-AGC" = 120, "Gly-GCC" = 95, "Val-CAC" = 40),
    supply = c("Ala" = 120, "Gly" = 95, "Val" = 40)
)

mock_iterations_output <- list( # Example of runIterations( ) output
    data.frame(
        optimal_cutoff = c(100, 100),
        type = c("anticodon", "supply"),
        stringsAsFactors = FALSE
    )
)

mock_cor_output <- data.frame( # Example of computeCorrelations() output
    cutoff = c(50, 51, 52),
    anticodon_spearman = c(0.96, 0.97, 0.98),
    supply_spearman = c(0.95, 0.96, 0.96),
    total_anticodon = c(255, 255, 255),
    total_supply = c(255, 255, 255),
    stringsAsFactors = FALSE
)

test_that("tRNASetCutoff executes full pipeline, prints verbose logs, and attaches plots", {

    mockery::stub(tRNASetCutoff, "referenceObject", mock_ref_output)
    mockery::stub(tRNASetCutoff, "transformCounts", data.frame(features = "dummy"))
    mockery::stub(tRNASetCutoff, "runIterations", mock_iterations_output)
    mockery::stub(tRNASetCutoff, "computeCorrelations", mock_cor_output)
    mockery::stub(tRNASetCutoff, "correlationCutoffPlot", "mocked_correlation_ggplot_object")
    mockery::stub(tRNASetCutoff, "selectionCutoffPlot", "mocked_histogram_ggplot_object")

    mock_counts <- matrix(
        c(100, 150, 200, 50, 80, 110),
        nrow = 3,
        ncol = 2,
        dimnames = list(c("tRNA-Ala-AGC", "tRNA-Gly-GCC", "tRNA-Val-CAC"), c("Sample_A", "Sample_B"))
    )

    expect_message({
        res <- tRNASetCutoff(
            data = mock_counts,
            num_iter = 2,
            cutoffs_limits = c(100, 105),
            generate_plot = TRUE,
            compute_aa = TRUE,
            verbose = TRUE
        )
    }, "1 . Computing reference tTEscanR object.")

    expect_type(res, "list")
    expect_named(res, c("optimal_cutoff", "correlation_plot", "histogram_plot"))
    expect_equal(res$correlation_plot, "mocked_correlation_ggplot_object")
    expect_equal(res$histogram_plot, "mocked_histogram_ggplot_object")
    expect_s3_class(res$optimal_cutoff, "data.frame")
    expect_true("iteration" %in% colnames(res$optimal_cutoff))
})

test_that("tRNASetCutoff drops plotting branches smoothly when generate_plot = FALSE", {

    mockery::stub(tRNASetCutoff, "referenceObject", mock_ref_output)
    mockery::stub(tRNASetCutoff, "transformCounts", data.frame(features = "dummy"))
    mockery::stub(tRNASetCutoff, "runIterations", mock_iterations_output)

    mock_counts <- matrix(c(100, 200), nrow = 2, ncol = 1)

    res <- tRNASetCutoff(
        data = mock_counts,
        num_iter = 1,
        cutoffs_limits = c(50, 52),
        generate_plot = FALSE,
        compute_aa = FALSE,
        verbose = FALSE
    )

    expect_type(res, "list")
    expect_named(res, "optimal_cutoff")
    expect_false("correlation_plot" %in% names(res))
    expect_false("histogram_plot" %in% names(res))
})

test_that("transformCounts converts a wide count matrix into an uncounted long data frame", {

    # CASE 1 - minimal count matrix example
    mock_matrix <- matrix(
        c(2, 0, 0, 1),
        nrow = 2,
        dimnames = list(c("tRNA-1", "tRNA-2"), c("Sample_A", "Sample_B"))
    )
    res <- transformCounts(mock_matrix)
    expect_equal(nrow(res), 3)
    expect_true(all(c("features", "conditions") %in% colnames(res)))
    expect_false("counts" %in% colnames(res))
    expect_equal(sum(res$features == "tRNA-1" & res$conditions == "Sample_A"), 2)
    expect_equal(sum(res$features == "tRNA-2" & res$conditions == "Sample_B"), 1)

    # CASE 2
    res <- transformCounts(tRNA_test_matrix)

    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), 10300)
    expect_true(all(c("features", "conditions") %in% colnames(res)))
    expect_false("counts" %in% colnames(res))
    expect_equal(sum(res$features == "tRNA-Ala" & res$conditions == "Sample1"), 100)
    expect_equal(sum(res$features == "tRNA-Gly" & res$conditions == "Sample2"), 6000)
    expect_equal(sum(res$features == "tRNA-1" & res$conditions == "Sample_A"), 0)
    expect_equal(sum(res$features == "tRNA-1" & res$conditions == "Sample1"), 0)
})


test_that("computeCorrelations accurately calculates spearman profiles against references", {
    ref_anticodon <- c("Ala-AGC" = 10, "Gly-GCC" = 20, "Val-CAC" = 30)
    ref_supply    <- c("Ala" = 10, "Gly" = 20, "Val" = 30)

    mock_filtered_obj <- list(
        anticodon = c("Ala-AGC" = 10, "Gly-GCC" = 20, "Val-CAC" = 30),
        supply    = c("Ala" = 10, "Gly" = 20, "Val" = 30)
    )
    mockery::stub(computeCorrelations, "filteringCutoffs", mock_filtered_obj)

    res <- computeCorrelations(
        data = "dummy_data",
        ref_anticodon = ref_anticodon,
        ref_supply = ref_supply,
        cutoffs = c(100, 200)
    )

    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), 2)
    expect_named(res, c("cutoff", "anticodon_spearman", "supply_spearman", "total_anticodon", "total_supply"))
    expect_equal(res$anticodon_spearman[1], 1.0)
    expect_equal(res$supply_spearman[1], 1.0)
    expect_equal(res$total_anticodon[1], 60)
    expect_equal(res$total_supply[1], 60)
})

mock_codon_table <- matrix(
    rpois(40, lambda = 50), nrow = 4, ncol = 10,
    dimnames = list(c("AUU", "AUC", "AUA", "AUG"), paste0("Gene", 1:10))
)

test_that("referenceObject builds structured list frameworks based on compute_aa", {
    mock_anticodon_vec <- c("AUU" = 50, "AUC" = 60)
    mock_supply_vec    <- c("Ile" = 110)

    mockery::stub(referenceObject, "CreateObject", "dummy_s4_obj")
    mockery::stub(referenceObject, "computeAnticodonUsage", "dummy_s4_obj")
    mockery::stub(referenceObject, "getAssay", "dummy_assay_obj")
    mockery::stub(referenceObject, "getReference", mock_anticodon_vec)

    # CASE 1 - compute_aa set to FALSE
    res_no_aa <- referenceObject(data = mock_codon_table, compute_aa = FALSE)
    expect_type(res_no_aa, "list")
    expect_named(res_no_aa, c("anticodon", "supply"))
    expect_equal(res_no_aa$anticodon, mock_anticodon_vec)
    expect_null(res_no_aa$supply)

    # CASE 2 - compute_aa set to TRUE
    mock_get_ref <- mockery::mock(mock_anticodon_vec, mock_supply_vec)
    mockery::stub(referenceObject, "getReference", mock_get_ref)
    mockery::stub(referenceObject, "computeAAUsage", "dummy_s4_obj")
    res_with_aa <- referenceObject(data = mock_codon_table, compute_aa = TRUE)
    expect_equal(res_with_aa$supply, mock_supply_vec)
})

test_that("runIterations loops properly and updates the progress bar", {
    mock_iter_output <- data.frame( # Single iteration output of iterateCutofftRNA
        optimal_cutoff = 100,
        type = "anticodon",
        stringsAsFactors = FALSE
    )

    mockery::stub(runIterations, "iterateCutofftRNA", mock_iter_output)
    mockery::stub(runIterations, "utils::txtProgressBar", "mock_pb")
    mockery::stub(runIterations, "utils::setTxtProgressBar", NULL)
    mockery::stub(runIterations, "close", NULL)

    res_iterations <- runIterations(
        num_iter = 3,
        data = "dummy_transformed_data",
        anticodon = c("AUU" = 50),
        supply = NULL,
        cuts = c(50, 100),
        slope = 0.001,
        rho = 0.95
    )

    expect_type(res_iterations, "list")
    expect_length(res_iterations, 3) # Length must match num_iter exactly
    expect_equal(res_iterations[[1]], mock_iter_output)
    expect_equal(res_iterations[[2]], mock_iter_output)
    expect_equal(res_iterations[[3]], mock_iter_output)
})

mock_cor_data <- data.frame(
    cutoff = c(50, 100),
    anticodon_spearman = c(0.95, 0.98),
    supply_spearman = c(0.91, 0.94),
    total_anticodon = c(500, 500),
    total_supply = c(500, 500),
    stringsAsFactors = FALSE
)

mock_boot_data <- data.frame(
    optimal_cutoff = c(100, 100, 150),
    type = c("anticodon", "supply", "anticodon"),
    stringsAsFactors = FALSE
)

test_that("correlationCutoffPlot builds a ggplot with accurate internal mappings", {

    mockery::stub(correlationCutoffPlot, "savePlot", NULL)

    p <- correlationCutoffPlot(data = mock_cor_data, add_titles = TRUE)
    expect_s3_class(p, "ggplot")
    layer_geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    expect_true("GeomPoint" %in% layer_geoms)

    expect_equal(rlang::as_label(p$mapping$x), "cutoff")
    expect_equal(rlang::as_label(p$mapping$y), "spearman_corr")
    expect_equal(rlang::as_label(p$mapping$colour), "type")
    expect_equal(p$labels$title, "Correlation vs Reference Across Cutoff Thresholds")
})

test_that("correlationCutoffPlot drops labels and executes export routes optionally", {
    mock_cor_data <- data.frame(cutoff = 50, anticodon_spearman = 0.9, supply_spearman = 0.8)

    mock_saver <- mockery::mock()
    mockery::stub(correlationCutoffPlot, "savePlot", mock_saver)

    p <- correlationCutoffPlot(
        data = mock_cor_data,
        add_titles = FALSE,
        save_format = "pdf",
        out_name = "test_plot"
    )

    expect_null(p$labels$title)
    mockery::expect_called(mock_saver, 1)
})

test_that("selectionCutoffPlot creates a faceted bar layout structured for distributions", {
    mockery::stub(selectionCutoffPlot, "savePlot", NULL)

    p <- selectionCutoffPlot(data = mock_boot_data, add_titles = TRUE)

    expect_s3_class(p, "ggplot")
    layer_geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    expect_true("GeomBar" %in% layer_geoms)
    expect_equal(rlang::as_label(p$mapping$x), "optimal_cutoff")
    expect_equal(rlang::as_label(p$mapping$fill), "type")
    expect_s3_class(p$facet, "FacetWrap")
    expect_equal(names(p$facet$params$facets), "type")

    expect_equal(p$labels$title, "Distribution of Cutoff Values")
})
