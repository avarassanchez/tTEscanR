data("default_tTEscanR_mRNA_data", package = "tTEscanR")
data("default_tTEscanR_metadata", package = "tTEscanR")

mRNA_data_test <- default_tTEscanR_mRNA_data[1:100, 1:100]

test_that("Differential expression analysis executes correctly (multiple scenarios)", {
    ## CASE 1: exploratory approach with defaults and custom dim_reduct
    expect_no_error(suppressWarnings(runDEAnalysis(list_data = list(mRNA = mRNA_data_test), metadata = default_tTEscanR_metadata, batch = "tissue", color_factor = "tissue", compute_pairwise = FALSE)))
    expect_no_error(suppressWarnings(runDEAnalysis(list_data = list(mRNA = mRNA_data_test), metadata = default_tTEscanR_metadata, batch = "tissue", dim_reduct = "PCA", compute_pairwise = FALSE)))

    ## CASE 2: Invalid input parameter guards
    expect_error(suppressWarnings(runDEAnalysis(list_data = list(mRNA = mRNA_data_test), metadata = default_tTEscanR_metadata, compute_pairwise = FALSE)), "specify the") # Missing batch
    expect_error(suppressWarnings(runDEAnalysis(list_data = mRNA_data_test, metadata = default_tTEscanR_metadata, batch = "tissue", compute_pairwise = FALSE)), "No matching Sample IDs") # Unnamed list / raw matrix input
    expect_error(runDEAnalysis(list_data = list(mRNA = mRNA_data_test), batch = "tissue", compute_pairwise = FALSE)) # Missing metadata

    ## CASE 3: Targeted approach with two classes vs single class failure
    targets_metadata <- default_tTEscanR_metadata
    targets_metadata$targets <- "other"
    targets_metadata$targets[grep("endothelial", targets_metadata$cell.type, fixed = TRUE)] <- "endothelial"

    expect_no_error(suppressWarnings(
        runDEAnalysis(
            list_data = list(mRNA = default_tTEscanR_mRNA_data), metadata = targets_metadata, target = "targets",
            batch = "tissue", color_factor = "tissue", fc_threshold = 1, padj_threshold = 0.05, compute_pairwise = FALSE
        )
    ))

    ## CASE 4: Single class target fails
    single_class_metadata <- default_tTEscanR_metadata
    single_class_metadata$targets <- "other"
    expect_error(suppressWarnings(runDEAnalysis(list_data = list(mRNA = mRNA_data_test), metadata = single_class_metadata, target = "targets", batch = "tissue", compute_pairwise = FALSE)))
})

DEA_results <- runDEAnalysis(list_data = list(mRNA = default_tTEscanR_mRNA_data), metadata = default_tTEscanR_metadata, batch = "tissue", color_factor = "tissue", compute_pairwise = FALSE)

test_that("Exploratory approach handles different dim_reduct methods", {
    ## CASE 1: Test valid reduction algorithms with heatmap = TRUE & FALSE
    for (method in c("PCA", "UMAP", "tSNE")) {
        expect_no_error(suppressWarnings(plotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = TRUE, dim_reduct = method, color_factor = "tissue")))
        expect_no_error(suppressWarnings(plotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = FALSE, dim_reduct = method, color_factor = "tissue")))
    }

    ## CASE 2: dim_reduct set to NULL
    expect_no_error(suppressWarnings(plotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = TRUE, dim_reduct = NULL, color_factor = "tissue")))
    expect_no_error(suppressWarnings(plotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = FALSE, dim_reduct = NULL, color_factor = "tissue")))
})

test_that("Validating dim_reduct arguments and numPC settings", {
    ## CASE 1: Valid numPC
    expect_no_error(suppressWarnings(plotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = TRUE, dim_reduct = "PCA", numPC = 3, color_factor = "tissue")))

    # CASE 2: Invalid numPC warning & fallback
    expect_warning(plotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = TRUE, dim_reduct = "PCA", numPC = 0, color_factor = "tissue"), "Invalid 'numPC' given")
    expect_warning(plotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = TRUE, dim_reduct = "PCA", numPC = NULL, color_factor = "tissue"), "Invalid 'numPC' given")

    # CASE 3: Invalid dim_reduct arguments (unquoted, invalid string, multiple values)
    expect_error(suppressWarnings(plotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = FALSE, dim_reduct = PCA, color_factor = "tissue")))
    expect_error(suppressWarnings(plotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = FALSE, dim_reduct = "invalid_method", color_factor = "tissue")), "'arg' should be one of")
    expect_error(suppressWarnings(plotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = FALSE, dim_reduct = c("PCA", "UMAP"), color_factor = "tissue")), "'arg' must be of length 1")})
