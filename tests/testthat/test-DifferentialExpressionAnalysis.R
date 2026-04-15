data(default_tTEscanR_mRNA_data, default_tTEscanR_metadata)

test_that("Differential expression analysis executes correctly (multiple scenarios)", {

  # CASE 1: exploratory approach, controls
  expect_no_error(suppressWarnings(RunDEAnalysis(list_data = list(mRNA = mRNA_data_test), metadata = default_tTEscanR_metadata, batch = "tissue", color_factor = "tissue")))

  # CASE 2: missing required parameters
  expect_error(suppressWarnings(RunDEAnalysis(list_data = list(mRNA = mRNA_data_test), metadata = default_tTEscanR_metadata)), "specify the") # batch missing
  expect_error(suppressWarnings(RunDEAnalysis(list_data = mRNA_data_test, metadata = default_tTEscanR_metadata, batch = "tissue")), "No matching Sample IDs") # list_data not a named list
  expect_no_error(suppressWarnings(RunDEAnalysis(list_data = list(mRNA = mRNA_data_test), metadata = default_tTEscanR_metadata, batch = "tissue"))) # Generates only the heatmap
  expect_error(suppressWarnings(RunDEAnalysis(list_data = list(mRNA = mRNA_data_test), metadata = default_tTEscanR_metadata, batch = "tissue", dim_reduct = "PCA")), "specify a suitable 'color_factor'") # color_factor missing

  # CASE 3: targeted approach with two classes
  default_tTEscanR_metadata$targets <- "other"
  default_tTEscanR_metadata$targets[grep("endothelial", default_tTEscanR_metadata$cell.type)] <- "endothelial"
  expect_no_error(suppressWarnings(RunDEAnalysis(list_data = list(mRNA = default_tTEscanR_mRNA_data), metadata = default_tTEscanR_metadata, condition = "targets",
                                                 batch = "tissue", color_factor = "tissue", fc_threshold = 1, padj_threshold = 0.05)))

  # CASE 4: single class should fail
  default_tTEscanR_metadata$targets <- "other"
  expect_error(suppressWarnings(RunDEAnalysis(list_data = list(mRNA = mRNA_data_test), metadata = default_tTEscanR_metadata, condition = "targets", batch = "tissue")))

  # CASE 5: invalid data argument
  expect_error(RunDEAnalysis(list_data = mRNA_data_test, metadata = default_tTEscanR_metadata, batch = "tissue", dim_reduct = "PCA"))

  # CASE 6: metadata missing
  expect_error(RunDEAnalysis(list_data = list(mRNA = mRNA_data_test), batch = "tissue"))
})

DEA_results <- RunDEAnalysis(list_data = list(mRNA = default_tTEscanR_mRNA_data), metadata = default_tTEscanR_metadata, batch = "tissue", color_factor = "tissue")

test_that("Exploratory approach handles different dim_reduct methods", {

  # CASE 1: heatmap = TRUE + dim_reduct method
  expect_no_error(suppressWarnings(PlotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = TRUE, dim_reduct = "PCA", color_factor = "tissue")))
  expect_no_error(suppressWarnings(PlotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = TRUE, dim_reduct = "UMAP", color_factor = "tissue")))
  expect_no_error(suppressWarnings(PlotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = TRUE, dim_reduct = "tSNE", color_factor = "tissue")))

  # CASE 2: heatmap = FALSE + dim_reduct method
  expect_no_error(suppressWarnings(PlotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = FALSE, dim_reduct = "PCA", color_factor = "tissue")))
  expect_no_error(suppressWarnings(PlotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = FALSE, dim_reduct = "UMAP", color_factor = "tissue")))
  expect_no_error(suppressWarnings(PlotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = FALSE, dim_reduct = "tSNE", color_factor = "tissue")))

  # CASE 4: dim_reduct is NULL
  expect_no_error(suppressWarnings(PlotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = TRUE, dim_reduct = NULL, color_factor = "tissue")))
  expect_no_error(suppressWarnings(PlotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = FALSE, dim_reduct = NULL, color_factor = "tissue")))
})

test_that("Validating the number of principal components when dim_reduct = PCA", {

  # CASE 1: valid numPC
  expect_no_error(suppressWarnings(PlotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = TRUE, dim_reduct = "PCA", numPC = 3, color_factor = "tissue")))
  expect_no_error(suppressWarnings(PlotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = TRUE, dim_reduct = "PCA", color_factor = "tissue"))) # Default numPC = 2

  # CASE 2: invalid numPC, expect message, and use default numPC = 2
  expect_message(suppressWarnings(PlotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = TRUE, dim_reduct = "PCA", numPC = 0, color_factor = "tissue")), "Invalid 'numPC' given")
  expect_message(suppressWarnings(PlotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = TRUE, dim_reduct = "PCA", numPC = NULL, color_factor = "tissue")), "Invalid 'numPC' given")

})

test_that("Assessing the dim_reduct parameter validity", {

  # CASE 1: invalid dim_reduct should error
  expect_error(suppressWarnings(PlotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = FALSE, dim_reduct = PCA, color_factor = "tissue")))
  expect_error(suppressWarnings(PlotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = FALSE, dim_reduct = "pca", color_factor = "tissue")), "'arg' should be one of")
  expect_error(suppressWarnings(PlotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = FALSE, dim_reduct = "invalid", color_factor = "tissue")), "'arg' should be one of")
  expect_error(suppressWarnings(PlotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = FALSE, dim_reduct = c("PCA", "UMAP"), color_factor = "tissue")), "'arg' must be of length 1")

  # CASE 2: valid dim_reduct values
  for (method in c("PCA", "UMAP", "tSNE")) {
    expect_no_error(suppressWarnings(PlotDEResults(DE_results_list = DEA_results$results, dataset_name = "mRNA", heatmap = FALSE, dim_reduct = method, color_factor = "tissue")))
  }
})
