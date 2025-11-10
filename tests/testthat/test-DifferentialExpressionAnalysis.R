test_that("Differential expression analysis executes correctly (multiple scenarios)", {

  data(subset_mRNA_data, metadata)

  # CASE 1: exploratory approach, controls
  expect_no_error(suppressWarnings(RunDEAnalysis(list_data = list(mRNA = subset_mRNA_data), metadata = metadata, corr_factor = "tissue")))
  expect_no_error(suppressWarnings(RunDEAnalysis(list_data = list(mRNA = subset_mRNA_data), metadata = metadata, corr_factor = "tissue", color_factor = "tissue")))

  # CASE 2: missing required parameters
  expect_error(suppressWarnings(RunDEAnalysis(list_data = list(mRNA = subset_mRNA_data), metadata = metadata))) # corr_factor missing
  expect_error(suppressWarnings(RunDEAnalysis(list_data = subset_mRNA_data, metadata = metadata, corr_factor = "tissue"))) # list_data not a named list

  # CASE 3: targeted approach with two classes
  metadata$targets <- "other"
  metadata$targets[grep("endothelial", metadata$cell.type)] <- "endothelial"
  expect_no_error(suppressWarnings(RunDEAnalysis(list_data = list(mRNA = subset_mRNA_data), metadata = metadata,
                                                 target_factor = "targets", corr_factor = "tissue", fc_threshold = 1, pval_threshold = 0.05)))

  # CASE 4: single class should fail
  metadata$targets <- "other"
  expect_error(suppressWarnings(RunDEAnalysis(list_data = list(mRNA = subset_mRNA_data), metadata = metadata, target_factor = "targets", corr_factor = "tissue")))

  # CASE 5: invalid data argument
  expect_error(RunDEAnalysis(list_data = subset_mRNA_data, metadata = metadata, corr_factor = "tissue", dim.reduct = "PCA"))

  # CASE 6: metadata missing
  expect_error(RunDEAnalysis(list_data = list(mRNA = subset_mRNA_data), corr_factor = "tissue"))
})

test_that("Exploratory approach combinations (dim.reduct and heatmap)", {

  data(subset_mRNA_data, metadata)

  # CASE 1: heatmap = TRUE, PCA
  expect_no_error(suppressWarnings(RunDEAnalysis(list_data = list(mRNA = subset_mRNA_data), metadata = metadata,
                                                 corr_factor = "tissue", heatmap = TRUE, dim.reduct = "PCA")))

  # CASE 2: heatmap = FALSE, UMAP
  expect_no_error(suppressWarnings(RunDEAnalysis(list_data = list(mRNA = subset_mRNA_data), metadata = metadata,
                                                 corr_factor = "tissue", heatmap = FALSE, dim.reduct = "UMAP", color_factor = "tissue")))

  # CASE 3: heatmap = TRUE, tSNE
  expect_no_error(suppressWarnings(RunDEAnalysis(list_data = list(mRNA = subset_mRNA_data), metadata = metadata,
                                                 corr_factor = "tissue", heatmap = TRUE, dim.reduct = "tSNE", color_factor = "tissue")))

  # CASE 4: heatmap = FALSE, no dim.reduct
  expect_no_error(suppressWarnings(RunDEAnalysis(list_data = list(mRNA = subset_mRNA_data), metadata = metadata,
                                                 corr_factor = "tissue", heatmap = FALSE, dim.reduct = NULL)))
})

test_that("PCA number of components validation", {

  data(subset_mRNA_data, metadata)

  # CASE 1: valid numPC
  expect_no_error(suppressWarnings(RunDEAnalysis(list_data = list(mRNA = subset_mRNA_data), metadata = metadata,
                                                 corr_factor = "tissue", dim.reduct = "PCA", numPC = 3, color_factor = "tissue")))

  # CASE 2: numPC = 0, expect message about invalid numPC
  expect_message(suppressWarnings(RunDEAnalysis(list_data = list(mRNA = subset_mRNA_data), metadata = metadata, corr_factor = "tissue",
                                                dim.reduct = "PCA", numPC = 0, color_factor = "tissue", verbose = FALSE)), "Invalid `numPC` given")
})

test_that("dim.reduct parameter validation", {

  data(subset_mRNA_data, metadata)

  # CASE 1: invalid dim.reduct should error
  expect_error(RunDEAnalysis(list_data = list(mRNA = subset_mRNA_data), metadata = metadata, corr_factor = "tissue", dim.reduct = "INVALID"))

  # CASE 2: valid dim.reduct values
  for (method in c("PCA", "UMAP", "tSNE")) {
    expect_no_error(suppressWarnings(RunDEAnalysis(list_data = list(mRNA = subset_mRNA_data),
                                                   metadata = metadata, corr_factor = "tissue", dim.reduct = method, color_factor = "tissue")))
  }
})
