test_that("The differential expression analysis is correctly executed (individual)", {

  data(subset_mRNA_data, metadata)

  # CASE 1 - controls
  expect_no_error(suppressWarnings(DESeq2runner(data = subset_mRNA_data, metadata = metadata, corr_factor = "tissue", heatmap = TRUE, PCA = TRUE, color_factor = "tissue")))
  expect_no_error(suppressWarnings(DESeq2runner(data = subset_mRNA_data, metadata = metadata, corr_factor = "tissue", numPC = 3, color_factor = "tissue", labels = TRUE, label_factor = "cell.type")))
  # CASE 2 - missing required parameters
  expect_no_error(suppressWarnings(DESeq2runner(data = subset_mRNA_data, metadata = metadata, corr_factor = "tissue", heatmap = FALSE, PCA = FALSE))) # color_factor missing
  expect_error(suppressWarnings(DESeq2runner(data = subset_mRNA_data, metadata = metadata))) # color_factor missing
  expect_error(suppressWarnings(DESeq2runner(data = subset_mRNA_data, corr_factor = "tissue"))) # metadata not included

  # CASE 3 - usage of the targets parameter
  targets_endothelial <- data.frame(search = c("endothelial"), class = c("endothelial"))
  expect_no_error(suppressWarnings(DESeq2runner(data = subset_mRNA_data, metadata = metadata, corr_factor = "tissue", targets = targets_endothelial, target_factor = "cell.type", fc_threshold = 1, pval_threshold = 0.05)))
  expect_error(suppressWarnings(DESeq2runner(data = subset_mRNA_data, metadata = metadata, corr_factor = "tissue", targets = targets_endothelial, fc_threshold = 1, pval_threshold = 0.05)))

  targets_neuron <- data.frame(search = c("neuron", "ENS neurons"), class = c("neuron", "other"))
  expect_no_error(suppressWarnings(DESeq2runner(data = subset_mRNA_data, metadata = metadata, corr_factor = "tissue", targets = targets_neuron, target_factor = "cell.type", fc_threshold = 1, pval_threshold = 0.05)))
  expect_error(suppressWarnings(DESeq2runner(data = subset_mRNA_data, metadata = metadata, corr_factor = "tissue", targets = targets_neuron, fc_threshold = 1, pval_threshold = 0.05)))
})

test_that("The differential expression analysis is correctly executed (multiple)", {

  data(subset_mRNA_data, metadata)

  # CASE 1 - controls
  expect_no_error(suppressWarnings(ExecuteDESeq2runner(list_data = list(mRNA = subset_mRNA_data), metadata = metadata, corr_factor = "tissue")))
  expect_no_error(suppressWarnings(ExecuteDESeq2runner(list_data = list(mRNA = subset_mRNA_data), metadata = metadata, corr_factor = "tissue", color_factor = "tissue")))
  expect_no_error(suppressWarnings(ExecuteDESeq2runner(list_data = list(mRNA = subset_mRNA_data), metadata = metadata, corr_factor = "tissue", color_factor = "tissue", labels = TRUE)))

  # CASE 2
  expect_no_error(suppressWarnings(ExecuteDESeq2runner(list_data = list(mRNA = subset_mRNA_data), metadata = metadata, corr_factor = "tissue"))) # color_factor missing
  expect_error(suppressWarnings(ExecuteDESeq2runner(list_data = list(mRNA = subset_mRNA_data), metadata = metadata))) # corr_factor missing
  expect_error(suppressWarnings(ExecuteDESeq2runner(list_data = list(subset_mRNA_data), metadata = metadata)))

  # CASE 3
  targets_endothelial <- data.frame(search = c("endothelial"), class = c("endothelial"))
  expect_no_error(suppressWarnings(ExecuteDESeq2runner(list_data = list(mRNA = subset_mRNA_data), metadata = metadata, corr_factor = "tissue",
                                                       targets = targets_endothelial, target_factor = "cell.type", fc_threshold = 1, pval_threshold = 0.05)))
  expect_error(suppressWarnings(ExecuteDESeq2runner(list_data = list(mRNA = subset_mRNA_data), metadata = metadata, corr_factor = "tissue",
                                                    targets = targets_endothelial, fc_threshold = 1, pval_threshold = 0.05)))

  # CASE 4
  targets_neuron <- data.frame(search = c("neuron", "ENS neurons"), class = c("neuron", "other"))
  results <- expect_no_error(suppressWarnings(ExecuteDESeq2runner(list_data = list(mRNA = subset_mRNA_data), metadata = metadata, corr_factor = "tissue",
                                                                  targets = targets_neuron, target_factor = "cell.type", fc_threshold = 1, pval_threshold = 0.05)))
  results <- expect_error(suppressWarnings(ExecuteDESeq2runner(list_data = list(mRNA = subset_mRNA_data), metadata = metadata, corr_factor = "tissue",
                                                               targets = targets_neuron, fc_threshold = 1, pval_threshold = 0.05)))

  # CASE 5 - incorrect data argument
  expect_error(ExecuteDESeq2runner(list_data = subset_mRNA_data, metadata = metadata, corr_factor = "tissue", heatmap = TRUE, PCA = TRUE, numPC = 3))
  expect_error(ExecuteDESeq2runner(data = subset_mRNA_data, metadata = metadata, corr_factor = "tissue", heatmap = TRUE, PCA = TRUE, numPC = 3))

  # CASE 6
  expect_error(ExecuteDESeq2runner(list_data = list(mRNA = subset_mRNA_data), corr_factor = "tissue")) # metadata not included

  # CASE 7
  expect_no_error(suppressWarnings(ExecuteDESeq2runner(list_data = list(subset_mRNA_data), metadata = metadata, corr_factor = "tissue", heatmap = TRUE,
                                                       PCA = TRUE, numPC = 3, color_factor = "tissue")))
})

test_that("The differential expression analysis is correctly executed (combinations PCA and heatmap)", {

  data(subset_mRNA_data, metadata)

  # CONTROLS
  expect_no_error(suppressWarnings(ExecuteDESeq2runner(list_data = list(mRNA = subset_mRNA_data), metadata = metadata, corr_factor = "tissue")))
  expect_no_error(suppressWarnings(DESeq2runner(data = subset_mRNA_data, metadata = metadata, corr_factor = "tissue", verbose = FALSE, targets = NULL)))

  # CASE 1 - heatmap = TRUE and PCA = TRUE
  expect_no_error(suppressWarnings(ExecuteDESeq2runner(list_data = list(mRNA = subset_mRNA_data), metadata = metadata, corr_factor = "tissue", heatmap = TRUE, PCA = TRUE)))
  expect_no_error(suppressWarnings(DESeq2runner(data = subset_mRNA_data, metadata = metadata, corr_factor = "tissue", heatmap = TRUE, PCA = TRUE)))

  # CASE 2 - heatmap = FALSE and PCA = FALSE
  expect_no_error(suppressWarnings(ExecuteDESeq2runner(list_data = list(mRNA = subset_mRNA_data), metadata = metadata, corr_factor = "tissue", heatmap = FALSE, PCA = FALSE)))
  expect_no_error(suppressWarnings(DESeq2runner(data = subset_mRNA_data, metadata = metadata, corr_factor = "tissue", heatmap = FALSE, PCA = FALSE)))

  # CASE 3 - heatmap = FALSE and PCA = TRUE
  expect_no_error(suppressWarnings(ExecuteDESeq2runner(list_data = list(mRNA = subset_mRNA_data), metadata = metadata, corr_factor = "tissue", heatmap = FALSE, PCA = TRUE)))
  expect_no_error(suppressWarnings(DESeq2runner(data = subset_mRNA_data, metadata = metadata, corr_factor = "tissue", heatmap = FALSE, PCA = TRUE)))

  # CASE 4 - heatmap = TRUE and PCA = FALSE
  expect_no_error(suppressWarnings(ExecuteDESeq2runner(list_data = list(mRNA = subset_mRNA_data), metadata = metadata, corr_factor = "tissue", heatmap = TRUE, PCA = FALSE)))
  expect_no_error(suppressWarnings(DESeq2runner(data = subset_mRNA_data, metadata = metadata, corr_factor = "tissue", heatmap = TRUE, PCA = FALSE)))
})

test_that("The differential expression analysis is correctly executed (evaluating the PCs)", {

  data(subset_mRNA_data, metadata)

  # CONTROLS
  expect_no_error(suppressWarnings(ExecuteDESeq2runner(list_data = list(mRNA = subset_mRNA_data), metadata = metadata, corr_factor = "tissue",)))
  expect_no_error(suppressWarnings(DESeq2runner(data = subset_mRNA_data, metadata = metadata, corr_factor = "tissue")))

  # CASE 1 - default numPC = 2
  expect_no_error(suppressWarnings(ExecuteDESeq2runner(list_data = list(mRNA = subset_mRNA_data), metadata = metadata, heatmap = FALSE, PCA = TRUE, corr_factor = "tissue", numPC = 2)))
  expect_no_error(suppressWarnings(DESeq2runner(data = subset_mRNA_data, metadata = metadata, corr_factor = "tissue", heatmap = FALSE, PCA = TRUE, numPC = 2)))

  # CASE 2 - forcing numPCs = 0
  expect_message(suppressWarnings(ExecuteDESeq2runner(list_data = list(mRNA = subset_mRNA_data), metadata = metadata, heatmap = FALSE, PCA = TRUE,
                                                      corr_factor = "tissue", numPC = 0, verbose = FALSE)), "Invalid `numPC` given")
  expect_message(suppressWarnings(DESeq2runner(data = subset_mRNA_data, metadata = metadata, corr_factor = "tissue", heatmap = FALSE, PCA = TRUE,
                                               numPC = 0, verbose = FALSE)), "Invalid `numPC` given")
})
