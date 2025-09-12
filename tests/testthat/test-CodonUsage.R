test_that("The different corr_method possibilities - inside ComputeCodonUsage()", {

  data(subset_mRNA_data, metadata)

  tTEobject <- Create_tTEscanR_Object(counts = subset_mRNA_data, assay = "mRNA", meta.data = list(metadata, "tissue"), meta.data.ids = list("ConditionsLabels", "CorrectionFactor"))

  # CASE 1: no error - default parameters without computing additional metrics
  expect_no_error(suppressWarnings(ComputeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", filter = "canonical",
                                                     additional.metrics = FALSE, verbose = FALSE)))
  expect_no_error(suppressWarnings(ComputeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", filter = "canonical",
                                                     additional.metrics = FALSE, corr_method = "spearman", verbose = FALSE))) # corr_method will be ignored as it is only required if additional.metrics = TRUE

  # CASE 2: no error - default parameters and computing additional metrics
  expect_no_error(suppressWarnings(ComputeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", filter = "canonical",
                                                     additional.metrics = TRUE, verbose = FALSE))) # the default corr_method will be used
  expect_no_error(suppressWarnings(ComputeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", filter = "canonical",
                                                     additional.metrics = TRUE, corr_method = "spearman", verbose = FALSE)))

  # CASE 3: no error - trying other corr_method parameters
  expect_no_error(suppressWarnings(ComputeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", filter = "canonical",
                                                     additional.metrics = TRUE, corr_method = "pearson", verbose = FALSE)))
  expect_no_error(suppressWarnings(ComputeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", filter = "canonical",
                                                     additional.metrics = TRUE, corr_method = "kendall", verbose = FALSE)))

  # CASE 4: error - trying non-existent corr_method parameters
  expect_error(suppressWarnings(ComputeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", filter = "canonical",
                                                  additional.metrics = TRUE, corr_method = "spearmans", verbose = FALSE)))
  expect_error(suppressWarnings(ComputeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", filter = "canonical",
                                                  additional.metrics = TRUE, corr_method = "pearsons", verbose = FALSE)))
  expect_error(suppressWarnings(ComputeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", filter = "canonical",
                                                  additional.metrics = TRUE, corr_method = "kendalls", verbose = FALSE)))

})

test_that("The different corr_method possibilities - inside ComputeAdditionalMetrics_CodonUsage()", {

  data(subset_mRNA_data, metadata)

  tTEobject <- Create_tTEscanR_Object(counts = subset_mRNA_data, assay = "mRNA", meta.data = metadata, meta.data.ids = "ConditionsLabels")
  tTEobject <- ComputeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", filter = "canonical", additional.metrics = FALSE)

  # CASE 1: no error - default
  expect_no_error(suppressWarnings(ComputeAdditionalMetrics_CodonUsage(codon_usage = tTEobject@assays$CodonUsage, codon_freq = codon_freq_table_hg38_canonical,
                                                                       metadata = tTEobject@meta.data$ConditionsLabels, corr_method = "spearman", corr_factor = "tissue", verbose = FALSE)))
  expect_no_error(suppressWarnings(ComputeAdditionalMetrics_CodonUsage(codon_usage = tTEobject@assays$CodonUsage, codon_freq = codon_freq_table_hg38_length,
                                                                       metadata = tTEobject@meta.data$ConditionsLabels, corr_method = "spearman", corr_factor = "tissue", verbose = FALSE)))

  # CASE 2: no error - trying other corr_method parameters
  expect_no_error(suppressWarnings(ComputeAdditionalMetrics_CodonUsage(codon_usage = tTEobject@assays$CodonUsage, codon_freq = codon_freq_table_hg38_canonical,
                                                                       metadata = tTEobject@meta.data$ConditionsLabels, corr_method = "pearson", corr_factor = "tissue", verbose = FALSE)))
  expect_no_error(suppressWarnings(ComputeAdditionalMetrics_CodonUsage(codon_usage = tTEobject@assays$CodonUsage, codon_freq = codon_freq_table_hg38_length,
                                                                       metadata = tTEobject@meta.data$ConditionsLabels, corr_method = "pearson", corr_factor = "tissue", verbose = FALSE)))

  expect_no_error(suppressWarnings(ComputeAdditionalMetrics_CodonUsage(codon_usage = tTEobject@assays$CodonUsage, codon_freq = codon_freq_table_hg38_canonical, metadata = tTEobject@meta.data$ConditionsLabels, corr_method = "kendall", corr_factor = "tissue", verbose = FALSE)))
  expect_no_error(suppressWarnings(ComputeAdditionalMetrics_CodonUsage(codon_usage = tTEobject@assays$CodonUsage, codon_freq = codon_freq_table_hg38_length,
                                                                       metadata = tTEobject@meta.data$ConditionsLabels, corr_method = "kendall", corr_factor = "tissue", verbose = FALSE)))

  # CASE 3: error - trying non-existent corr_method parameters
  expect_error(suppressWarnings(ComputeAdditionalMetrics_CodonUsage(codon_usage = tTEobject@assays$CodonUsage, codon_freq = codon_freq_table_hg38_canonical,
                                                                    metadata = tTEobject@meta.data$ConditionsLabels, corr_method = "spearmans", verbose = FALSE)))
  expect_error(suppressWarnings(ComputeAdditionalMetrics_CodonUsage(codon_usage = tTEobject@assays$CodonUsage, codon_freq = codon_freq_table_hg38_canonical,
                                                                    metadata = tTEobject@meta.data$ConditionsLabels, corr_method = "pearsons", verbose = FALSE)))
  expect_error(suppressWarnings(ComputeAdditionalMetrics_CodonUsage(codon_usage = tTEobject@assays$CodonUsage, codon_freq = codon_freq_table_hg38_canonical,
                                                                    metadata = tTEobject@meta.data$ConditionsLabels, corr_method = "kendalls", verbose = FALSE)))
  # CASE 4: error - no specification of corr_factor
  expect_error(suppressWarnings(ComputeAdditionalMetrics_CodonUsage(codon_usage = tTEobject@assays$CodonUsage, codon_freq = codon_freq_table_hg38_canonical,
                                                                    metadata = tTEobject@meta.data$ConditionsLabels, corr_method = "spearman", verbose = FALSE)))
  expect_error(suppressWarnings(ComputeAdditionalMetrics_CodonUsage(codon_usage = tTEobject@assays$CodonUsage, codon_freq = codon_freq_table_hg38_length,
                                                                    metadata = tTEobject@meta.data$ConditionsLabels, verbose = FALSE)))
})

test_that("The different corr_method possibilities - inside ComputeCorrelationBackground()", {

  data(subset_mRNA_data, metadata)

  tTEobject <- Create_tTEscanR_Object(counts = subset_mRNA_data, assay = "mRNA", meta.data = list(metadata, "tissue"), meta.data.ids = list("ConditionsLabels", "CorrectionFactor"))
  tTEobject <- ComputeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", filter = "canonical", additional.metrics = FALSE)

  mean_codon_usage <- suppressWarnings(ComputeMeanUsage(data = tTEobject, assay = "CodonUsage", mode = "raw"))
  exonic_background <- suppressWarnings(ComputeCodonExonicBackground(data = tTEobject@assays$CodonUsage))

  # CASE 1: no error - default
  expect_no_error(suppressWarnings(ComputeCorrelationBackground(mean = mean_codon_usage, background = exonic_background, corr_method = "spearman")))
  expect_no_error(suppressWarnings(ComputeCorrelationBackground(mean = mean_codon_usage, background = exonic_background)))

  # CASE 2: no error - trying other corr_method parameters
  expect_no_error(suppressWarnings(ComputeCorrelationBackground(mean = mean_codon_usage, background = exonic_background, corr_method = "pearson")))
  expect_no_error(suppressWarnings(ComputeCorrelationBackground(mean = mean_codon_usage, background = exonic_background, corr_method = "kendall")))

  # CASE 3: error - trying non-existent corr_method parameters
  expect_error(suppressWarnings(ComputeCorrelationBackground(mean = mean_codon_usage, background = exonic_background, corr_method = "spearmans")))
  expect_error(suppressWarnings(ComputeCorrelationBackground(mean = mean_codon_usage, background = exonic_background, corr_method = "pearsons")))
  expect_error(suppressWarnings(ComputeCorrelationBackground(mean = mean_codon_usage, background = exonic_background, corr_method = "kendalls")))
})

test_that("The different parameters to perform the mean codon usage - ComputeMeanUsage()", {

  data(subset_mRNA_data, metadata)

  tTEobject <- Create_tTEscanR_Object(counts = subset_mRNA_data, assay = "mRNA", meta.data = list(metadata, "tissue"), meta.data.ids = list("ConditionsLabels", "CorrectionFactor"))
  tTEobject <- ComputeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", filter = "canonical", additional.metrics = FALSE)

  # CASE 1: no error - using object
  expect_no_error(suppressWarnings(ComputeMeanUsage(data = tTEobject, assay = "CodonUsage", mode = "raw")))
  expect_no_error(suppressWarnings(ComputeMeanUsage(data = tTEobject, assay = "CodonUsage", mode = "raw", metadata = tTEobject@meta.data$ConditionsLabels, corr_factor = "tissue"))) # metadata and corr_factor will be ignored

  # CASE 2: error - using object wrong assay
  expect_error(ComputeMeanUsage(data = tTEobject, assay = CodonUsage, mode = "raw", metadata = tTEobject@meta.data$ConditionsLabels, corr_factor = "tissue")) # the assay needs to be a string
  expect_error(ComputeMeanUsage(data = tTEobject, assay = "CodonUse", mode = "raw", metadata = tTEobject@meta.data$ConditionsLabels)) # non-existent assay
  expect_error(ComputeMeanUsage(data = tTEobject, mode = "raw")) # assay not specified

  # CASE 3: no error - using a dataset
  expect_no_error(suppressWarnings(ComputeMeanUsage(data = tTEobject@assays$CodonUsage, mode = "raw", metadata = tTEobject@meta.data$ConditionsLabels, corr_factor = "tissue")))

  # CASE 4: error - no specifying the corr_factor or metadata
  expect_error(suppressWarnings(ComputeMeanUsage(data = tTEobject@assays$CodonUsage, mode = "raw", metadata = tTEobject@meta.data$ConditionsLabels)))
  expect_error(suppressWarnings(ComputeMeanUsage(data = tTEobject@assays$CodonUsage, mode = "raw", corr_factor = "tissue")))
})

test_that("The mean usage across different features - ComputeMeanUsage()", {

  data(subset_mRNA_data, subset_tRNA_data, metadata)

  tTEobject <- Create_tTEscanR_Object(counts = list(subset_mRNA_data, subset_tRNA_data), assay = list("mRNA", "tRNA"), meta.data = list(metadata, "tissue"), meta.data.ids = list("ConditionsLabels", "CorrectionFactor"))
  tTEobject <- ComputeCodonUsage(object = tTEobject, codon_freq = NULL, species = "hg38", filter = "canonical", additional.metrics = FALSE)
  tTEobject <- ComputeAnticodonUsage(object = tTEobject)
  tTEobject <- ComputeAAUsage(object = tTEobject, level = "both")

  # CASE 1: no error - CODON
  expect_no_error(suppressWarnings(ComputeMeanUsage(data = tTEobject, assay = "CodonUsage", mode = "raw")))
  expect_no_error(suppressWarnings(ComputeMeanUsage(data = tTEobject@assays$CodonUsage, mode = "raw", metadata = tTEobject@meta.data$ConditionsLabels, corr_factor = "tissue")))

  # CASE 2: no error - ANTICODON
  expect_no_error(suppressWarnings(ComputeMeanUsage(data = tTEobject, assay = "AnticodonUsage", mode = "raw")))
  expect_no_error(suppressWarnings(ComputeMeanUsage(data = tTEobject@assays$AnticodonUsage, mode = "raw", metadata = tTEobject@meta.data$ConditionsLabels, corr_factor = "tissue")))

  # CASE 3: no error - AA DEMAND
  expect_no_error(suppressWarnings(ComputeMeanUsage(data = tTEobject, assay = "AADemand", mode = "raw")))
  expect_no_error(suppressWarnings(ComputeMeanUsage(data = tTEobject@assays$AADemand, mode = "raw", metadata = tTEobject@meta.data$ConditionsLabels, corr_factor = "tissue")))

  # CASE 4: no error - AA SUPPLY
  expect_no_error(suppressWarnings(ComputeMeanUsage(data = tTEobject, assay = "AASupply", mode = "raw")))
  expect_no_error(suppressWarnings(ComputeMeanUsage(data = tTEobject@assays$AASupply, mode = "raw", metadata = tTEobject@meta.data$ConditionsLabels, corr_factor = "tissue")))
})
