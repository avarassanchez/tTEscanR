data("default_tTEscanR_mRNA_data", package = "tTEscanR")
data("default_tTEscanR_tRNA_data", package = "tTEscanR")

test_that("The tTEscanR object is correctly generated", {
    cell_types <- colnames(default_tTEscanR_mRNA_data)

    ## CASE 1: Valid creation formats (single matrix, lists, default assay)
    tTEobject <- createObject(counts = default_tTEscanR_mRNA_data, assay = "mRNA")
    expect_s4_class(object = tTEobject, class = "MultiAssayExperiment")
    expect_equal(SummarizedExperiment::assay(tTEobject, "mRNA"), as.matrix(default_tTEscanR_mRNA_data))

    expect_no_error(createObject(counts = default_tTEscanR_mRNA_data, assay = "mRNA"))
    expect_no_error(createObject(counts = list(mRNA = default_tTEscanR_mRNA_data)))
    expect_no_error(createObject(counts = list(default_tTEscanR_mRNA_data), assay = list("mRNA")))
    expect_no_error(createObject(counts = default_tTEscanR_mRNA_data)) # default assay

    ## CASE 2: Metadata integration
    expect_no_error(createObject(counts = default_tTEscanR_mRNA_data, assay = "mRNA", params = list("cell.types" = cell_types)))
    expect_no_error(createObject(counts = list(mRNA = default_tTEscanR_mRNA_data), params = list("cell.types" = cell_types)))
    expect_no_error(createObject(counts = list(default_tTEscanR_mRNA_data), assay = list("mRNA"), params = list("cell.types" = cell_types)))

    ## CASE 3: Constructor parameter and input format errors
    expect_error(createObject(), "argument \"counts\" is missing, with no default")
    expect_error(createObject(counts = matrix(1:4, nrow = 2), assay = "mRNA"), "'counts' must have column names when 'meta_data' is NULL.")
    expect_error(createObject(counts = c(1, 2, 3, 4, 5), assay = "mRNA"))
    expect_error(createObject(counts = c(default_tTEscanR_mRNA_data), assay = "mRNA"))
    expect_error(createObject(counts = default_tTEscanR_mRNA_data, assay = mRNA))
    expect_error(createObject(counts = default_tTEscanR_mRNA_data, assay = c("mRNA", "mRNA_2")))
    expect_error(createObject(counts = list(default_tTEscanR_mRNA_data, default_tTEscanR_mRNA_data), assay = list("mRNA", "mRNA"), params = list("cell.types" = cell_types)))
})

test_that("updateObject adds new counts and metadata correctly", {
    initial_counts <- matrix(1:4, nrow = 2, dimnames = list(c("gene1", "gene2"), c("s1", "s2")))
    obj <- createObject(counts = initial_counts, assay = "mRNA")

    ## CASE 1: Invalid assay name & missing overwrite flag checks
    new_counts <- matrix(5:8, nrow = 2, dimnames = list(c("gene1", "gene2"), c("s1", "s2")))
    expect_error(updateObject(object = obj, counts = new_counts, assay = "protein"), "Invalid 'assay' name")
    expect_error(updateObject(object = obj, counts = new_counts, assay = "mRNA"), "overwrite = TRUE")

    ## CASE 2: Valid assay overwrite and metadata update
    expect_no_error(updated_obj <- updateObject(object = obj, counts = new_counts, assay = "mRNA", overwrite = TRUE))
    new_metadata <- c(20, 30)
    updated_obj <- updateObject(object = updated_obj, params = list("patient_age" = new_metadata))

    ## CASE 3: Assert integrity of updated assays and metadata
    expect_true("mRNA" %in% names(updated_obj))
    expect_equal(SummarizedExperiment::assay(updated_obj, "mRNA"), new_counts)
    expect_true("patient_age" %in% names(S4Vectors::metadata(updated_obj)))
    expect_equal(S4Vectors::metadata(updated_obj)$patient_age, new_metadata)
})

test_that("updateObject handles overwrite correctly", {
    toy_counts <- matrix(1:4, nrow = 2, dimnames = list(NULL, c("V1", "V2")))
    obj <- createObject(counts = toy_counts, assay = "mRNA")

    ## CASE 1: Attempt overwrite without setting overwrite = TRUE
    expect_error(updateObject(object = obj, counts = toy_counts, assay = "mRNA"), "overwrite = TRUE")

    ## CASE 2: Attempt overwrite with unaligned column names
    new_RNA_counts <- matrix(9:12, nrow = 2)
    expect_error(updateObject(object = obj, counts = new_RNA_counts, assay = "mRNA", overwrite = TRUE), "colData rownames and ExperimentList colnames are empty")

    ## CASE 3: Successful overwrite with aligned column names
    colnames(new_RNA_counts) <- c("V3", "V4")
    updated_obj <- updateObject(object = obj, counts = new_RNA_counts, assay = "mRNA", overwrite = TRUE)
    expect_equal(SummarizedExperiment::assay(updated_obj, "mRNA"), new_RNA_counts)
})

test_that("updateObject handles named and unnamed metadata lists", {
    initial_counts <- matrix(1:4, nrow = 2, dimnames = list(c("gene1", "gene2"), c("s1", "s2")))
    obj <- createObject(counts = initial_counts, assay = "mRNA")

    ## CASE 1: Named metadata update
    named_meta <- list(age = c(25, 30), height = c(170, 180))
    expect_error(updateObject(object = obj, meta_data = named_meta), "The 'meta_data' requires matching row names to the columns in 'counts'")
    updated_obj_named <- updateObject(object = obj, params = named_meta)
    expect_equal(S4Vectors::metadata(updated_obj_named), named_meta)

    ## CASE 2: Unnamed metadata list rejection
    unnamed_meta <- list(c("m", "f"), c("urban", "rural"))
    expect_error(updateObject(object = obj, params = unnamed_meta), "'params' must be a named list.")
})

test_that("The tTEscanR object updates across diverse inputs and parameters", {
    tTEobject <- createObject(counts = default_tTEscanR_mRNA_data, assay = "mRNA")
    cell_types <- colnames(default_tTEscanR_mRNA_data)

    ## CASE 1: Proper update execution paths
    expect_no_error(updateObject(object = tTEobject, counts = default_tTEscanR_tRNA_data, assay = "tRNA"))
    expect_no_error(updateObject(object = tTEobject, params = list("cell.types" = cell_types)))
    expect_no_error(updateObject(object = tTEobject, counts = default_tTEscanR_tRNA_data, assay = "tRNA", params = list("cell.types" = cell_types)))

    ## CASE 2: Overwrite behavior on populated objects
    tTEobject <- updateObject(object = tTEobject, counts = default_tTEscanR_tRNA_data, assay = "tRNA")
    expect_s4_class(object = tTEobject, class = "MultiAssayExperiment")
    expect_equal(SummarizedExperiment::assay(tTEobject, "tRNA"), default_tTEscanR_tRNA_data)
    expect_error(updateObject(object = tTEobject, counts = default_tTEscanR_tRNA_data, assay = "tRNA"))
    expect_no_error(updateObject(object = tTEobject, counts = default_tTEscanR_tRNA_data, assay = "tRNA", overwrite = TRUE))
    expect_error(updateObject(object = tTEobject, counts = list(default_tTEscanR_mRNA_data, default_tTEscanR_tRNA_data), assay = list("mRNA", "tRNA")))
    expect_no_error(updateObject(object = tTEobject, counts = list(default_tTEscanR_mRNA_data, default_tTEscanR_tRNA_data), assay = list("mRNA", "tRNA"), overwrite = TRUE))

    ## CASE 3: Flexible list vs single input handling
    expect_no_error(updateObject(object = tTEobject, counts = list(default_tTEscanR_mRNA_data), assay = list("mRNA"), overwrite = TRUE))
    expect_no_error(updateObject(object = tTEobject, counts = list(default_tTEscanR_mRNA_data), assay = "mRNA", overwrite = TRUE))
    expect_no_error(updateObject(object = tTEobject, counts = default_tTEscanR_mRNA_data, assay = list("mRNA"), overwrite = TRUE))
    expect_error(updateObject(object = tTEobject, counts = c(default_tTEscanR_mRNA_data, default_tTEscanR_tRNA_data), assay = c("mRNA", "tRNA"), overwrite = TRUE))
    expect_error(updateObject(object = tTEobject, params = c("cell.types" = cell_types)), "'params' must be a named list")

    ## CASE 4: Missing parameters exceptions
    expect_error(updateObject(object = tTEobject))
    expect_error(updateObject(object = tTEobject, assay = "tRNA"))
    expect_error(updateObject(object = tTEobject, meta.data.ids = "cell_types"))
    expect_error(updateObject(object = tTEobject, counts = default_tTEscanR_tRNA_data))
    expect_error(updateObject(object = tTEobject, meta.data = default_tTEscanR_tRNA_data))
    expect_error(updateObject(object = tTEobject, meta.data = cell_types))
})
