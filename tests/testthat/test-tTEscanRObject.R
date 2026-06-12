test_that("UpdateObject adds new counts and metadata correctly", {
    initial_counts <- matrix(1:4, nrow = 2, dimnames = list(c("gene1", "gene2"), c("s1", "s2")))
    obj <- CreateObject(counts = initial_counts, assay = "mRNA")

    # CASE 1: Add new counts
    new_counts <- matrix(5:8, nrow = 2, dimnames = list(c("gene1", "gene2"), c("s1", "s2")))

    expect_error(updated_obj <- UpdateObject(object = obj, counts = new_counts, assay = "protein"), "Invalid 'assay' names detected:")
    expect_error(updated_obj <- UpdateObject(object = obj, counts = new_counts, assay = "mRNA"), "overwrite = TRUE")
    expect_no_error(updated_obj <- UpdateObject(object = obj, counts = new_counts, assay = "mRNA", overwrite = TRUE))
    updated_obj <- UpdateObject(object = obj, counts = new_counts, assay = "mRNA", overwrite = TRUE)

    # CASE 2: add new metadata
    new_metadata <- c(20, 30)
    updated_obj <- UpdateObject(object = updated_obj, meta.data = new_metadata, meta.data.ids = "patient_age")

    # Check if assays and metadata were added correctly
    expect_true("mRNA" %in% names(updated_obj@assays))
    expect_equal(updated_obj@assays$mRNA, new_counts)
    expect_true("patient_age" %in% names(updated_obj@meta.data))
    expect_equal(updated_obj@meta.data$patient_age, new_metadata)
})

test_that("UpdateObject handles overwrite correctly", {
    obj <- CreateObject(counts = matrix(1:4, nrow = 2), assay = "mRNA")

    # CASE 1: Attempt to overwrite without setting overwrite = TRUE (should fail)
    expect_error(UpdateObject(object = obj, counts = matrix(5:8, nrow = 2), assay = "mRNA"), "overwrite = TRUE")

    # CASE 2: Overwrite with overwrite = TRUE (should succeed)
    new_RNA_counts <- matrix(9:12, nrow = 2)
    updated_obj <- UpdateObject(object = obj, counts = new_RNA_counts, assay = "mRNA", overwrite = TRUE)

    expect_equal(updated_obj@assays$mRNA, new_RNA_counts)
})

test_that("UpdateObject handles named and unnamed metadata lists", {
    expect_error(obj <- CreateObject(), "argument \"counts\" is missing, with no default")
    initial_counts <- matrix(1:4, nrow = 2, dimnames = list(c("gene1", "gene2"), c("s1", "s2")))
    obj <- CreateObject(counts = initial_counts, assay = "mRNA")

    # CASE 1: Add metadata with a named list
    named_meta <- list(age = c(25, 30), height = c(170, 180))
    updated_obj_named <- UpdateObject(object = obj, meta.data = named_meta)

    expect_equal(updated_obj_named@meta.data, named_meta)

    # CASE 2: Add metadata with an unnamed list and meta.data.ids
    unnamed_meta <- list(c("m", "f"), c("urban", "rural"))
    updated_obj_unnamed <- UpdateObject(object = obj, meta.data = unnamed_meta, meta.data.ids = c("gender", "location"))

    expected_unnamed_meta <- list(gender = c("m", "f"), location = c("urban", "rural"))
    expect_equal(updated_obj_unnamed@meta.data, expected_unnamed_meta)
})

test_that("UpdateObject handles overwrite correctly", {
    obj <- new("tTEscanR_Object", assays = list(mRNA = matrix(1:4, nrow = 2)))

    # Test Case 2: Attempt to overwrite without setting overwrite = TRUE (should fail)
    expect_error(UpdateObject(object = obj, counts = matrix(5:8, nrow = 2), assay = "mRNA"), "overwrite = TRUE")

    # Test Case 3: Overwrite with overwrite = TRUE (should succeed)
    new_RNA_counts <- matrix(9:12, nrow = 2)
    updated_obj <- UpdateObject(object = obj, counts = new_RNA_counts, assay = "mRNA", overwrite = TRUE)

    expect_equal(updated_obj@assays$mRNA, new_RNA_counts)
})

test_that("The tTEscanR object is correctly generated", {
    tTEobject <- CreateObject(counts = mRNA_data_test, assay = "mRNA")
    cell_types <- colnames(mRNA_data_test)

    expect_s4_class(object = tTEobject, class = "tTEscanR_Object")
    expect_equal(tTEobject@assays$mRNA, mRNA_data_test)
    expect_no_error(IsInObject(tTEobject, "assays", "mRNA"))

    # CASE 1: no error - counts can be a list or a single input
    expect_no_error(CreateObject(counts = mRNA_data_test, assay = "mRNA"))
    expect_no_error(CreateObject(counts = list(mRNA = mRNA_data_test)))
    expect_no_error(CreateObject(counts = list(mRNA_data_test), assay = list("mRNA")))

    # CASE 2: error - counts has a wrong data format
    expect_error(CreateObject(counts = mRNA_data_test))
    expect_error(CreateObject(counts = c(1, 2, 3, 4, 5), assay = "mRNA"))
    expect_error(CreateObject(counts = c(mRNA_data_test), assay = "mRNA"))

    # CASE 3: error - wrong usage of the assay parameter
    expect_error(CreateObject(counts = mRNA_data_test, assay = mRNA))
    expect_error(CreateObject(counts = mRNA_data_test, assay = c("mRNA", "mRNA_2")))

    # CASE 4: no error - incorporation of metadata
    expect_no_error(CreateObject(counts = mRNA_data_test, assay = "mRNA", meta.data = cell_types, meta.data.ids = "cell.types"))
    expect_no_error(CreateObject(counts = list(mRNA = mRNA_data_test), meta.data = cell_types, meta.data.ids = "cell.types"))
    expect_no_error(CreateObject(counts = list(mRNA_data_test), assay = list("mRNA"), meta.data = cell_types, meta.data.ids = "cell.types"))
    expect_error(CreateObject(counts = mRNA_data_test, meta.data = list(cell_types), meta.data.ids = list("cell.types")))
    expect_no_error(CreateObject(counts = list(mRNA_data_test), assay = "mRNA", meta.data = cell_types, meta.data.ids = "cell.types"))

    # CASE 5: error - repeated ids
    expect_no_error(CreateObject(counts = list(mRNA = mRNA_data_test), assay = list("mRNA"), meta.data = cell_types, meta.data.ids = "cell.types")) # The name from the assay will be used
    expect_no_error(CreateObject(counts = list(mRNA_data_test), assay = list("mRNA"), meta.data = list(cell_types = cell_types), meta.data.ids = list("cell.types")))
    expect_error(CreateObject(counts = list(mRNA_data_test, mRNA_data_test), assay = list("mRNA", "mRNA"), meta.data = cell_types, meta.data.ids = "cell.types"))
})

test_that("The tTEscanR object is correctly updated", {
    tTEobject <- CreateObject(counts = mRNA_data_test, assay = "mRNA")
    cell_types <- colnames(mRNA_data_test)

    expect_s4_class(object = tTEobject, class = "tTEscanR_Object")
    expect_equal(tTEobject@assays$mRNA, mRNA_data_test)

    # CASE 0: check if individual sections are in the tTEscanR object
    expect_no_error(IsInObject(tTEobject, slot = "assays", section = "mRNA"))
    expect_error(IsInObject(tTEobject, slot = "assays", section = "tRNA"))

    # CASE 1: proper way to update the tTEscanR object
    expect_no_error(UpdateObject(object = tTEobject, counts = default_tTEscanR_tRNA_data, assay = "tRNA")) # just counts
    expect_no_error(UpdateObject(object = tTEobject, meta.data = cell_types, meta.data.ids = "cell_types")) # just metadata
    expect_no_error(UpdateObject(object = tTEobject, counts = default_tTEscanR_tRNA_data, assay = "tRNA", meta.data = cell_types, meta.data.ids = "cell_types")) # counts and metadata

    # CASE 2: use of overwrite
    tTEobject <- UpdateObject(object = tTEobject, counts = default_tTEscanR_tRNA_data, assay = "tRNA")
    expect_s4_class(object = tTEobject, class = "tTEscanR_Object")
    expect_equal(tTEobject@assays$tRNA, default_tTEscanR_tRNA_data)
    expect_no_error(IsInObject(tTEobject, slot = "assays", section = "mRNA"))
    expect_no_error(IsInObject(tTEobject, slot = "assays", section = "tRNA"))
    expect_error(UpdateObject(object = tTEobject, counts = default_tTEscanR_tRNA_data, assay = "tRNA")) # tRNA already in tTEobject
    expect_no_error(UpdateObject(object = tTEobject, counts = default_tTEscanR_tRNA_data, assay = "tRNA", overwrite = TRUE))
    expect_error(UpdateObject(object = tTEobject, counts = list(default_tTEscanR_mRNA_data, default_tTEscanR_tRNA_data), assay = list("mRNA", "tRNA")))
    expect_no_error(UpdateObject(object = tTEobject, counts = list(default_tTEscanR_mRNA_data, default_tTEscanR_tRNA_data), assay = list("mRNA", "tRNA"), overwrite = TRUE))

    # CASE 3: different input formats (list or single)
    expect_no_error(UpdateObject(object = tTEobject, counts = list(default_tTEscanR_mRNA_data), assay = list("mRNA"), overwrite = TRUE))
    expect_no_error(UpdateObject(object = tTEobject, counts = list(default_tTEscanR_mRNA_data), assay = "mRNA", overwrite = TRUE))
    expect_no_error(UpdateObject(object = tTEobject, counts = default_tTEscanR_mRNA_data, assay = list("mRNA"), overwrite = TRUE))
    expect_error(UpdateObject(object = tTEobject, counts = c(default_tTEscanR_mRNA_data, default_tTEscanR_tRNA_data), assay = c("mRNA", "tRNA"), overwrite = TRUE)) # using vectors

    cell_types_tRNA <- colnames(default_tTEscanR_tRNA_data)
    expect_no_error(UpdateObject(object = tTEobject, meta.data = list(cell_types), meta.data.ids = list("cell_types")))
    expect_no_error(UpdateObject(object = tTEobject, meta.data = list(cell_types, cell_types_tRNA), meta.data.ids = list("cell_types", "cell_types_tRNA")))
    expect_error(UpdateObject(object = tTEobject, meta.data = list(cell_types, cell_types_tRNA), meta.data.ids = list("cell_types", "cell_types"))) # same ids
    expect_error(UpdateObject(object = tTEobject, meta.data = c(cell_types, cell_types_tRNA), meta.data.ids = c("cell_types", "cell_types_tRNA"))) # using vectors

    # CASE 4: error - missing parameters
    expect_error(UpdateObject(object = tTEobject)) # nothing to add
    expect_error(UpdateObject(object = tTEobject, assay = "tRNA")) # nothing to add
    expect_error(UpdateObject(object = tTEobject, meta.data.ids = "cell_types")) # nothing to add
    expect_error(UpdateObject(object = tTEobject, counts = default_tTEscanR_tRNA_data)) # no ids
    expect_error(UpdateObject(object = tTEobject, meta.data = default_tTEscanR_tRNA_data)) # no ids
    expect_error(UpdateObject(object = tTEobject, meta.data = cell_types)) # no ids
})

test_that("The function to find sections in a tTEscanR object works", {
    tTEobject <- CreateObject(counts = default_tTEscanR_mRNA_data, assay = "mRNA")
    tTEobject <- UpdateObject(object = tTEobject, counts = default_tTEscanR_tRNA_data, assay = "tRNA")
    tTEobject <- UpdateObject(object = tTEobject, meta.data = colnames(default_tTEscanR_mRNA_data), meta.data.ids = "cell.types")

    expect_s4_class(object = tTEobject, class = "tTEscanR_Object")
    expect_equal(tTEobject@assays$mRNA, default_tTEscanR_mRNA_data)
    expect_equal(tTEobject@assays$tRNA, default_tTEscanR_tRNA_data)

    # CASE 1: no error - right slot and section
    expect_no_error(IsInObject(object = tTEobject, slot = "assays", section = "mRNA"))
    expect_no_error(IsInObject(object = tTEobject, slot = "assays", section = "tRNA"))
    expect_no_error(IsInObject(object = tTEobject, slot = "meta.data", section = "cell.types"))

    expect_no_error(IsInObject(object = tTEobject, slot = "assays", section = "mRNA", update_assay = FALSE, overwrite = FALSE))
    expect_no_error(IsInObject(object = tTEobject, slot = "assays", section = "tRNA", update_assay = FALSE, overwrite = FALSE))
    expect_no_error(IsInObject(object = tTEobject, slot = "meta.data", section = "cell.types", update_assay = FALSE, overwrite = FALSE))

    # CASE 2: error - wrong slot
    expect_error(IsInObject(object = tTEobject, slot = "meta.data", section = "mRNA"))
    expect_error(IsInObject(object = tTEobject, slot = "meta.data", section = "tRNA"))
    expect_error(IsInObject(object = tTEobject, slot = "assays", section = "cell.types"))

    # CASE 3: error - non available section
    expect_error(IsInObject(object = tTEobject, slot = "assays", section = mRNA))
    expect_error(IsInObject(object = tTEobject, slot = "assays", section = "CodonUsage"))

    # CASE 4: control of the slots
    expect_no_error(IsInObject(object = tTEobject, slot = "assay", section = "mRNA")) # with the implementation of match.arg() it is considered correctly as assays
    expect_error(IsInObject(object = tTEobject, slot = "counts", section = "tRNA"))
    expect_error(IsInObject(object = tTEobject, slot = "metadata", section = "cell.types"))
    expect_no_error(IsInObject(object = tTEobject, slot = "meta", section = "cell.types"))
})

test_that("The metadata is correctly updated", {
    tTEobject <- CreateObject(counts = mRNA_data_test, assay = "mRNA")

    # CASE 1: no error - single addition with and without lists
    expect_no_error(UpdateObject(object = tTEobject, meta.data = ENSG_gene_names_mRNA_data, meta.data.ids = "ENSG_gene_names_mRNA"))
    expect_no_error(UpdateObject(object = tTEobject, meta.data = list(ENSG_gene_names_mRNA_data), meta.data.ids = list("ENSG_gene_names_mRNA")))

    # CASE 2: no error - multiple addition
    expect_no_error(UpdateObject(object = tTEobject, meta.data = list(ENSG_gene_names_mRNA_data, short_gene_names_mRNA_data), meta.data.ids = list("ENSG_gene_names_mRNA", "short_gene_names_mRNA")))
    expect_no_error(UpdateObject(object = tTEobject, meta.data = list(ENSG = ENSG_gene_names_mRNA_data, short = short_gene_names_mRNA_data)))

    # CASE 3: error - inconsistent input format
    expect_no_error(UpdateObject(object = tTEobject, meta.data = list(short_gene_names_mRNA_data), meta.data.ids = "short_gene_names_mRNA"))
    expect_no_error(UpdateObject(object = tTEobject, meta.data = list(ENSG = ENSG_gene_names_mRNA_data), meta.data.ids = "ENSG_gene_names_mRNA_2"))

    # CASE 4: using the overwrite parameter
    tTEobject <- UpdateObject(object = tTEobject, meta.data = list(ENSG_gene_names_mRNA_data), meta.data.ids = list("ENSG_gene_names_mRNA"))
    expect_error(UpdateObject(object = tTEobject, meta.data = list(ENSG_gene_names_mRNA_data), meta.data.ids = "ENSG_gene_names_mRNA"))
    expect_no_error(UpdateObject(object = tTEobject, meta.data = list(ENSG_gene_names_mRNA_data), meta.data.ids = "ENSG_gene_names_mRNA", overwrite = TRUE))
    expect_no_error(UpdateObject(object = tTEobject, meta.data = ENSG_gene_names_mRNA_data, meta.data.ids = "ENSG_gene_names_mRNA", overwrite = TRUE))
    expect_no_error(UpdateObject(object = tTEobject, meta.data = list(ENSG_gene_names_mRNA_data), meta.data.ids = list("ENSG_gene_names_mRNA"), overwrite = TRUE))
})
