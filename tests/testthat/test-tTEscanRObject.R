test_that("Update_tTEscanR_Object adds new counts and metadata correctly", {
  initial_counts <- matrix(1:4, nrow = 2, dimnames = list(c("gene1", "gene2"), c("s1", "s2")))
  obj <- Create_tTEscanR_Object(counts = initial_counts, assay = "mRNA")

  # CASE 1: Add new counts
  new_counts <- matrix(5:8, nrow = 2, dimnames = list(c("gene1", "gene2"), c("s1", "s2")))

  expect_error(updated_obj <- Update_tTEscanR_Object(object = obj, counts = new_counts, assay = "protein"), "Invalid assay names detected:")
  expect_error(updated_obj <- Update_tTEscanR_Object(object = obj, counts = new_counts, assay = "mRNA"), "change the input parameter `overwrite` to TRUE")
  expect_no_error(updated_obj <- Update_tTEscanR_Object(object = obj, counts = new_counts, assay = "mRNA", overwrite.assay = TRUE))
  updated_obj <- Update_tTEscanR_Object(object = obj, counts = new_counts, assay = "mRNA", overwrite.assay = TRUE)

  # CASE 2: add new metadata
  new_metadata <- c(20, 30)
  updated_obj <- Update_tTEscanR_Object(object = updated_obj, meta.data = new_metadata, meta.data.ids = "patient_age")

  # Check if assays and metadata were added correctly
  expect_true("mRNA" %in% names(updated_obj@assays))
  expect_equal(updated_obj@assays$mRNA, new_counts)
  expect_true("patient_age" %in% names(updated_obj@meta.data))
  expect_equal(updated_obj@meta.data$patient_age, new_metadata)
})

test_that("Update_tTEscanR_Object handles overwrite.assay correctly", {

  obj <- Create_tTEscanR_Object(counts = matrix(1:4, nrow = 2), assay = "mRNA")

  # CASE 1: Attempt to overwrite without setting overwrite = TRUE (should fail)
  expect_error( Update_tTEscanR_Object(object = obj, counts = matrix(5:8, nrow = 2), assay = "mRNA"), "Specify another name for mRNA or change the input parameter `overwrite` to TRUE.")

  # CASE 2: Overwrite with overwrite = TRUE (should succeed)
  new_RNA_counts <- matrix(9:12, nrow = 2)
  updated_obj <- Update_tTEscanR_Object(object = obj, counts = new_RNA_counts, assay = "mRNA", overwrite.assay = TRUE)

  expect_equal(updated_obj@assays$mRNA, new_RNA_counts)
})

test_that("Update_tTEscanR_Object handles named and unnamed metadata lists", {
  expect_error(obj <- Create_tTEscanR_Object(), "argument \"counts\" is missing, with no default")

  initial_counts <- matrix(1:4, nrow = 2, dimnames = list(c("gene1", "gene2"), c("s1", "s2")))
  obj <- Create_tTEscanR_Object(counts = initial_counts, assay = "mRNA")

  # CASE 1: Add metadata with a named list
  named_meta <- list(age = c(25, 30), height = c(170, 180))
  updated_obj_named <- Update_tTEscanR_Object(object = obj, meta.data = named_meta)

  expect_equal(updated_obj_named@meta.data, named_meta)

  # CASE 2: Add metadata with an unnamed list and meta.data.ids
  unnamed_meta <- list(c("m", "f"), c("urban", "rural"))
  updated_obj_unnamed <- Update_tTEscanR_Object(object = obj, meta.data = unnamed_meta, meta.data.ids = c("gender", "location"))

  expected_unnamed_meta <- list(gender = c("m", "f"), location = c("urban", "rural"))
  expect_equal(updated_obj_unnamed@meta.data, expected_unnamed_meta)
})

test_that("Update_tTEscanR_Object handles overwrite.assay correctly", {
  obj <- new("tTEscanR_Object", assays = list(mRNA = matrix(1:4, nrow = 2)))

  # Test Case 2: Attempt to overwrite without setting overwrite = TRUE (should fail)
  expect_error(
    Update_tTEscanR_Object(object = obj, counts = matrix(5:8, nrow = 2), assay = "mRNA"),
    "Specify another name for mRNA or change the input parameter `overwrite` to TRUE."
  )

  # Test Case 3: Overwrite with overwrite = TRUE (should succeed)
  new_RNA_counts <- matrix(9:12, nrow = 2)
  updated_obj <- Update_tTEscanR_Object(object = obj, counts = new_RNA_counts, assay = "mRNA", overwrite.assay = TRUE)

  expect_equal(updated_obj@assays$mRNA, new_RNA_counts)
})

test_that("The tTEscanR object is correctly generated", {

  data(subset_mRNA_data)

  tTEobject <- Create_tTEscanR_Object(counts = subset_mRNA_data, assay = "mRNA")
  cell_types <- colnames(subset_mRNA_data)

  expect_s4_class(object = tTEobject, class = "tTEscanR_Object")
  expect_equal(as.data.frame(tTEobject@assays$mRNA), subset_mRNA_data)
  expect_no_error(IsIn_tTEscanR_Object(tTEobject, "assays", "mRNA"))

  # CASE 1: no error - counts can be a list or a single input
  expect_no_error(Create_tTEscanR_Object(counts = subset_mRNA_data, assay = "mRNA"))
  expect_no_error(Create_tTEscanR_Object(counts = list(mRNA = subset_mRNA_data)))
  expect_no_error(Create_tTEscanR_Object(counts = list(subset_mRNA_data), assay = list("mRNA")))

  # CASE 2: error - counts has a wrong data format
  expect_error(Create_tTEscanR_Object(counts = subset_mRNA_data))
  expect_error(Create_tTEscanR_Object(counts = c(1, 2, 3, 4, 5), assay = "mRNA"))
  expect_error(Create_tTEscanR_Object(counts = c(subset_mRNA_data), assay = "mRNA"))

  # CASE 3: error - wrong usage of the assay parameter
  expect_error(Create_tTEscanR_Object(counts = subset_mRNA_data, assay = mRNA))
  expect_error(Create_tTEscanR_Object(counts = subset_mRNA_data, assay = c("mRNA", "mRNA_2")))

  # CASE 4: no error - incorporation of metadata
  expect_no_error(Create_tTEscanR_Object(counts = subset_mRNA_data, assay = "mRNA", meta.data = cell_types, meta.data.ids = "cell.types"))
  expect_no_error(Create_tTEscanR_Object(counts = list(mRNA = subset_mRNA_data), meta.data = cell_types, meta.data.ids = "cell.types"))
  expect_no_error(Create_tTEscanR_Object(counts = list(subset_mRNA_data), assay = list("mRNA"), meta.data = cell_types, meta.data.ids = "cell.types"))
  expect_error(Create_tTEscanR_Object(counts = subset_mRNA_data, meta.data = list(cell_types), meta.data.ids = list("cell.types")))
  expect_no_error(Create_tTEscanR_Object(counts = list(subset_mRNA_data), assay = "mRNA", meta.data = cell_types, meta.data.ids = "cell.types"))

  # CASE 5: error - repeated ids
  expect_error(Create_tTEscanR_Object(counts = list(mRNA = subset_mRNA_data), assay = list("mRNA"), meta.data = cell_types, meta.data.ids = "cell.types"))
  expect_error(Create_tTEscanR_Object(counts = list(subset_mRNA_data), assay = list("mRNA"), meta.data = list(cell_types = cell_types), meta.data.ids = list("cell.types")))
  expect_error(Create_tTEscanR_Object(counts = list(subset_mRNA_data, subset_mRNA_data), assay = list("mRNA", "mRNA"), meta.data = cell_types, meta.data.ids = "cell.types"))
})

test_that("The tTEscanR object is correctly updated", {

  tTEobject <- Create_tTEscanR_Object(counts = subset_mRNA_data, assay = "mRNA")
  cell_types <- colnames(subset_mRNA_data)

  expect_s4_class(object = tTEobject, class = "tTEscanR_Object")
  expect_equal(as.data.frame(tTEobject@assays$mRNA), subset_mRNA_data)

  # CASE 0: check if individual sections are in the tTEscanR object
  expect_no_error(IsIn_tTEscanR_Object(tTEobject, slot = "assays", section = "mRNA"))
  expect_error(IsIn_tTEscanR_Object(tTEobject, slot = "assays", section = "tRNA"))

  # CASE 1: proper way to update the tTEscanR object
  expect_no_error(Update_tTEscanR_Object(object = tTEobject, counts = subset_tRNA_data, assay = "tRNA")) # just counts
  expect_no_error(Update_tTEscanR_Object(object = tTEobject, meta.data = cell_types, meta.data.ids = "cell_types")) # just metadata
  expect_no_error(Update_tTEscanR_Object(object = tTEobject, counts = subset_tRNA_data, assay = "tRNA", meta.data = cell_types, meta.data.ids = "cell_types")) # counts and metadata

  # CASE 2: use of overwrite
  tTEobject <- Update_tTEscanR_Object(object = tTEobject, counts = subset_tRNA_data, assay = "tRNA")
  expect_s4_class(object = tTEobject, class = "tTEscanR_Object")
  expect_equal(tTEobject@assays$tRNA, subset_tRNA_data)
  expect_no_error(IsIn_tTEscanR_Object(tTEobject, slot = "assays", section = "mRNA"))
  expect_no_error(IsIn_tTEscanR_Object(tTEobject, slot = "assays", section = "tRNA"))
  expect_error(Update_tTEscanR_Object(object = tTEobject, counts = subset_tRNA_data, assay = "tRNA")) # tRNA already in tTEobject
  expect_no_error(Update_tTEscanR_Object(object = tTEobject, counts = subset_tRNA_data, assay = "tRNA", overwrite.assay = TRUE))
  expect_error(Update_tTEscanR_Object(object = tTEobject, counts = list(subset_mRNA_data, subset_tRNA_data), assay = list("mRNA", "tRNA")))
  expect_no_error(Update_tTEscanR_Object(object = tTEobject, counts = list(subset_mRNA_data, subset_tRNA_data), assay = list("mRNA", "tRNA"),  overwrite.assay = TRUE))

  # CASE 3: different input formats (list or single)
  expect_no_error(Update_tTEscanR_Object(object = tTEobject, counts = list(subset_mRNA_data), assay = list("mRNA"),  overwrite.assay = TRUE))
  expect_no_error(Update_tTEscanR_Object(object = tTEobject, counts = list(subset_mRNA_data), assay = "mRNA",  overwrite.assay = TRUE)) # inconsistencies list-single
  expect_error(Update_tTEscanR_Object(object = tTEobject, counts = subset_mRNA_data, assay = list("mRNA"),  overwrite.assay = TRUE)) # inconsistencies single-list
  expect_error(Update_tTEscanR_Object(object = tTEobject, counts = c(subset_mRNA_data, subset_tRNA_data), assay = c("mRNA", "tRNA"),  overwrite.assay = TRUE)) # using vectors

  cell_types_tRNA <- colnames(subset_tRNA_data)
  expect_no_error(Update_tTEscanR_Object(object = tTEobject, meta.data = list(cell_types), meta.data.ids = list("cell_types")))
  expect_no_error(Update_tTEscanR_Object(object = tTEobject, meta.data = list(cell_types, cell_types_tRNA), meta.data.ids = list("cell_types", "cell_types_tRNA")))
  expect_error(Update_tTEscanR_Object(object = tTEobject, meta.data = list(cell_types, cell_types_tRNA), meta.data.ids = list("cell_types", "cell_types"))) # same ids
  expect_error(Update_tTEscanR_Object(object = tTEobject, meta.data = c(cell_types, cell_types_tRNA), meta.data.ids = c("cell_types", "cell_types_tRNA"))) # using vectors

  # CASE 4: error - missing parameters
  expect_error(Update_tTEscanR_Object(object = tTEobject)) # nothing to add
  expect_error(Update_tTEscanR_Object(object = tTEobject, assay = "tRNA")) # nothing to add
  expect_error(Update_tTEscanR_Object(object = tTEobject, meta.data.ids = "cell_types")) # nothing to add
  expect_error(Update_tTEscanR_Object(object = tTEobject, counts = subset_tRNA_data)) # no ids
  expect_error(Update_tTEscanR_Object(object = tTEobject, meta.data = subset_tRNA_data)) # no ids
  expect_error(Update_tTEscanR_Object(object = tTEobject, meta.data = cell_types)) # no ids
})

test_that("The function to find sections in a tTEscanR object works", {

  tTEobject <- Create_tTEscanR_Object(counts = subset_mRNA_data, assay = "mRNA")
  tTEobject <- Update_tTEscanR_Object(object = tTEobject, counts = subset_tRNA_data, assay = "tRNA")
  tTEobject <- Update_tTEscanR_Object(object = tTEobject, meta.data = colnames(subset_mRNA_data), meta.data.ids = "cell.types")

  expect_s4_class(object = tTEobject, class = "tTEscanR_Object")
  expect_equal(as.data.frame(tTEobject@assays$mRNA), subset_mRNA_data)
  expect_equal(tTEobject@assays$tRNA, subset_tRNA_data)

  # CASE 1: no error - right slot and section
  expect_no_error(IsIn_tTEscanR_Object(object = tTEobject, slot = "assays", section = "mRNA"))
  expect_no_error(IsIn_tTEscanR_Object(object = tTEobject, slot = "assays", section = "tRNA"))
  expect_no_error(IsIn_tTEscanR_Object(object = tTEobject, slot = "meta.data", section = "cell.types"))

  expect_no_error(IsIn_tTEscanR_Object(object = tTEobject, slot = "assays", section = "mRNA", update.assay = FALSE, overwrite = FALSE))
  expect_no_error(IsIn_tTEscanR_Object(object = tTEobject, slot = "assays", section = "tRNA", update.assay = FALSE, overwrite = FALSE))
  expect_no_error(IsIn_tTEscanR_Object(object = tTEobject, slot = "meta.data", section = "cell.types", update.assay = FALSE, overwrite = FALSE))

  # CASE 2: error - wrong slot
  expect_error(IsIn_tTEscanR_Object(object = tTEobject, slot = "meta.data", section = "mRNA"))
  expect_error(IsIn_tTEscanR_Object(object = tTEobject, slot = "meta.data", section = "tRNA"))
  expect_error(IsIn_tTEscanR_Object(object = tTEobject, slot = "assays", section = "cell.types"))

  # CASE 3: error - non available section
  expect_error(IsIn_tTEscanR_Object(object = tTEobject, slot = "assays", section = mRNA))
  expect_error(IsIn_tTEscanR_Object(object = tTEobject, slot = "assays", section = "CodonUsage"))

  # CASE 4: error - non available slot
  expect_error(IsIn_tTEscanR_Object(object = tTEobject, slot = "assay", section = "mRNA"))
  expect_error(IsIn_tTEscanR_Object(object = tTEobject, slot = "counts", section = "tRNA"))
  expect_error(IsIn_tTEscanR_Object(object = tTEobject, slot = "metadata", section = "cell.types"))
})

test_that("The metadata is correctly updated", {

  tTEobject <- Create_tTEscanR_Object(counts = subset_mRNA_data, assay = "mRNA")

  # CASE 1: no error - single addition with and without lists
  expect_no_error(Update_tTEscanR_Object(object = tTEobject, meta.data = ENSG_gene_names_mRNA_data, meta.data.ids = "ENSG_gene_names_mRNA"))
  expect_no_error(Update_tTEscanR_Object(object = tTEobject, meta.data = list(ENSG_gene_names_mRNA_data), meta.data.ids = list("ENSG_gene_names_mRNA")))

  # CASE 2: no error - multiple addition
  expect_no_error(Update_tTEscanR_Object(object = tTEobject, meta.data = list(ENSG_gene_names_mRNA_data, short_gene_names_mRNA_data), meta.data.ids = list("ENSG_gene_names_mRNA", "short_gene_names_mRNA")))
  expect_no_error(Update_tTEscanR_Object(object = tTEobject, meta.data = list(ENSG = ENSG_gene_names_mRNA_data, short = short_gene_names_mRNA_data)))

  # CASE 3: error - inconsistent input format
  expect_no_error(Update_tTEscanR_Object(object = tTEobject, meta.data = list(short_gene_names_mRNA_data), meta.data.ids = "short_gene_names_mRNA"))
  expect_error(Update_tTEscanR_Object(object = tTEobject, meta.data = list(ENSG = ENSG_gene_names_mRNA_data), meta.data.ids = "ENSG_gene_names_mRNA_2"))

  # CASE 4: using the overwrite parameter
  tTEobject <- Update_tTEscanR_Object(object = tTEobject, meta.data = list(ENSG_gene_names_mRNA_data), meta.data.ids = list("ENSG_gene_names_mRNA"))
  expect_error(Update_tTEscanR_Object(object = tTEobject, meta.data = list(ENSG_gene_names_mRNA_data), meta.data.ids = "ENSG_gene_names_mRNA"))
  expect_error(Update_tTEscanR_Object(object = tTEobject, meta.data = list(ENSG_gene_names_mRNA_data), meta.data.ids = "ENSG_gene_names_mRNA", overwrite.assay = TRUE))
  expect_no_error(Update_tTEscanR_Object(object = tTEobject, meta.data = ENSG_gene_names_mRNA_data, meta.data.ids = "ENSG_gene_names_mRNA", overwrite.metadata = TRUE))
  expect_no_error(Update_tTEscanR_Object(object = tTEobject, meta.data = list(ENSG_gene_names_mRNA_data), meta.data.ids = list("ENSG_gene_names_mRNA"), overwrite.metadata = TRUE))
})
