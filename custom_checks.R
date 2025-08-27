# Load required devtools tools
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocCheck", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("BiocCheck")
}
if (!requireNamespace("spelling", quietly = TRUE)) install.packages("spelling")

pkg <- "."  # package root (adjust if needed)

cat("\nRunning devtools::check() ...\n")
devtools::check(pkg, cran = TRUE, document = TRUE)

cat("\nChecking Bioconductor compliance ...\n")
BiocCheck::BiocCheck(pkg)

cat("\nChecking DESCRIPTION fields ...\n")
desc <- read.dcf("DESCRIPTION")
required_fields <- c("Package","Title","Description","Version","Authors@R","License")
missing <- setdiff(required_fields, colnames(desc))
if (length(missing)) {
  warning("Missing fields in DESCRIPTION: ", paste(missing, collapse = ", "))
} else {
  message("DESCRIPTION has all required fields")
}

cat("\nChecking spelling ...\n")
spelling::spell_check_package(pkg, ignore = readLines("inst/WORDLIST", warn = FALSE))

cat("\nChecking installed size ...\n")
devtools::install(pkg, quiet = TRUE, upgrade = "never")
libpath <- .libPaths()[1]
pkgpath <- file.path(libpath, desc[1, "Package"])
size <- sum(file.info(list.files(pkgpath, recursive = TRUE, full.names = TRUE))$size) / 1e6
message("Installed size: ", round(size, 2), " MB")
if (size > 5) warning("Installed size > 5 MB (CRAN may complain)")

cat("\nChecking number of imports ...\n")
imports <- strsplit(desc[1, "Imports"], ",")[[1]]
imports <- trimws(imports)
message("Number of Imports: ", length(imports))
if (length(imports) > 20) warning("Too many Imports (CRAN may suggest moving to Suggests)")

# --- Extra: check .rda object sizes ---
if (dir.exists(file.path(pkg, "data"))) {
  cat("\nChecking dataset sizes (data/*.rda)...\n")
  rda_files <- list.files(file.path(pkg, "data"), pattern = "\\.rda$", full.names = TRUE)
  for (f in rda_files) {
    e <- new.env()
    load(f, envir = e)
    objs <- ls(e)
    for (obj in objs) {
      sz <- object.size(e[[obj]]) / 1e6
      message(basename(f), " -> ", obj, " (", round(sz, 2), " MB)")
    }
  }
}

cat("\nPrecheck complete!\n")
