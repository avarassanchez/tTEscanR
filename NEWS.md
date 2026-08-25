# tTEscanR 0.99.0

## Initial features 
* Initial submission to Bioconductor.
* Added core functionality for the computation of translation efficiency
from sequencing data.
* Included vignettes and documentation.

# tTEscanR 0.99.1

## New features
* Redesign the tTEscanR object as a `MultiAssayExperiment()`.

## Bug fixes
* Updated `Authors@R` including funding details.
* Removed the `as.character()` calls around numeric variables in `message()`.
* Optimized test cases to avoid redundancy.
* Updated `style.css` file for the vignettes.

## Refactoring
* Deprecated custom `savePlot()` wrapper in favor of direct `ggplot2::ggsave()` call.
* Removed own color palette stored internally to be used in `getSafeColorScale()`.
* Substitute `ggradar()` by standard `ggplot2` functions.
