#' mRNA Expression Data Subset
#'
#' This dataset contains an example of mRNA gene expression data with
#' protein-coding genes as rows and cell types as columns.
#'
#' @name default_tTEscanR_mRNA_data
#' @docType data
#' @format A matrix where rows represent genes and columns represent individual
#' samples or cell types.
#' @usage data(default_tTEscanR_mRNA_data)
"default_tTEscanR_mRNA_data"

#' tRNA Expression Data Subset
#'
#' This dataset contains an example of tRNA gene expression data with tRNA
#' genes as rows and cell types as columns.
#'
#' @name default_tTEscanR_tRNA_data
#' @docType data
#' @format A matrix or data frame where rows represent tRNA gene names
#' and columns represent individual samples or cell types.
#' @usage data(default_tTEscanR_tRNA_data)
"default_tTEscanR_tRNA_data"

#' Metadata of mRNA and tRNA data
#'
#' This dataset contains the extra information to indicate the conditions of
#' the data
#'
#' @format A data frame with samples as rows and annotation columns such as:
#' \describe{
#'   \item{tissue}{The sequencing batch or sample origin}
#'   \item{cell.type}{The experimental group or cell type}
#'   \item{conditions}{The combination of the previous items as will be referred
#'   in the columns of the count matrices}
#' }
#' @name default_tTEscanR_metadata
#' @docType data
#' @usage data(default_tTEscanR_metadata)
"default_tTEscanR_metadata"
