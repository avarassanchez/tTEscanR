PerformTranslation <- function(data_to_translate, notation.from, notation.to, verbose){

  ###
  # CALL: FeaturesToAA()
  # DESCRIPTION: This function takes a vector of features and translates them based on the notation.X parameters.
  # The accepted features are codons and anticodons that can be translated into anticodons, codons or amino acids.
  ###

  data_to_translate <- as.data.frame(data_to_translate)

  # Check that the features are not already amino acids
  if (any(data_to_translate %in% feature_to_aa$aa)) stop("The selected features are already amino acids.")

  # Perform the translation from codons or anticodons to amino acids
  if (notation.from == "codon" || notation.from == "anticodon") {
    if (any(data_to_translate %in% feature_to_aa[[notation.from]])) stop(paste("No overlapping features have been identified belonging to", notation.from, "."))
    if (verbose) message(paste("-", sum(data_to_translate %in% feature_to_aa[[notation.from]]), "out of", length(data_to_translate), "will be translated."))
    colnames(data_to_translate) <- notation.from
    translated_features <- dplyr::left_join(data_to_translate, feature_to_aa, by = notation.from)
  }

  # Keep the columns that we are requested by the notation parameter
  translated_features <- data.frame(translated_features[[notation.from]], translated_features[[notation.to]])
  colnames(translated_features) <- c(notation.from, notation.to)

  return (translated_features) # Return a table with the input and the translated features
}

RetrieveTranslation <- function(format_input, position, translated_features, data_to_translate, notation.to, notation.from){

  ###
  # CALL: FeaturesToAA()
  # DESCRIPTION: This function takes the original and the translated features and properly adjusts the original data and returns the translated data.
  ###

  if (format_input == "table"){ # data_to_translate is a TABLE
    # Update the row names with the translated features
    if (position == "row"){
      indexes <- which(translated_features[[notation.from]] %in% rownames(data_to_translate))
      translated_features <- translated_features[indexes, ]
      rownames(data_to_translate) <- translated_features[[notation.to]]
    # Update the column names with the translated features
    } else if (position == "column"){
      indexes <- which(translated_features[[notation.from]] %in% colnames(data_to_translate))
      translated_features <- translated_features[indexes, ]
      colnames(data_to_translate) <- translated_features[[notation.to]]
    # Update the content of a column
    } else {
      translated_features <- translated_features[translated_features[[notation.from]] %in% data_to_translate[[position]], ]
      index <- match(data_to_translate[[position]], translated_features[[notation.from]])
      data_to_translate[[position]] <- translated_features[[notation.to]][index]
    }
  } else { # data_to_translate is a VECTOR
    indexes <- which(translated_features[[notation.from]] %in% data_to_translate)
    translated_features <- translated_features[indexes, ]
    data_to_translate <- translated_features[[notation.to]]
  }
  return(data_to_translate) # Returns translated data
}

#' Translates Genes Between Annotations or Codons, Anticodons and Amino Acids.
#' @description
#' This function converts: codons and anticodons into anticodons, codons or amino acids based on the genetic code.
#'
#' @param data_to_translate A \code{character} vector or \code{data.frame} containing gene names or features to be translated.
#' @param position Optional; either \code{"row"}, \code{"column"} or \code{<column_name>} to specify the location of the genes or features in \code{data_to_translate}. Required if \code{data_to_transalte} is a \code{data.frame}.
#' @param notation.to Either \code{"codon"}, \code{"anticodon"} or \code{"aa"} to select the output format of the features in \code{data_to_translate}.
#' @param notation.from Either \code{"codon"} or \code{"anticodon"} to select the input format of the features in \code{data_to_translate}.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return Translated features (codons, anticodons or amino acids) from \code{data_to_translate}.
#' @export
#'
#' @examples
#' data(subset_mRNA_data, metadata)
#' tTEscanR_obj <- Create_tTEscanR_Object(counts = subset_mRNA_data, assay = "mRNA",
#'                                        meta.data = list(metadata, "tissue"),
#'                                        meta.data.ids = list("ConditionsLabels", "CorrectionFactor"))
#' tTEscanR_obj <- ComputeCodonUsage(object = tTEscanR_obj, codon_freq = NULL,
#'                                   species = "hg38", filter = "canonical",
#'                                   additional.metrics = FALSE, reduce = 1000)
#' codons <- rownames(tTEscanR_obj@assays$CodonUsage)
#' codons_to_AA <- FeaturesToAA(data_to_translate = codons, notation.from = "codon",
#'                              notation.to = "aa")
#' codons_to_anticodons <- FeaturesToAA(data_to_translate = codons, notation.from = "codon",
#'                                      notation.to = "anticodon")

FeaturesToAA <- function(data_to_translate, position = NULL, notation.from, notation.to, verbose = TRUE){

  ###
  # CALL: User and ConsistencyWithCodonFreq()
  # DESCRIPTION: This function translates codons or anticodons into anticodons, codons or amino acids
  ###

  message("1 . Checking the format of the input data.") # The different input parameters enable to deal with diverse input sources of data
  if (notation.from == notation.to) stop("The same `notation.from` and `notation.to` parameters have been input.")
  if (!(notation.from %in% c("codon", "anticodon"))) stop("Incorrect `notation.from` input value.\n", "Supported formats: codon or anticodon.")
  if (!(notation.to %in% c("codon", "anticodon", "aa"))) stop("Incorrect `notation.to` input value.\n", "Supported formats: codon, anticodon, aa.")
  if (is.null(data_to_translate)) stop("Incorrect `data_to_translate` input parameter.")
  format_input <- IdentifyInputFormat(data = data_to_translate) # Check if the data_to_translate is a table or a vector

  if (format_input == "table"){ # data_to_translate is a TABLE
    CheckDataFrame(data = data_to_translate, names = TRUE) # Ensure that the data has proper column or row names
    # Check the position parameter
    if (is.null(position) && (!(position %in% c("row", "column"))) && (!(position %in% colnames(data_to_translate)))) stop("Specify a value for the `position` input.\n",
                                                                                                                           "Supported formats: row, column or column name in data_to_translate.")
    message("   1 . COMPLETED\n", "2 . Extracting the vector of features to translate.")
    features <- if (position == "row") rownames(data_to_translate) else if (position == "column") colnames(data_to_translate) else data_to_translate[[position]]

  } else{ # data_to_translate is a VECTOR
    message("   1 . COMPLETED\n", "2 . Extracting the vector of features to translate.")
    features <- data_to_translate
  }

  # Check that the features have been properly retrieved
  if (length(features) == 0 || all(is.na(features))) stop("The selected data to translate is empty or contains only missing values.")

  message("   2 . COMPLETED\n", "3 . Translating the features.")
  translated_features <- PerformTranslation(data_to_translate = features, notation.from = notation.from, notation.to = notation.to, verbose = verbose)

  message("   3 . COMPLETED\n", "4 . Retrieving the translated input data.")
  data_to_translate <- RetrieveTranslation(format_input = format_input, position = position, data_to_translate = data_to_translate,
                                           translated_features = translated_features, notation.to = notation.to, notation.from = notation.from)
  message("   4 . COMPLETED")
  return (data_to_translate) # Return the translated input data
}

IsEnsemblID <- function(gene_vector, action, notation = NULL) {

  ###
  # CALL: CheckGeneAnnotation()
  # DESCRIPTION: This function evaluates if the format of a gene vector matches the Ensembl-like annotation format.
  ###

  # Checking input parameters
  format <- IdentifyInputFormat(data = gene_vector, mode = "flexible") # Format of the input gene_vector
  if (format != "vector") stop("Please, provide a proper `gene_vector` parameter.\n", "Supported format: vector")
  if ((length(gene_vector) != 1) && (any(!sapply(gene_vector, is.character)))) stop("Element in `gene_vector` is not a string.") # Check that all elements are strings

  # Returns TRUE if the gene follows the Ensembl annotation and FALSE otherwise - gives one output per gene in the gene_vector
  annotation <- grepl("^ENS[A-Z]*[0-9]+$", gene_vector)

  # Evaluates if the Ensembl format contains the gene version
  if (any(grepl("^ENS[A-Z]*G\\d+\\.\\d+$", gene_vector))) stop("There are genes containing the gene version. Those cases can not be handled.")

  if (length(unique(annotation)) > 1){ # Implies that there were TRUE and FALSE outcomes
    if (action == "check") stop("There are gene annotation inconsistencies in the gene vector.")

    # If the translate action is enabled it will report an error in those cases that makes incompatible to modify the gene annotation format
    # Otherwise it will return the format in which the gene_vector would be translated to
    if (action == "translate"){
      if (is.null(notation)) stop("There are annotation inconsistencies in the gene vector.\n", "Please, specify a suitable `notation` parameter.\n", "Supported formats: id for Ensembl ID and symbol for gene name.")
      if (!(notation %in% c("id", "symbol"))) stop("Please, specify a suitable `notation` parameter.\n", "Supported formats: id for Ensembl ID and symbol for gene name.")
      if (notation == "id") return("Gene Symbol") else return("Ensembl ID")
    }

  } else { # All the genes are in the same format
    if (all(annotation)) return("Ensembl ID") else return("Gene Symbol") # Returns the format in which the genes are
  }
}
