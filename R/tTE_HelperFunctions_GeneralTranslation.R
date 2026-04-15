PerformTranslation <- function(data_to_translate, notation_from, notation_to, genetic_code, verbose){
  ###
  # CALL: FeaturesToAA()
  # DESCRIPTION: This function takes a vector of features and translates them based on the notation_X parameters.
  # The accepted features are codons and anticodons that can be translated into anticodons, codons or amino acids.
  ###

  data_to_translate <- gsub("U", "T", data_to_translate)

  col_from <- switch(notation_from, "codon" = "Codon", "anticodon" = "Anticodon")
  col_to <- switch(notation_to, "codon" = "Codon", "anticodon" = "Anticodon", "aa" = genetic_code)
  if (notation_to == "aa" && !(col_to %in% colnames(final_matrix_genetic_code))) {
    stop("The genetic code '", genetic_code, "' is not a valid column in final_matrix_genetic_code.")
  }

  dict_from <- final_matrix_genetic_code[[col_from]]
  dict_to <- final_matrix_genetic_code[[col_to]]

  index <- match(data_to_translate, dict_from)
  valid_index <- !is.na(index)
  matches <- sum(valid_index)

  if (matches == 0) stop("No overlapping features identified for ", notation_from)
  if (verbose) message(sprintf("- %d out of %d will be translated.", matches, length(data_to_translate)))

  translated_features <- data.frame(from = data_to_translate, to = dict_to[index], stringsAsFactors = FALSE)
  colnames(translated_features) <- c(notation_from, notation_to)

  return (translated_features) # Return a table with the input and the translated features
}

RetrieveTranslation <- function(format_input, position, translated_features, data_to_translate, notation_to, notation_from){

  ###
  # CALL: FeaturesToAA()
  # DESCRIPTION: This function takes the original and the translated features and properly adjusts the original data and returns the translated data.
  ###

  new_values <- translated_features[[notation_to]]
  has_translation <- !is.na(new_values) # Update only the names that have a valid translation

  if (format_input){
    if (position == "row"){

      original_names <- rownames(data_to_translate) # Get the original row names
      original_names[has_translation] <- new_values[has_translation]
      # if (any(duplicated(original_names))) warning("Translating to amino acids resulted in duplicate row names.
      #                                              Depending on the object type, R may append unique numbers as indexes or throw an error.")
      rownames(data_to_translate) <- original_names

    } else if (position == "column"){
      original_names <- colnames(data_to_translate)
      original_names[has_translation] <- new_values[has_translation]
      colnames(data_to_translate) <- original_names

    } else { # Update specific column in-place
      original_names <- data_to_translate[[position]]
      data_to_translate[[position]] <- ifelse(is.na(new_values), original_names, new_values)
    }

  } else { # Working with a vector
    data_to_translate <- ifelse(is.na(new_values), data_to_translate, new_values)
  }

  return(data_to_translate)
}

#' Translates Genes Between Annotations or Codons, Anticodons and Amino Acids.
#' @description
#' This function converts: codons and anticodons into anticodons, codons or amino acids based on the genetic code.
#'
#' @param data_to_translate A \code{character} vector or \code{data.frame} containing gene names or features to be translated.
#' @param position Optional; either \code{"row"}, \code{"column"} or \code{<column_name>} to specify the location of the genes or features in \code{data_to_translate}. Required if \code{data_to_transalte} is a \code{data.frame}.
#' @param genetic_code A \code{character} string to specify the genetic code to be used. Defaults to \code{"Standard"}.
#' @param notation_to Either \code{"codon"}, \code{"anticodon"} or \code{"aa"} to select the output format of the features in \code{data_to_translate}.
#' @param notation_from Either \code{"codon"} or \code{"anticodon"} to select the input format of the features in \code{data_to_translate}.
#' @param verbose Logical; if \code{TRUE}, displays information messages. Defaults to \code{TRUE}.
#'
#' @return Translated features (codons, anticodons or amino acids) from \code{data_to_translate}.
#' @export
#'
#' @examples
#' data(default_tTEscanR_mRNA_data, default_tTEscanR_metadata)
#' tTEscanR_obj <- Create_tTEscanR_Object(counts = default_tTEscanR_mRNA_data, assay = "mRNA",
#'                                        meta.data = list(default_tTEscanR_metadata, "tissue"),
#'                                        meta.data.ids = list("ConditionsLabels", "CorrectionFactor"))
#' tTEscanR_obj <- ComputeCodonUsage(object = tTEscanR_obj, codon_freq = NULL,
#'                                   species = "hg38", additional_metrics = FALSE)
#' codons <- rownames(tTEscanR_obj@assays$CodonUsage)
#' codons_to_AA <- FeaturesToAA(data_to_translate = codons, notation_from = "codon",
#'                              notation_to = "aa")
#' codons_to_anticodons <- FeaturesToAA(data_to_translate = codons, notation_from = "codon",
#'                                      notation_to = "anticodon")

FeaturesToAA <- function(data_to_translate, position = NULL, genetic_code = "Standard", notation_from = c("codon", "anticodon"), notation_to = c("aa", "anticodon", "codon"), verbose = TRUE){

  ###
  # CALL: User and ConsistencyWithCodonFreq()
  # DESCRIPTION: This function translates codons or anticodons into anticodons, codons or amino acids
  ###

  message("--- Translating the input features ---", "\n1 . Checking the format of the input data.")
  notation_from <- match.arg(notation_from)
  notation_to <- match.arg(notation_to)
  valid_positions <- c("row", "column", colnames(data_to_translate))
  if (notation_from == notation_to) stop("The same 'notation_from' and 'notation_to' parameters have been input.")

  is_table <- is.data.frame(data_to_translate) || is.matrix(data_to_translate)

  message("1 . COMPLETED\n", "2 . Extracting the features to translate.")
  if (is_table) {
    if (is.null(position)) stop("Parameter 'position' required for tabular data.")
    if (!position %in% valid_positions) stop("Incorrect value for 'position' input.\n", "Supported formats: 'row', 'column', or a valid column name.")
    features <- switch(position, "row" = rownames(data_to_translate), "column" = colnames(data_to_translate), data_to_translate[[position]])
  } else {
    features <- data_to_translate
  }

  message("2 . COMPLETED\n", "3 . Translating the features.")
  translated_features <- PerformTranslation(data_to_translate = features, notation_from = notation_from, notation_to = notation_to, genetic_code = genetic_code, verbose = verbose)
  data_to_translate <- RetrieveTranslation(format_input = is_table, position = position, data_to_translate = data_to_translate,
                                           translated_features = translated_features, notation_to = notation_to, notation_from = notation_from)
  message("3 . COMPLETED\n", "--- The features have been successfully translated ---")
  return (data_to_translate)
}

IsEnsemblID <- function(gene_vector) {

  ###
  # CALL: CheckGeneAnnotation()
  # DESCRIPTION: This function evaluates if the format of a gene vector matches the Ensembl-like annotation format.
  # Returns the format in which the gene_vector would be translated to.
  ###

  if (!is.character(gene_vector) || !is.vector(gene_vector)) stop("Input 'gene_vector' must be a character vector.")  # Checking input parameters

  matches <- grepl("^ENS[A-Z]*G[0-9]+(\\.[0-9]+)?$", gene_vector, perl = TRUE)

  if (any(matches)) if (any(grepl("\\.", gene_vector[matches]))) stop ("Gene versions detected. Strip them first.") # Evaluates if the Ensembl format contains the gene version

  all_ens <- all(matches)
  any_ens <- any(matches)

  if (any_ens && !all_ens) stop ("There are gene annotation inconsistencies in the gene vector.")

  return(if (all_ens) "Ensembl ID" else "Gene Symbol") # Same annotation in all genes
}
