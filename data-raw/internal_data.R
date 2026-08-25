library(tidyverse)
library(devtools)
library(tRNAscanImport)

###
# GENERATE THE CODONS x GENETIC CODE TABLE
###
bases <- c("T", "C", "A", "G")
codon_table <- expand.grid(Third = bases, Second = bases, First = bases) %>%
    mutate(Codon = paste0(First, Second, Third)) %>%
    pull(Codon)

raw_data <- read.delim("/Users/anavaras/Downloads/tTEscanR_default_data/gencode.dmp",
    sep = "|",
    header = TRUE,
    strip.white = TRUE,
    stringsAsFactors = FALSE
)

gen_codes <- raw_data[, c(1, 3, 4)]
colnames(gen_codes) <- c("ID", "Name", "AA_String")

final_matrix_genetic_code <- gen_codes %>%
    rowwise() %>%
    mutate(AA_List = list(strsplit(AA_String, "", fixed = TRUE)[[1]])) %>%
    ungroup() %>%
    unnest(AA_List) %>%
    mutate(Codon = rep(codon_table, times = nrow(gen_codes))) %>%
    select(Codon, Name, AA_List) %>%
    pivot_wider(names_from = Name, values_from = AA_List)

final_matrix_genetic_code$Anticodon <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(final_matrix_genetic_code$Codon)))

###
# EXTRACT THE DEFAULT CODON FREQUENCY MATRICES
###
codon_freq_results_canonical_hg38 <- getCodonFreq(
    dataset_name = "hsapiens_gene_ensembl",
    filter = "canonical",
    retain_geneversion = TRUE,
    retain_mitochondrial = FALSE,
    out_format = "external_gene_name"
)
codon_freq_table_canonical_hg38 <- codon_freq_results_canonical_hg38$codon_freq_per_gene_matrix

codon_freq_results_canonical_mm39 <- getCodonFreq(
    dataset_name = "mmusculus_gene_ensembl",
    filter = "canonical",
    retain_geneversion = TRUE,
    retain_mitochondrial = FALSE,
    out_format = "external_gene_name"
)
codon_freq_table_canonical_mm39 <- codon_freq_results_canonical_mm39$codon_freq_per_gene_matrix

###
# GENERAL VARIABLES
###
amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
three_letter_amino_acids <- c("Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr")

assay_map_AA <- list(
    demand = list(assay_id = list("AADemand"), check_assays = list("CodonUsage"), AAfunction = "sense"),
    supply = list(assay_id = list("AASupply"), check_assays = list("AnticodonUsage"), AAfunction = "reverse"),
    both = list(assay_id = list("AADemand", "AASupply"), check_assays = list("CodonUsage", "AnticodonUsage"), AAfunction = list("sense", "reverse"))
)

assay_map_TE <- list(
    codon = list(mRNA = "CodonUsage", tRNA = "AnticodonUsage", id = "tTEresults_codon"),
    aa = list(mRNA = "AADemand", tRNA = "AASupply", id = "tTEresults_AA")
)

###
# tRNA SPECIFIC PROCESSING DATA
###
# Downloaded from GtRNAdb: http://gtrnadb.ucsc.edu/index.html
hg38_set_path <- "/Users/anavaras/Downloads/tTEscanR_default_data/hg38-tRNAs-confidence-set.ss"
hg38_set <- import.tRNAscanAsGRanges(hg38_set_path)

hg38_map_path <- "/Users/anavaras/Downloads/tTEscanR_default_data/hg38-tRNAs_name_map.txt"
hg38_map <- read.delim(file = hg38_map_path)

mm39_set_path <- "/Users/anavaras/Downloads/tTEscanR_default_data/mm39-tRNAs-confidence-set.ss"
mm39_set <- import.tRNAscanAsGRanges(mm39_set_path)

mm39_map_path <- "/Users/anavaras/Downloads/tTEscanR_default_data/mm39-tRNAs_name_map.txt"
mm39_map <- read.delim(file = mm39_map_path)

###
# SAVE THE VARIABLES
###
usethis::use_data(final_matrix_genetic_code,
    codon_freq_table_canonical_hg38, codon_freq_table_canonical_mm39,
    amino_acids, three_letter_amino_acids,
    assay_map_AA, assay_map_TE,
    hg38_set, hg38_map, mm39_set, mm39_map,
    overwrite = TRUE, internal = TRUE
)
