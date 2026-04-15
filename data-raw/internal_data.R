library(tidyverse)
library(devtools)

###
# GENERATE THE CODONS x GENETIC CODE TABLE
###
bases <- c("T", "C", "A", "G")
codon_table <- expand.grid(Third = bases, Second = bases, First = bases) %>%
  mutate(Codon = paste0(First, Second, Third)) %>%
  pull(Codon)

raw_data <- read.delim("/Users/anavarassanchez/Downloads/taxdmp/gencode.dmp",
                       sep = "|",
                       header = TRUE,
                       strip.white = TRUE,
                       stringsAsFactors = FALSE)

gen_codes <- raw_data[, c(1, 3, 4)]
colnames(gen_codes) <- c("ID", "Name", "AA_String")

final_matrix_genetic_code <- gen_codes %>% rowwise() %>% mutate(AA_List = list(strsplit(AA_String, "")[[1]])) %>%
  ungroup() %>% unnest(AA_List) %>% mutate(Codon = rep(codon_table, times = nrow(gen_codes))) %>% select(Codon, Name, AA_List) %>%
  pivot_wider(names_from = Name, values_from = AA_List)

final_matrix_genetic_code$Anticodon <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(final_matrix_genetic_code$Codon)))
###
# EXTRACT THE DEFAULT CODON FREQUENCY MATRICES
###
codon_freq_results_canonical_hg38 <- GetCodonFreq(dataset_name = "hsapiens_gene_ensembl", filter = "canonical", retain_geneversion = TRUE,
                                                  retain_mitochondrial = FALSE, out_format = "external_gene_name")
codon_freq_table_canonical_hg38 <- codon_freq_results_canonical_hg38$codon_freq_per_gene_matrix

codon_freq_results_canonical_mm39 <- GetCodonFreq(dataset_name = "mmusculus_gene_ensembl", filter = "canonical", retain_geneversion = TRUE,
                                                  retain_mitochondrial = FALSE, out_format = "external_gene_name")
codon_freq_table_canonical_mm39 <- codon_freq_results_canonical_mm39$codon_freq_per_gene_matrix

###
# DEFINE mRNA - tRNA TEST DATA
###
genes <- c("ADCY9", "YIPF1", "TMEM106A", "GFPT2", "GOT1L1", "CCHCR1", "CAMKK2",
           "ROGDI", "CEP89", "HSDL2", "NOC3L", "COL5A2", "ANO4", "CCRL2")

genes_translated <- c("ENSG00000162104", "ENSG00000058799", "ENSG00000184988", "ENSG00000131459", "ENSG00000169154", "ENSG00000204536", "ENSG00000110931",
                      "ENSG00000067836", "ENSG00000121289", "ENSG00000119471", "ENSG00000173145", "ENSG00000204262", "ENSG00000151572", "ENSG00000121797")

cell_types <- c("Adrenal-Adrenocortical cells", "Adrenal-Chromaffin cells", "Adrenal-CSH1_CSH2 positive cells",
                "Adrenal-Erythroblasts", "Adrenal-Lymphoid cells", "Adrenal-Megakaryocytes",
                "Adrenal-Myeloid cells", "Adrenal-Schwann cells", "Adrenal-SLC26A4_PAEP positive cells",
                "Adrenal-Stromal cells", "Adrenal-Sympathoblasts", "Adrenal-Vascular endothelial cells",
                "Kidney-Erythroblasts", "Kidney-Lymphoid cells", "Kidney-Megakaryocytes",
                "Kidney-Mesangial cells", "Kidney-Metanephric cells", "Kidney-Myeloid cells",
                "Kidney-Stromal cells", "Kidney-Ureteric bud cells")

mRNA_values <- c(3186, 939, 1, 6, 3, 5, 5, 18, 0, 436, 145, 391, 1, 9, 2, 4752, 7770, 20, 125, 583,
                 3587, 161, 1, 17, 6, 14, 35, 10, 3, 64, 14, 97, 0, 11, 2, 664, 1742, 75, 10, 181,
                 833, 58, 0, 13, 3, 8, 126, 5, 1, 118, 17, 84, 1, 10, 6, 788, 1477, 525, 24, 80,
                 456, 31, 0, 0, 0, 1, 1, 4, 1, 55, 3, 11, 0, 0, 0, 1532, 2230, 4, 42, 36,
                 8, 8, 0, 0, 0, 0, 0, 0, 0, 1, 9, 0, 0, 0, 0, 55, 13, 0, 0, 1,
                 1784, 91, 0, 19, 1, 6, 10, 4, 3, 74, 12, 106, 3, 3, 0, 758, 2736, 33, 22, 218,
                 2103, 316, 0, 23, 2, 13, 36, 12, 1, 159, 20, 226, 0, 20, 6, 2003, 6574, 122, 73, 310,
                 8731, 318, 0, 44, 7, 23, 32, 12, 0, 106, 38, 127, 9, 8, 4, 534, 1400, 101, 29, 116,
                 5535, 476, 0, 20, 1, 10, 26, 15, 2, 280, 115, 313, 4, 12, 2, 3595, 9981, 98, 103, 820,
                 29417, 354, 1, 54, 15, 9, 102, 34, 6, 467, 79, 756, 9, 30, 8, 2572, 10182, 245, 111, 1765,
                 5298, 240, 0, 24, 11, 19, 20, 12, 0, 160, 35, 169, 9, 8, 0, 1267, 3069, 70, 19, 256,
                 2006, 116, 1, 10, 4, 5, 11, 614, 0, 4178, 29, 135, 0, 5, 0, 32163, 3299, 14, 1797, 1343,
                 7904, 194, 2, 7, 0, 1, 2, 3, 0, 46, 217, 27, 0, 0, 0, 219, 752, 12, 8, 26,
                 107, 4, 0, 25, 0, 7, 5, 0, 0, 3, 1, 57, 1, 6, 1, 5, 18, 24, 0, 3)

mRNA_data_test <- matrix(mRNA_values, nrow = 14, byrow = TRUE, dimnames = list(genes, cell_types))
mRNA_data_test_translated <- matrix(mRNA_values, nrow = 14, byrow = TRUE, dimnames = list(genes_translated, cell_types))

###
# DEFINE COLOR PALETTES
###
# "#ffffbf",
aa_colors <- c(Phe = "#7f3b08", Ile = "#b35806", Met = "#e08214", Val = "#fdb863", Ser = "#fee0b6",
               Thr = "#2d004b", Gln = "#542788", Asn = "#8073ac", Pro = "#b2abd2",
               Gly ="#4575b4", Tyr = "#abd9e9", His = "#00441b", Trp = "#1b7837", Cys = "#5aae61",
               Asp = "#a50026", Glu = "#d73027", Arg = "#c51b7d", Lys = "#f1b6da")

gradual_groups_35 <- c("#e08214", "#35978f", "#de77ae", "#7fbc41", "#9970ab", "#d6604d", "#4393c3",
                       "#fdb863", "#80cdc1", "#f1b6da", "#b8e186", "#c2a5cf", "#f4a582", "#92c5de",
                       "#b35806", "#01665e", "#c51b7d", "#4d9221", "#762a83", "#b2182b", "#2166ac",
                       "#fee0b6", "#c7eae5", "#fde0ef", "#e6f5d0", "#e7d4e8", "#fddbc7", "#d1e5f0",
                       "#7f3b08", "#003c30", "#8e0152", "#276419", "#40004b", "#67001f", "#053061")

###
# GENERAL VARIABLES
###
amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
three_letter_amino_acids <- c("Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr")
short_gene_names_mRNA_data <- rownames(mRNA_data_test)
ENSG_gene_names_mRNA_data <- rownames(mRNA_data_test_translated)

assay_map_AA <- list(demand = list(assay_id = list("AADemand"), check_assays = list("CodonUsage"), AAfunction = "sense"),
                     supply = list(assay_id = list("AASupply"), check_assays = list("AnticodonUsage"), AAfunction = "reverse"),
                     both = list(assay_id = list("AADemand", "AASupply"), check_assays = list("CodonUsage", "AnticodonUsage"), AAfunction = list("sense", "reverse")))

assay_map_TE <- list(codon = list(mRNA = "CodonUsage", tRNA = "AnticodonUsage", id = "tTEresults_codon"),
                     aa = list(mRNA = "AADemand", tRNA = "AASupply", id = "tTEresults_AA"))

###
# SAVE THE VARIABLES
###
usethis::use_data(final_matrix_genetic_code,
                  codon_freq_table_canonical_hg38, codon_freq_table_canonical_mm39,
                  mRNA_data_test, mRNA_data_test_translated,
                  aa_colors, gradual_groups_35,
                  amino_acids, three_letter_amino_acids, short_gene_names_mRNA_data, ENSG_gene_names_mRNA_data,
                  assay_map_AA, assay_map_TE,
                  overwrite = TRUE, internal = TRUE)
