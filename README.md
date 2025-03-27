# tTEscanR <a href="https://your-link-here"><img src="https://github.com/user-attachments/assets/213de690-12cd-4fa8-862a-cadbd8872bf4" alt="tTEscanR_logo" align="right" width="120"></a>
<img src="https://img.shields.io/badge/Language-R-blue.svg" style="zoom:100%;" />

## 1. What is tTEscanR?

**tTEscanR** is a versatile and user-friendly **R package** designed to quantify and analyze the relationship between codon usage in mRNA and the availability of corresponding anticodons in tRNA. The package computes a **theoretical translation efficiency (tTE)** score as a proxy of translation efficiency. As a result, both gene expression and chromatin accessibility profile data are required.

The analysis can be carried out across three hierarchical layers of information: gene expression, codon and anticodon pool, and amino acid level. This multi-layered approach provides a comprehensive view of translation efficiency. The package is optimized for both bulk and single-cell datasets.

With its modular structure, the **tTEscanR** package allows to run specific components independently or as part of a comprehensive pipeline, offering flexibility to enhance and complement the analysis of **codon-anticodon dynamics** across various biological contexts.

## 2. Setup

The **tTEscanR** source code is entirely written in R. The latest full distribution release can be downloaded from GitHub:

```{r}
install.packages("/avarassanchez/tTEscanR")
library(tTEscanR)
```

### 2.1. Data specifications

The **input** data is expected to be **already pre-processed** following the steps outlined in [Gao et al., 2022](#5-references) or other desired user-custom pipelines to obtain gene expression count matrices. Both mRNA and tRNA inputs should be given as gene expression count matrices with features as rows and conditions as columns. In this context, features represent mRNA transcripts or tRNA genes. The conditions will differ based on the data source.  In **bulk** datasets, conditions typically represent a combination of model and replicate, whereas in **single-cell** datasets, conditions correspond to specific tissue and cell type designations. For the features **tTEscanR** considers mRNA genes and confidently predicted tRNA genes. For simplicity, some functions have been incorporated into the pre-processing module, as described in the [Helper functions](#35-helper-functions) section.

## 3. Workflow functionalities

![tTEscanR_workflow](https://github.com/user-attachments/assets/2cf26e1c-5e48-40d6-b81b-300e6293c98a)

### 3.1. Defining the tTEscanR object

The **tTEscanR** object is dynamically updated to store matrices and metadata at each analysis step, ensuring efficient tracking and organization of inputs and outputs throughout the pipeline. It can be created using a single dataset and later updated with additional datasets, or all desired datasets can be in a list and passed directly to the main function. 

**Function:** `Create_tTEscanR_Object()`. 

Creates a **tTEscanR** object that will contain an "assays" and "meta.data" slots and in each individual sections as specified in the parameters. 

**Parameters:**
- `counts`: Count matrix (or list of matrices). 
- `assay`: Label (or list of labels) to identify the `counts`. 
- `meta.data`: Additional data to include in the object. 
- `meta.data.ids`: Label(s) to identify the `meta.data`. 
- `verbose`: Logical, if TRUE, displays information messages. 

<hr>

**Function:** `Update_tTEscanR_Object()`. 

Updates an existing **tTEscanR** object by overwriting existing sections or adding new ones, as specified in the parameters.

**Parameters:**
- `object`: The existing tTEscanR object to be updated.
- `counts`: Count matrix (or list of matrices).
- `assay`: Label(s) to identify the `counts`.
- `meta.data`: Additional data  to include in the object.
- `meta.data.ids`: Label(s) to identify the `meta.data`.
- `overwrite.assay`: Logical, if TRUE, overwrites the assays (repeated `assay`) present in the object.
- `overwrite.metadata`: Logical, if TRUE, overwrites the metadata (repeated `meta.data.ids`) present in the object.
- `verbose`: Logical, if TRUE, displays information messages. 

```{r}
data(mRNA_data)
data(tRNA_data)

# Adding the mRNA dataset to the object
tTEobject <- Create_tTEscanR_Object(counts = mRNA_data, assay = "mRNA") 

# Updating the object created before with the tRNA dataset
tTEobject <- Update_tTEscanR_Object(object = tTEobject, counts = tRNA_data, assay = "tRNA", 
                                    meta.data = NULL, meta.data.ids = NULL,
                                    overwrite.assay = FALSE, overwrite.metadata = FALSE)

# If we wanted to add at the same time both datasets we can use the following sintaxis:
tTEobject_paired <- Create_tTEscanR_Object(counts = list(mRNA_data, tRNA_data), 
                                           assay = list("mRNA", "tRNA"))
```

### 3.2. Codon-anticodon usage assessment

**Codon usage** is computed through matrix multiplication between the mRNA gene expression data (genes as rows and conditions as columns) and a **codon frequency per gene** table. This table represents the codon distribution of each protein-coding gene in a reference genome, and includes   protein-coding genes from **Ensembl** formatted as a gene matrix, with the 61 sense codons as rows and each gene as columns. As a result, the codon usage calculation produces a matrix where **codons are represented as rows and conditions as columns**. This transformation enables downstream analysis of codon preferences across different conditions.

**Function:** `ComputeCodonUsage()`. 

Updates the **tTEscanR** object by computing the codon usage of the mRNA data using matrix multiplication with a reference codon frequency table. The function allows for customization using a user-provided codon frequency per gene table (`codon_freq`) or default species-specific tables available in **tTEscanR**. 

**Parameters:**
- `object`: The existing **tTEscanR** object containing a mRNA assay.
- `codon_freq`: User-provided codon frequency per gene table.
- `species`: Label species-specific for a default codon frequency table ("hg38" for human, "mm39" for mouse).
- `filter`: Label to define which transcript to use if several are available for a gene ("canonical" or "length" for the longest).
- `translate`: Logical, if TRUE, changes gene annotation to match between codon frequency table and mRNA assay.
- `reduce`: Integer specifying a factor to divide codon usage values if they exceed R's maximum allowed values.
- `compute_codon_exonic_background`: Logical, if TRUE, computes the codon exonic background.
- `compute_mean_codon_usage`: Logical, if TRUE, computes the mean codon usage across conditions.
- `conditions`:
- `name_sep`:
- `formula`:
- `compute_correlation_background_mean`: Logical, if TRUE, computes the exonic background and meand codon usage correlation.
- `overwrite.assay`: Logical, if TRUE, overwrites the codon usage assay present in the object.
- `overwrite.metadata`: Logical, if TRUE, overwrites the metadata present in the object.
- `verbose`: Logical, if TRUE, displays information messages.

```{r}
# Here we use: 
# - Default human hg38 codon frequency per gene table
# - Canonical setting to filter potential mRNA transcript repetitions
# - Enable to translate the genes if no matching formats are detected

tTEobject <- ComputeCodonUsage(object = tTEobject, 
                               codon_freq = NULL, 
                               species = "hg38", 
                               filter = "canonical",
                               translate = TRUE)
```

<hr>

**Anticodon usage** is computed from a tRNA gene expression matrix including high-confidence tRNA genes predicted using [tRNAscan-SE](https://trna.ucsc.edu/tRNAscan-SE/) as rows, and conditions as columns. Details on how to pre-process the tRNA data described in [Helper functions](#35-helper-functions). To create a more concise representation, tRNA genes sharing the same anticodon are grouped together by **tRNA isotypes**, producing a reduced matrix. 

**Function:** `ComputeAnticodonUsage()`.

Updates the **tTEscanR** object by computing the anticodon usage from the tRNA expression data. It groups tRNA genes that share the same anticodon, generating a reduced expression matrix where tRNA isotypes are represented as rows and conditions as columns.

**Parameters:**
- `object`: The existing **tTEscanR** object containing a tRNA assay.
- `overwrite.assay`: Logical, if TRUE, overwrites the codon usage assay present in the object.
- `overwrite.metadata`: Logical, if TRUE, overwrites the metadata present in the object.
- `verbose`: Logical, if TRUE, displays information messages. 

```{r}
tTEobject <- ComputeAnticodonUsage(object = tTEobject)
```

### 3.3. Amino acid level assessment

The **amino acid level** assessment involves calculating amino acid **demand** and **supply**. Demand is determined by grouping codons that correspond to the same amino acid (based on the genetic code), while supply is measured by pooling anticodon families based on the amino acid they recognize. These calculations can be performed separately using "demand" or "supply" or together" by specifying "both".

**Function:** `ComputeAAUsage()`.

Computes the amino acid usage, either by assessing the `demand` based on codon usage, or the `supply` based on anticodon usage. You can compute both demand and supply simultaenously, depending on the specified action.

**Parameters:**
- `object`: The existing **tTEscanR** object containing a codon usage and/or anticodon usage assay.
- `action`: Label action-specific:
    - "demand", computes amino acid demand from the codon usage matrix.</li>
    - "supply", computes amino acid supply from the anticodon usage matrix.</li>
    - "both", computes both amino acid demand and supply.</li>
- `overwrite.assay`: Logical, if TRUE, overwrites the codon usage assay present in the object.
- `overwrite.metadata`: Logical, if TRUE, overwrites the metadata present in the object.
- `verbose`: Logical, if TRUE, displays information messages.

```{r}
tTEobject <- ComputeAAUsage(object = tTEobject, 
                            action = "demand")

tTEobject <- ComputeAAUsage(object = tTEobject, 
                            action = "supply")

# If we wanted to compute the AA demand and supply simultaneously
tTEobject <- ComputeAAUsage(object = tTEobject, 
                            action = "both", 
                            overwrite.assay = TRUE, 
                            overwrite.metadata = TRUE)
```
  
### 3.4. tTE computation

The **Theoretical Translation Efficiency (tTE)** is assessed by calculating Spearman's rank correlation coefficient between amino acid demand, derived from mRNA codon usage, and amino acid supply, determined by tRNA anticodon usage. Ensuring **matching conditions** between both datasets is crucial for this step, as only shared conditions are retained for downstream analysis. Consequently, the assessment is highly sensitive to annotation variability. 

**Function:** `Compute_tTE()`

Computes the tTE score across matching conditions at the codon-anticodon usage level and/or amino acid demand and supply level.

**Parameters:**
- `object`: The existing **tTEscanR** object containing both codon-anticodon or amino acid demand-supply assays.
- `conditions`:
- `name_sep`:
- `formula`:

```{r}
tTEobject <- Compute_tTE(object = tTEobject, 
                         conditions = c("tissue", "cell.type"),
                         name_sep = "-",
                         formula = ~ tissue)
```

### 3.5. Helper functions

**tTEscanR** provides supplementary functions designed to enhance and streamline analysis. These functions are organized into modules and can be used independently or as part of the core functions described above.

-------- **Pre-processing module** --------

The tRNAscanR pipeline requires pre-processed files, as outlined in [Gao et al., 2022](#5-references). However, this module includes several functions designed to assist the user with pre-processing tasks. 
- For **tRNA data**, the pre-processing steps necessary for running tTEscanR involve: (i) create of a tRNA expression matrix, (ii) filter out low tRNA cuts, and (iii) prediction of high-confidence tRNA genes using [tRNAscan-SE](https://trna.ucsc.edu/tRNAscan-SE/).
- For **mRNA data**, the funcitons focus on: (i) trimming protein-coding genes (or any user-defined list of targeted genes), and (ii) ensuring consistent gene annotations across different nomenclature formats.

**Function:** `Get_tRNAexpressionMatrix()`

Extracts tRNA expression data from a Chromatin Assay or a Seurat object and converts it into a structured expression matrix. Optionally, the output can be saved to a specified directory.

**Parameters:**
- `chrom`: ChromatinAssay or Seurat object containing tRNA data.
- `assay`: The name of the assay to extract if `chrom` is a Seurat object.
- `tRNA_annotations` List of tRNA gene annotations.
- `species`: The reference genome version of the species, required if tRNA_annotations is not provided.
- `name_sep`: The delimeter used to format tRNA gene names in the output matrix.
- `save`: Logical, if TRUE, saves the tRNA expression matrix into a .rds file.
- `out_name`: String parameter specifying the output name (if save is TRUE).
- `out_directory`: String parameter specifying the output directory name if save is TRUE).
- `verbose`: Logical, if TRUE, displays information messages.

```{r}
# Here we use: 
# - The Seurat object obtained after integrating the fragment files to the scATAC-seq peak matrix
#   - The file has been extracted from Gao et al., 2022
# - Default human hg38 tRNA annotation file
# - The output matrix will not be saved

chromatin_object <- readRDS()
tRNA_expression_matrix <- Get_tRNAexpressionMatrix(chrom = chromatin_object,
                                                   assay = "peaks",
                                                   species = "hg38",
                                                   name_sep = c("-", "-"),
                                                   save = FALSE)
# The generated tRNA_expression_matrix corresponds to the tRNA_data available in tTEscanR
```

<hr>

**Function:** `tRNACutsFilter()`

Removes conditions where tRNA gene expression falls below a specified threshold. It helps ensure that only conditions with sufficient tRNA counts are included in downstream analyses. 

**Parameters:**
- `tRNA_data`: The tRNA gene expression data containing tRNA cuts as rows and conditions as columns.
- `cutoff`: Integer specifying the minimum expression level required for a condition to be retained.

```{r}
tRNA_data_filtered <- tRNACutsFilter(tRNA_data = tRNA_data, cutoff = 5000)
```

<hr>

**Function:** `TranslateGeneName()`

Verifies and, if necessary, translates gene annotations across multiple vectors or data frames, from Ensembl IDs to gene names, or vice-versa. Enables to maintan **consistent annotations** throught the analysis.

**Parameters:** 
- `data_to_translate`: Vector or dataframe containing gene names to be translated.
- `translator_table`: User-provided reference table for gene translation (first column: Ensembl IDs, second column: gene symbols).
- `species`: The reference genome version of the species (required if translator_table is not provided).
- `position`: Location of the genes ("row" or "column") if data_to_translate is a dataframe.
- `notation`:  The gene annotation format ("symbol" or "id") to resolve inconsistencies in data_to_translate.
- `verbose`: Logical, if TRUE, displays information messages.

```{r}
translated_mRNA_genes <- TranslateGeneName(data_to_translate = mRNA_data, species = "hg38", position = "row", notation = "id")
```

-------- **Codon frequency table module** --------

To complement the codon usage analysis, we have developed an independent module to generate the **codon frequency per gene matrix** for reference organisms. This matrix can be obtained using Ensembl (through **biomaRt**) to access the reference sequence of an organism or by directly providing a reference genome or sequence files of interest. This flexibility allows users to retrieve codon freqeuncy per gene matrices for any organim available in Ensembl, enabling broad applicability of the computations. However, it is worth mentioning that accessing Ensembl may occasionally result in connection errors due to platform-related issues unrelated to **tTEscanR**. To facilitate a smoother analysis, we have included pre-computed codon frequency per gene matrices for human ("hg38") and mouse ("mm39") as default options. 

**Function:** `ObtainCodonComposition()`

Based on a Ensembl dataset name retrieves the reference genome of the organism and computes the codon frequency per gene matrix. This function can also be used in a targeted mode by using the "transcripts" parameter to give a vector of gene ids to retrieve. 

**Parameters:**
- `dataset_name`: String specifying the Ensembl dataset name for the species.
- `genome_file`: File containing the full DNA sequence of a genome.
- `genes_file`: File containing gene IDs along with their corresponding sequences.
- `transcripts`: Vector of genes to subset the analyses. 
- `filter`: String defining the criteria ("canonical" or "length") to select a transcript if several are available for the same gene.
- `verbose`: Logical, if TRUE, displays information messages.

```{r}
# To check the list of Ensembl dataset names
datasets <- biomaRt::listDatasets(useEnsembl(biomart = "ensembl"))
datasets$dataset

# Retrieving the codon frequency per gene matrix from the human reference genome
human_codon_freq_per_gene_matrix <- ObtainCodonFreqPerGene(dataset_name = "hsapiens_gene_ensembl", filter = "canonical")

# Using a targeted approach to get the codon frequency of the genes included in the mRNA data 
targeted_genes <- c("ENSG00000059588", "ENSG00000052841", "ENSG00000173153", "ENSG00000058799", "ENSG00000071203")
targeted_genes_codon_freq_per_gene_matrix <- ObtainCodonFreqPerGene(dataset_name = "hsapiens_gene_ensembl", transcripts = targeted_genes)
```

<hr>

**Function:** `ExtractCodonComposition()`

Computes the codon composition of a given set of DNA sequences, analyzing the frequency of codons within each sequence.  

**Parameters:**
- `sequences`: A vector or list of DNA sequences to analyze.  
- `verbose`: Logical, if TRUE, displays information messages. 

```{r}
codon_composition <- extract_codon_composition(c("ATGCGTACG", "TTAAGGCCG"))
```

-------- **Differential expression analysis module** --------

The differential expression analysis module provides a comprehensive ser of functions for identifying and analyzing differentially expressed features across multiple biological levels, including gene expression, codon-anticodon interactions, and amino acid usage. It supports various statistical frameworks, normalization techniques, and visualization tools to ensure accurate and comprehensive analysis.

**Function:** `DESeq2runner()`

Performs differential expression analysis using the DESeq2 package and provides several output visualizations (heatmap, PCA plot and/or volcano plot), as well as normalizes the input data and gives the corrected data object. 

**Parameters:**
- `data`: Matrix with features as rows and conditions as columns.
- `conditions`: A vector with the labels for each section of the conditions labels.
- `name_sep`: Specification of the notation separation used in the condition labels to run a DESeq2 analysis.
- `formula`: Design formula to use for the DESeq2 analysis. The names need to be consistent with the columns in `data`.
- `targets`: Dataset with two columns: (i)conditions to select, and (ii) labels to use for the comparisons.
- `fc_threshold`: Fold change threshold to be used in the volcano plot (applicable for targeted approach).
- `pval_threshold`: P-value threshold to be used in the volcano plot (applicable for targeted approach).
- `reduce`:  Integer specifying a factor to divide codon usage values if they exceed R's maximum allowed values.
- `heatmap`: Logical, if TRUE, a heatmap should will be displayed (not applicable for targeted approach).
- `PCA`: Logical, if TRUE, the principal component analysis will be computed (not applicable for targeted approach).
- `numPC`: Number of principal components to consider in the PCA analysis (if PCA is TRUE).
- `color_factor`: Column name of `data` to use to define the groups (if PCA is TRUE.
- `labels Boolean variable to indicate if labels should be included in the PCA plots (if \code{PCA} is TRUE).
- `verbose`: Logical, if TRUE, displays information messages. 

```{r}
differential_expression_mRNA_analysis <- DESeq2runner(data = mRNA_data,
                                                      conditions = c("tissue", "cell.type"),
                                                      name_sep = "-",
                                                      formula = ~ tissue,
                                                      heatmap = TRUE,
                                                      PCA = TRUE,
                                                      numPC = 5,
                                                      color_factor = "tissue",
                                                      labels = FALSE)
# List of outputs: (i) Heatmap, (ii) PCA plot, and (iii) Size corrected data.
# Each element stored in the list of outputs can be accessed using '$'.

# Targeted approach
targets_neuron <- data.frame(search = c("neuron", "ENS neurons"),
                                      class = c("neuron", "other"))
differential_expression_mRNA_analysis_targeted <- DESeq2runner(data = mRNA_data,
                                                               conditions = c("tissue", "cell.type"),
                                                               name_sep = "-",
                                                               formula = ~ class,
                                                               targets = targets_neuron,
                                                               fc_threshold = 1,
                                                               pval_threshold = 0.05)
```


## 4. Usage

```{r}
library(tTEscanR)
data(mRNA_data)
data(tRNA_data)

# Adding the mRNA dataset to the object
tTEobject <- Create_tTEscanR_Object(counts = mRNA_data, assay = "mRNA") 

# Updating the object created before with the tRNA dataset
tTEobject <- Update_tTEscanR_Object(object = tTEobject, counts = tRNA_data, assay = "tRNA", 
                                    meta.data = NULL, meta.data.ids = NULL,
                                    overwrite.assay = FALSE, overwrite.metadata = FALSE)

# Adding metadata to the object
matching_celltypes <- intersect(colnames(mRNA_data), colnames(tRNA_data)) # We are defining a vector of strings
tTEobject <- Update_tTEscanR_Object(object = tTEobject, 
                                    meta.data = list(matching_celltypes), 
                                    meta.data.ids = list("matching_celltypes"), 
                                    overwrite.metadata = TRUE)

# Here we use: 
# - Default human hg38 codon frequency per gene table
# - Canonical setting to filter potential mRNA transcript repetitions
# - Enable to translate the genes if no matching formats are detected
tTEobject <- ComputeCodonUsage(object = tTEobject, 
                               codon_freq = NULL, 
                               species = "hg38", 
                               filter = "canonical",
                               translate = TRUE)

tTEobject <- ComputeAnticodonUsage(object = tTEobject)

# The both parameter allows to perform at the same time the AA demand and supply assessment
tTEobject <- ComputeAAUsage(object = tTEobject, action = "both")

tTEobject <- Compute_tTE(object = tTEobject, 
                         conditions = c("tissue", "cell.type"),
                         name_sep = "-",
                         formula = ~ tissue)
```

## 5. References

Gao W, Gallardo-Dodd CJ, Kutter C. *Cell type-specific analysis by single-cell profiling identifies a stable mammalian tRNA-mRNA interface and increased translation efficiency in neurons.* Genome Res. 2022;32(1):97-110. [doi:10.1101/gr.275944.121](https://doi.org/10.1101/gr.275944.121)
