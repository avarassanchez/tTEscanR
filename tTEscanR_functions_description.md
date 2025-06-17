# tTEscanR - Description of the funcitons <a href="https://your-link-here"><img src="https://github.com/user-attachments/assets/8cca530a-f0a5-4284-bf2e-cf030a8193fa" alt="logo_tTEscanR" align="right" width="120"></a>
<img src="https://img.shields.io/badge/Language-R-blue.svg" style="zoom:100%;" />

**tTEscanR** has a modular structure that allows to run specific components independently or as part of a comprehensive pipline, offering flexibility to enhance and complement the analysis of **codon-anticodon dynamics** across various biological contexts.

## 1. Core functions

### 1.1. Defining the tTEscanR object

The **tTEscanR** object is dynamically updated to store assays and metadata at each analysis step, ensuring efficient tracking and organization of inputs and outputs throughout the pipeline. 

**Function:** `Create_tTEscanR_Object()`. 

Initializes a **tTEscanR** object, a strucutred container designed to hold data. The object consists of two main components; *"assays"*, which stores data matrices, and *"meta.data"*, which stores associated information and additional intermediate calculations.  A **tTEscanR** object can be created using a single dataset, or initialized with a list of datasets provided at once. Additional assays and metadata layers can be appended later using `Update_tTEscanR_Object()`. 

**Parameters:**
- `counts`: A count matrix (or list of matrices) to be stored in the *"assays"* slot. 
- `assay`: A character string (or vector of characters) to identify the name of the `counts`. 
- `meta.data`: Optional; a list with additional information that to be stored in the *"meta.data"* slot.
- `meta.data.ids`: Optional; a list of labels to identify the `meta.data`. 
- `verbose`: Logical; if TRUE, displays information messages. 

<hr>

**Function:** `Update_tTEscanR_Object()`. 

Updates an existing **tTEscanR** object by overwriting existing sections or adding new ones, as specified in the parameters.

**Parameters:**
- `object`: An existing **tTEscanR** object to be updated.
- `counts`: Optional; a count matrix (or list of matrices) to be stored in the *"assays"* slot.
- `assay`:  Optional; a character string to identify the name of the `counts`.
- `meta.data`: Optional; a list with additional information that to be stored in the *"meta.data"* slot.
- `meta.data.ids`: Optional; a list of labels to identify the `meta.data`. 
- `overwrite.assay`: Logical; if TRUE, overwrites any existing *"assays"* in the `object` if the `assay` label coincides.
- `overwrite.metadata`: Logical; if TRUE, overwrites any existing *"meta.data"* in the `object` if the `meta.data.ids` label coincides.
- `verbose`: Logical; if TRUE, displays information messages.

```{r}
data(mRNA_data)
data(tRNA_data)

# Adding the mRNA dataset to the object
tTEobject <- Create_tTEscanR_Object(counts = mRNA_data, assay = "mRNA") 

# Updating the object created before with the tRNA dataset
tTEobject <- Update_tTEscanR_Object(object = tTEobject, counts = tRNA_data, assay = "tRNA", 
                                    meta.data = NULL, meta.data.ids = NULL,
                                    overwrite.assay = FALSE, overwrite.metadata = FALSE)

# If we want to add at the same time both datasets we can use the following sintaxis:
tTEobject_paired <- Create_tTEscanR_Object(counts = list(mRNA_data, tRNA_data), 
                                           assay = list("mRNA", "tRNA"))
```

### 1.2. Codon-anticodon usage assessment

The **codon usage** is computed through matrix multiplication between the mRNA gene expression data (genes as rows and conditions as columns) and a **codon frequency-per-gene** table. This table represents the codon distribution of each protein-coding gene in a reference genome, with the 61 sense codons as rows and each gene as columns. As a result, the codon usage calculation produces a matrix where **codons are represented as rows and conditions as columns**. This transformation enables downstream analysis of codon preferences across different conditions.

**Function:** `ComputeCodonUsage()`. 

Estimates codon usage profiles based on gene-level mRNA expression data. It optionally accepts pre-computed user-prodived codon-frequency tables, or uses internally generated tables when not provided. Currently the default organisms are human (hg38) and mouse (mm39). When enabled, it can evaluate the correlation between background codon composition and observed mean codon usage. 

**Parameters:**
- `object`: An existing **tTEscanR** object containing a mRNA assay.
- `codon_freq`: Optional; a user-provided codon frequency-per-gene table.
- `translate`: Logical; if TRUE, resolves inconsistent gene annotations.
- `species`: Either *"hg38"* (human) or *"mm39"* (mouse) to specify which default codon frequency-per-gene table to use. Required if `codon_freq` is not provided or if `translate` is enabled.
- `filter`: Either *"canonical"* (default) or *"length"* (longest transcript) to specify which transcript to choose if several are available for the same gene.
- `reduce`: Numeric; a scaling factor used to normalize large expression values that exceed R's handling capacity.
- `compute_codon_exonic_background`: Logical; if TRUE, computes the codon exonic background.
- `compute_mean_codon_usage`: Logical; if TRUE, computes the mean codon usage across conditions.
- `conditions`: A character vector specifying the labels' structure of the conditions (columns in the matrix). Required if `compute_mean_codon_usage` is enabled.
- `name_sep`: A string delimiter used in `conditions` to separate each part of the labels. Required if `compute_mean_codon_usage` is enabled.
- `compute_correlation_background_mean`: Logical; if TRUE, computes the exonic background and meand codon usage correlation. Requires `compute_codon_exonic_background` and `compute_mean_codon_usage` to be enabled.
- `overwrite.assay`: Logical; if TRUE, overwrites any existing *"assays"* in the `object`.
- `overwrite.metadata`: Logical; if TRUE, overwrites any existing *"meta.data"* in the `object`.
- `verbose`: Logical; if TRUE, displays information messages.

```{r}
# Here we use: 
# - Default human hg38 codon frequency per gene table
# - Canonical setting to filter potential mRNA transcript repetitions
# - Enable to translate the genes if no matching formats are detected

tTEobject <- ComputeCodonUsage(object = tTEobject, 
                               codon_freq = NULL, 
                               species = "hg38", 
                               filter = "canonical",
                               translate = TRUE,
                               conditions = c("tissue", "cell.type"),
                               name_sep = "-")
```

<hr>

The **anticodon usage** is computed from a tRNA gene expression matrix including high-confidence tRNA genes predicted using [tRNAscan-SE](https://trna.ucsc.edu/tRNAscan-SE/) as rows, and conditions as columns. Details on how to pre-process the tRNA data are described in [Helper functions](#2-helper-functions). To create a more concise representation, tRNA genes sharing the same anticodon are grouped together by **tRNA isoacceptors**, producing a reduced matrix. 

**Function:** `ComputeAnticodonUsage()`.

Calculates anticodon usage profiles from tRNA gene expression data. It summarizes the expression of tRNAs by their anticodon isdentity, which can be used to estimate the tRNA supply landscape. As a result it produces a matrix where **anticodons are represented as rows and condition as columns**.

**Parameters:**
- `object`: An existing **tTEscanR** object containing a tRNA assay.
- `overwrite.assay`: Logical; if TRUE, overwrites any existing *"assays"* in the `object`.
- `overwrite.metadata`: Logical; if TRUE, overwrites any existing *"meta.data"* in the `object`.
- `verbose`: Logical; if TRUE, displays information messages.

```{r}
tTEobject <- ComputeAnticodonUsage(object = tTEobject)
```

### 1.3. Amino acid level assessment

The **amino acid level** assessment involves calculating amino acid **demand** and **supply**. Demand is determined by grouping codons that correspond to the same amino acid (based on the genetic code), while supply is measured by pooling anticodon families based on the amino acid they recognize. These calculations can be performed separately using "demand" or "supply" or together" by specifying "both".

**Function:** `ComputeAAUsage()`.

Computes the amino acid (AA) demand and/or supply from a given codon or anticodon usage matrix, respectively. It extracts codon and/or anticodon usage matrices and aggregates their contributions based on the standard genetic code, mapping each codon or anticodon to its corresponding amino acid. The resulting values reflect total usage (demand) or availability (supply) of each AA, depending on the input type.

**Parameters:**
- `object`: An existing **tTEscanR** object containing a codon and/or anticodon usage assay.
- `level`: Either *"demand"*, *"supply"* or *"both"* to indicate which analysis to perform.
    - "demand", computes amino acid demand from the codon usage matrix.</li>
    - "supply", computes amino acid supply from the anticodon usage matrix.</li>
    - "both", computes both amino acid demand and supply.</li>
- `overwrite.assay`: Logical; if TRUE, overwrites any existing *"assays"* in the `object`. 
- `verbose`: Logical; if TRUE, displays information messages.

```{r}
tTEobject <- ComputeAAUsage(object = tTEobject, 
                            level = "demand")

tTEobject <- ComputeAAUsage(object = tTEobject, 
                            level = "supply")

# If we want to compute the AA demand and supply simultaneously
tTEobject <- ComputeAAUsage(object = tTEobject, 
                            level = "both", 
                            overwrite.assay = TRUE)
```

### 1.4. tTE computation

The **Theoretical Translation Efficiency (tTE)** is assessed by calculating Spearman's rank correlation coefficient between amino acid demand, derived from mRNA codon usage, and amino acid supply, determined by tRNA anticodon usage. Ensuring **matching conditions** between both datasets is crucial for this step, as only shared conditions are retained for downstream analysis. Consequently, the assessment is highly sensitive to annotation variability. 

**Function:** `Compute_tTE()`.

Computes the tTE score across matching conditions at the codon-anticodon usage level and/or amino acid demand and supply level.

**Parameters:**
- `object`: An existing **tTEscanR** object containing both codon-anticodon and/or amino acid demand-supply assays.
- `level`: Either *"codon"*, *"AA"* or *"both"* to specify the level at which to compute the tTE score.
    - "codon", computes tTE from codon-anticodon usage.</li>
    - "AA", computes tTE from amino acid demand and supply.</li>
    - "both", computes tTE from codon-anticodon usage and amino acid demand and supply.</li>
- `conditions`: A character vector specifying the labels' structure of the conditions (columns in the matrix). Required if `compute_mean_codon_usage` is enabled.
- `name_sep`: A string delimiter used in `conditions` to separate each part of the labels. Required if `compute_mean_codon_usage` is enabled.
- `overwrite`: Logical; if TRUE, overwrites any existing *"meta.data"* in the `object`.
- `verbose`: Logical; if TRUE, displays information messages.

```{r}
tTEobject <- Compute_tTE(object = tTEobject,
                         level = "AA", 
                         conditions = c("tissue", "cell.type"),
                         name_sep = "-")
```

## 2. Helper functions

### 2.1. Pre-processing module

**Function:** `Get_tRNAexpressionMatrix()`.

Extracts tRNA expression data from a Chromatin Assay or a Seurat object and converts it into a structured expression matrix. Optionally, the output can be saved into a specified directory.

**Parameters:**
- `chrom`: A ChromatinAssay or Seurat object containing tRNA data.
- `assay`: Optional; q character string specifying the name of the assay to retrive from `chrom`it if is a Seurat object.
- `tRNA_annotations`: A list of tRNA gene annotations.
- `species`:  Either *"hg38"* (human) or *"mm39"* (mouse) to specify which default codon frequency-per-gene table to use. Required if `tRNA_annotations` is not provided.
- `name_sep`: A string delimiter to format the tRNA gene names in the output matrix.
- `save`: Logical, if TRUE, saves the tRNA expression matrix into a *".rds"* file.
- `out_name`: Optional; name for the output file. Required is `save` is enabled.
- `out_directory`: Optional; path to the directory where the output file will be saved. Required if `save` is enabled.
- `verbose`: Logical; if TRUE, displays information messages.

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

**Function:** `tRNACutsFilter()`.

Filters a tRNA expression matrix by removing tRNA genes (rows) that fall below a specific total read count, specified as a cutoff. It is useful to eliminate low-quality or poorly seqeunced conditions that may bias downstream analyses. 

**Parameters:**
- `tRNA_data`: The tRNA gene expression data containing tRNA cuts as rows and conditions as columns.
- `cutoff`: Minimum expression level required for a condition to be retained.

```{r}
tRNA_data_filtered <- tRNACutsFilter(tRNA_data = tRNA_data, cutoff = 5000)
```

<hr>

**Function:** `TranslateNames()`.

Verifies and, if necessary, translates gene annotations across multiple vectors or data frames, from Ensembl IDs to gene names, or vice-versa, ensuring **consistent annotation** throughout the analysis. Additionally, it enables bidirectional translation between codon and anticodons using the genetic code table, allowing assingment of their corresponding amino acids and supporting hierrchical annotation across levels.

**Parameters:** 
- `data_to_translate`: A vector or dataframe containing gene names or features (codons or anticodons) to be translated.
- `mode`: Either *"gene* or *"feature"* to specify the conversion to perform.
- `translator_table`: Optional; a two-column (Ensembl id and gene name) user-provided reference table for gene translation. Required if `mode` is *"gene"*.
- `species`: Either *"hg38"* (human) or *"mm39"* (mouse) to specify which default gene translator table to use. Required if `mode` is *"gene"*.
- `position`: Either *"row"*, *"column"* or *"column_name"* to specify the location of the genes or features in `data_to_translate`. Required if `data_to_translate` is a data frame. 
- `notation`: Either *"symbol"* or *"id"* for the genes, and *"codon"*, *"anticodon"* or *"AA"* for the features, to indicate the translation format to implement.
- `verbose`: Logical; if TRUE, displays information messages.

```{r}
translated_mRNA_genes <- TranslateNames(data_to_translate = mRNA_data,
                                        mode = "gene, 
                                        species = "hg38",
                                        position = "row",
                                        notation = "id")
```

### 2.2 Codon frequency table module

**Function:** `ObtainCodonComposition()`.

Based on an Ensembl reference genome or a user-provided gene sequence file, this function computes a codon frequency-per-gene matrix. It can optionally subset the analysis to specific transcripts and apply filtering criteria when multiple transcripts are available per gene.

**Parameters:**
- `dataset_name`: A character string specifying the Ensembl species dataset name.
- `genes_file`: A file path to a table containign gene IDs and their corresponding nucleotide sequences.
- `transcripts`: Optional; a vector of transcripts or gene IDs to subset the analysis.
- `filter`:  Either *"canonical"* (default) or *"length"* (longest transcript) to specify which transcript to choose if several are available for the same gene.
- `verbose`: Logical; if TRUE, displays information messages.

```{r}
# To check the list of Ensembl dataset names (datasets$dataset)
datasets <- biomaRt::listDatasets(useEnsembl(biomart = "ensembl"))

# Retrieving the codon frequency per gene matrix from the human reference genome
human_codon_freq_per_gene_matrix <- ObtainCodonFreqPerGene(dataset_name = "hsapiens_gene_ensembl",
                                                           filter = "canonical")

# Using a targeted approach to get the codon frequency of the genes included in the mRNA data 
targeted_genes <- c("ENSG00000059588", "ENSG00000052841", "ENSG00000173153",
                    "ENSG00000058799", "ENSG00000071203")

targeted_genes_codon_freq_per_gene_matrix <- ObtainCodonFreqPerGene(dataset_name = "hsapiens_gene_ensembl",
                                                                    transcripts = targeted_genes)
```

<hr>

**Function:** `ExtractCodonComposition()`.

Analyzes a given set of nucleotide sequences and computes the count of each codon present. 

**Parameters:**
- `sequences`: A list of nucleotide sequences from which to extract the codon composition.  
- `verbose`: Logical, if TRUE, displays information messages. 

```{r}
codon_composition <- extract_codon_composition(c("ATGCGTACG", "TTAAGGCCG"))
```

### 2.3. Differential expression analysis module

**Function:** `ExecuteDESeq2runner()`.

Performs differential expression analysis using the DESeq2 framework on a matrix or list of matrices of expression values. It supports both explorstory visualizations (heatmap and PCA) and targeted comparisons using a custom contrast table.

**Parameters:**
- `list_data`: A list of matrices with features as rows and conditions as columns.
- `conditions`: A character vector specifying the labels' structure of the conditions (columns in the matrix).
- `name_sep`: A string delimiter used in `conditions` to separate each part of the labels.
- `targets`: Optional; a data frame with one column indicating the conditions to select, and another column with the labels for comparisons.
- `fc_threshold`: Numeric; fold change threshold used for highlighting significant features in the volcano plot. Required if `targets` is specified. 
- `pval_threshold`: Numeric; p-value threshold used for highlighting significant features in the volcano plot. Required if `targets` is specified.
- `reduce`: Numeric; a scaling factor used to normalize large expression values that exceed R's handling capacity.
- `heatmap`: Logical, if TRUE, generates a heatmap for exploratory analysis. Not applicable if `targets` is specified.
- `PCA`: Logical, if TRUE, generates a principal component analysis. Not applicable if `targets` is specified. 
- `numPC`: Numeric; number of principal components to consider in the PCA analysis. Required if `PCA` is enabled.
- `color_factor`: A factor based on `conditions` to define the colors in the PCA plot. Required if `PCA` is enabled.
- `labels`: Logical; if TRUE, displays the data points labels in PCA plot. Required if `PCA` is enabled.
- `verbose`: Logical; if TRUE, displays information messages.

```{r}
datasets_to_analyze <- list(mRNA = mRNA_data, tRNA = tRNA_data)
differential_expression_mRNA_analysis <- ExecuteDESeq2runner(list_data = datasets_to_analyze,
                                                             conditions = c("tissue", "cell.type"),
                                                             name_sep = "-",
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
differential_expression_mRNA_analysis_targeted <- ExecuteDESeq2runner(data = list(target_mRNA = mRNA_data),
                                                                      conditions = c("tissue", "cell.type"),
                                                                      name_sep = "-",
                                                                      targets = targets_neuron,
                                                                      fc_threshold = 1,
                                                                      pval_threshold = 0.05)
```
