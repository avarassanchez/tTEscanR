# tTEscanR <a href="https://your-link-here"><img src="https://github.com/user-attachments/assets/213de690-12cd-4fa8-862a-cadbd8872bf4" alt="tTEscanR_logo" align="right" width="120"></a>
<img src="https://img.shields.io/badge/Language-R-blue.svg" style="zoom:100%;" />

## 1. What is tTEscanR?

**tTEscanR** is a versatile and user-friendly **R package** designed to quantify and analyze the relationship between codon usage in mRNA and the availability of corresponding anticodons in tRNA. The package computes a **theoretical translation efficiency (tTE)** score as a proxy of translation efficiency. As a result, both gene expression and chromatin accessibility profile data are required.

The analysis can be carried out across three hierarchical layers of informaiton: gene expression, codon and anticodon pool, and amino acid levle. This multi-layered approach provides a comprehensive view of translation efficiency. The package is optimized for both bulk and single-cell datasets.

With its modular structure, the **tTEscanR** package allows to run specific compoennts independently or as part of a comprehensive pipeline, offering flexibility to enhance and complement the analysis of **codon-anticodon dynamics** across various biological contexts.

## 2. Setup

The **tTEscanR** source code is entirely written in R. The latest full distribution release can be downloaded from GitHub:

```{r}
install.packages("(avarassanchez/tTEscanR")
library(tTEscanR)
```

### 2.1. Data specifications

The **input** data is expected to be **already pre-processed** following the steps outlined in [Add reference] or other desired user-custom pipelines to obtain gene expression count matrices. Both mRNA and tRNA inputs should be given as gene expression count matrices with features as rows and conditions as columns. In this context, features represent mRNA transcripts or tRNA genes. The conditions will differ based on the data source.  In **bulk** datasets, conditions typically represent a combination of model and replicate, whereas in **single-cell** datasets, conditions correspond to specific tissue and cell type designations. For the features **tTEscanR** considers mRNA genes and confidently predicted tRNA genes.

Further details can be found in the article referenced.

## 3. Workflow functionalities

![tTEscanR_workflow](https://github.com/user-attachments/assets/b74af6d7-3d89-46df-8885-98104b448d36)

### 3.1. Defining the tTEscanR object

The **tTEscanR** object is dynamically updated to store matrices and metadata at each analysis step, ensuring efficient tracking and organization of inputs and outputs throughout the pipeline. It can be created using a single dataset and later updated with additional datasets, or all desired datasets can be in a list and passed directly to the main function. 

**Function:** `Create_tTEscanR_Object()`. 

Creates a **tTEscanR** object that will contain an "assays" and "meta.data" slot and in each individual sections as specified in the parameters. 
**Parameters:**
- `counts`: Count matrix (or list of matrices). 
- `assay`: Label(s) to identify the `counts`. 
- `meta.data`: Additional data to include in the object. 
- `meta.data.ids`: Label(s) to identify the `meta.data`. 
- `verbose`: Logical, if TRUE, displays information messages. 
    
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

**Codon usage** is computed through matrix multiplication between the mRNA data (genes as rows and conditions as columns) and a reference **codon frequency per gene** table. This reference table consists of **Ensembl** protein-coding genes formatted as a gene matrix, with the 61 sense codons as rows and each  gene as column. As a result, the codon usage calculation produces a matrix where **codons are represented as rows and conditions as columns**. This transformation enables downstream analysis of codon preferences across different conditions.

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

**Anticodon usage** is computed from a tRNA gene expression matrix including high-confidence tRNA genes predicted using tRNAscan-SE [Addd reference] as rows, and conditions as columns. To create a more concise representation, tRNA genes sharing the same anticodon are grouped together by **tRNA isotypes**, producing a reduced matrix. 

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

Description: Computes the amino acid usage, either by assessing the `demand` based on codon usage, or the `supply` based on anticodon usage. You can compute both demand and supply simultaenously, depending on the specified action.
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

Description: Computes the tTE score across matching conditions at the codon-anticodon usage level and/or amino acid demand and supply level.

```{r}
tTEobject <- Compute_tTE(object = tTEobject, 
                         conditions = c("tissue", "cell.type"),
                         name_sep = "-",
                         formula = ~ tissue)
```

### 3.5. Helper functions

**tTEscanR** provides supplementary functions designed to enhance and streamline analysis. These functions are organized into modules and can be used independently or as part of the core functions described above.

mRNA transcripts are typically annotated using either Ensembl IDs or gene names. To maintain **consistent annotation** throughout the analysis, a dedicated function has been developed to verify and, if necessary, translate gene annotations across multiple vectors or data frames. This function is particularly important for codon usage assessment, where genes must be consistently annotated to enable matrix multiplication between the codon frequency per gene table and the mRNA gene expression count matrix.

<hr>

**Function:** `TranslateGeneName()`

Description: Translates genes from the Ensembl id annotation to the gene name format, and vice-versa.
Parameters: data_to_translate, translator_table, species, position, notation, verbose

```{r}
transalted_mRNA_genes <- TranslateGeneName(data_to_translate = mRNA_data, species = "hg38", position = "row", notation = "id")
```

<hr>

To complement the codon usage analysis we have developed an independent module to retrieve the **codon frequency per gene matrix** from reference organisms. This table can be computed using Ensembl (though biomaRt) to access the reference sequence of an organisms or by directly imputing the reference genome or sequence of interest files. It is worth mentioning that the approach using Ensembl can give connection errors to the platform itself, unrelated to **tTEscanR**. In order to enable a more straight forward analysis we have incorporated as default the human ("hg38") and mouse ("mm39") codon frequency per gene matrices. 

**Function:** `ObtainCodonComposition()`

Description: Based on a Ensembl dataset name retrieves the reference genome of the organism and computes the codon frequency per gene matrix. This function can also be used in a targeted mode by using the "transcripts" parameter to give a vector of gene ids to retrieve. 

```{r}
# To check the list of Ensembl dataset names
datasets <- biomaRt::listDatasets(useEnsembl(biomart = "ensembl"))
datasets$dataset

# Retrieving the codon frequency per gene matrix from the human reference genome
human_codon_freq_per_gene_matrix <- ObtainCodonFreqPerGene(dataset_name = "hsapiens_gene_ensembl")

# Using a targeted approach to get the codon frequency of the genes included in the mRNA data 
targeted_genes <- head(rownames(mRNA_data))
targeted_genes_codon_freq_per_gene_matrix <- ObtainCodonFreqPerGene(dataset_name = "hsapiens_gene_ensembl", transcripts = targeted_genes)
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
