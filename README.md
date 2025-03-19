# tTEscanR <p align="right"><img width="120" alt="tTEscanR_logo" src="https://github.com/user-attachments/assets/213de690-12cd-4fa8-862a-cadbd8872bf4"></p> 

## 1. What is tTEscanR?

**tTEscanR** is a versatile and user-friendly **R package** designed to quantify and analyze the relationship between codon usage in mRNA and the availability of corresponding anticodons in tRNA. The package computes a **theoretical translation efficiency (tTE)** score as a proxy of translation efficiency. As a result, both gene expression and chromatin accessibility profile data are required.

The analysis can be carried out across three hierarchical layers of informaiton: gene expression, codon and anticodon pool, and amino acid levle. This multi-layered approach provides a comprehensive view of translation efficiency. The package is optimized for both bulk and single-cell datasets.

With its modular structure, the **tTEscanR** package allows to run specific compoennts independently or as part of a comprehensive pipeline, offering flexibility to enhance and complement the analysis of **codon-anticodon dynamics** across various biological contexts.

![tTEscanR_workflow](https://github.com/user-attachments/assets/b74af6d7-3d89-46df-8885-98104b448d36)

## 2. Setup

The tTEscanR source code is entirely written in R. The latest full distribution release can be downloaded from GitHub:

```{r}
install.packages("(avarassanchez/tTEscanR")
library(tTEscanR)
```

### 2.1. Data specifications

The **input** data is expected to be **already pre-processed** following the steps outlined in [Add reference] or other desired user-custom pipelines to obtain gene expression count matrices. Both mRNA and tRNA inputs should be given as gene expression count matrices with features as rows and conditions as columns. In this context, features represent mRNA transcripts or tRNA genes. The conditions will differ based on the data source.  In **bulk** datasets, conditions typically represent a combination of model and replicate, whereas in **single-cell** datasets, conditions correspond to specific tissue and cell type designations. For the features tTEscanR considers mRNA genes and confidently predicted tRNA genes.

Further details can be found in the article referenced.

## 3. Workflow functionalities

### 3.1. Defining the tTEscanR object

The **tTEscanR** object is dynamically updated to store matrices and metadata at each analysis step, ensuring efficient tracking and organization of inputs and outputs throughout the pipeline. It can be created using a single dataset and later updated with additional datasets, or all desired datasets can be in a list and passed directly to the main function. 

**Function:** `Create_tTEscanR_Object()`. 

Creates a **tTEscanR** object that will contain an assays and meta.data slot and individual slots as specified in the parameters. 
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

### 3.2. Codon-anticodon usage assessment

**Codon usage** is computed through matrix multiplication between the mRNA data (genes as rows and conditions as columns) and a reference **codon frequency per gene** table. This reference table consists of **Ensembl** protein-coding genes formatted as a gene matrix, with the 61 sense codons as rows and each  gene as column. As a result, the codon usage calculation produces a matrix where **codons are represented as rows and conditions as columns**. This transformation enables downstream analysis of codon preferences across different conditions.

**Function:** `ComputeCodonUsage()`. 

Updates the **tTEscanR** object by computing the codon usage of the mRNA data using matrix multiplication with a reference codon frequency table. The function allows for customization using a user-provided codon frequency per gene table (`codon_freq`) or default species-specific tables available in **tTEscanR**. 
**Parameters:**
- `object`: The existing **tTEscanR** object containing a mRNA assay.
- `codon_freq`: User-provided codon frequency per gene table.
- `species`: Label specific for a default codon frequency table ("hg38" for human, "mm39" for mouse).
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

<hr>

**Anticodon usage** is computed from a tRNA gene expression matrix including high-confidence tRNA genes predicted using tRNAscan-SE [Addd reference] as rows, and conditions as columns. To create a more concise representation, tRNA genes sharing the same anticodon are grouped together by **tRNA isotypes**, producing a reduced matrix. 

**Function:** `ComputeAnticodonUsage()`.

Updates the **tTEscanR** object by computing the anticodon usage from the tRNA expression data. It groups tRNA genes that share the same anticodon, generating a reduced expression matrix where tRNA isotypes are represented as rows and conditions as columns.
**Parameters:**
- `object`: The existing **tTEscanR** object containing a tRNA assay.
- `overwrite.assay`: Logical, if TRUE, overwrites the codon usage assay present in the object.
- `overwrite.metadata`: Logical, if TRUE, overwrites the metadata present in the object.
- `verbose`: Logical, if TRUE, displays information messages. 

### 3.3. Amino acid level assessment
### 3.4. tTE computation
### 3.5. Helper functions

**Gene annotation**

**Function:** `TranslateGeneName()`

Description: Translates genes from the Ensembl id annotation to the gene name format, and vice-versa.
Parameters: data_to_translate, translator_table, species, position, notation, verbose

## 4. Cheat sheet
## 5. Usage

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

# If we wanted to add at the same time both datasets we can use the following sintaxis:
tTEobject_paired <- Create_tTEscanR_Object(counts = list(mRNA_data, tRNA_data), 
                                           assay = list("mRNA", "tRNA"))

# Adding metadata to the object
matching_celltypes <- intersect(colnames(mRNA_data), colnames(tRNA_data)) # We are defining a vector of strings
tTEobject <- Update_tTEscanObject(object = tTEobject, 
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

## 6. References
