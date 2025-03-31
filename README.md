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

**`Create_tTEscanR_Object()`**. Creates a **tTEscanR** object that will contain an "assays" and "meta.data" slots and in each individual sections as specified in the parameters. 

**`Update_tTEscanR_Object()`**. Updates an existing **tTEscanR** object by overwriting existing sections or adding new ones, as specified in the parameters.

### 3.2. Codon-anticodon usage assessment

The **codon usage** is computed through matrix multiplication between the mRNA gene expression data (genes as rows and conditions as columns) and a **codon frequency per gene** table. This table represents the codon distribution of each protein-coding gene in a reference genome, and includes protein-coding genes from **Ensembl** formatted as a gene matrix, with the 61 sense codons as rows and each gene as columns. As a result, the codon usage calculation produces a matrix where **codons are represented as rows and conditions as columns**. This transformation enables downstream analysis of codon preferences across different conditions.

**`ComputeCodonUsage()`**. Updates the **tTEscanR** object by computing the codon usage of the mRNA data using matrix multiplication with a reference codon frequency table. The function allows for customization using a user-provided codon frequency per gene table (`codon_freq`) or default species-specific tables available in **tTEscanR**. 

<hr>

The **anticodon usage** is computed from a tRNA gene expression matrix including high-confidence tRNA genes predicted using [tRNAscan-SE](https://trna.ucsc.edu/tRNAscan-SE/) as rows, and conditions as columns. Details on how to pre-process the tRNA data described in [Helper functions](#35-helper-functions). To create a more concise representation, tRNA genes sharing the same anticodon are grouped together by **tRNA isotypes**, producing a reduced matrix. 

**`ComputeAnticodonUsage()`**. Updates the **tTEscanR** object by computing the anticodon usage from the tRNA expression data. It groups tRNA genes that share the same anticodon, generating a reduced expression matrix where tRNA isotypes are represented as rows and conditions as columns.

### 3.3. Amino acid level assessment

The **amino acid level** assessment involves calculating amino acid **demand** and **supply**. Demand is determined by grouping codons that correspond to the same amino acid (based on the genetic code), while supply is measured by pooling anticodon families based on the amino acid they recognize. These calculations can be performed separately using "demand" or "supply" or together" by specifying "both".

**`ComputeAAUsage()`**. Computes the amino acid usage, either by assessing the `demand` based on codon usage, or the `supply` based on anticodon usage. You can compute both demand and supply simultaenously, depending on the specified action.
  
### 3.4. tTE computation

The **Theoretical Translation Efficiency (tTE)** is assessed by calculating Spearman's rank correlation coefficient between amino acid demand, derived from mRNA codon usage, and amino acid supply, determined by tRNA anticodon usage. Ensuring **matching conditions** between both datasets is crucial for this step, as only shared conditions are retained for downstream analysis. Consequently, the assessment is highly sensitive to annotation variability. 

**`Compute_tTE()`**. Computes the tTE score across matching conditions at the codon-anticodon usage level and/or amino acid demand and supply level.

### 3.5. Helper functions

**tTEscanR** provides supplementary functions designed to enhance and streamline analysis. These functions are organized into modules and can be used independently or as part of the core functions described above.

-------- **Pre-processing module** --------

The tRNAscanR pipeline requires pre-processed files, as outlined in [Gao et al., 2022](#5-references). However, this module includes several functions designed to assist the user with pre-processing tasks. 
- For **tRNA data**, the pre-processing steps necessary for running tTEscanR involve: (i) create of a tRNA expression matrix, (ii) filter out low tRNA cuts, and (iii) prediction of high-confidence tRNA genes using [tRNAscan-SE](https://trna.ucsc.edu/tRNAscan-SE/).
- For **mRNA data**, the funcitons focus on: (i) trimming protein-coding genes (or any user-defined list of targeted genes), and (ii) ensuring consistent gene annotations across different nomenclature formats.

**`Get_tRNAexpressionMatrix()`**. Extracts tRNA expression data from a Chromatin Assay or a Seurat object and converts it into a structured expression matrix. Optionally, the output can be saved to a specified directory.

**`tRNACutsFilter()`**.Removes conditions where tRNA gene expression falls below a specified threshold. It helps ensure that only conditions with sufficient tRNA counts are included in downstream analyses. 

**`TranslateGeneName()`**. Verifies and, if necessary, translates gene annotations across multiple vectors or data frames, from Ensembl IDs to gene names, or vice-versa. Enables to maintan **consistent annotations** throught the analysis.

<hr>

-------- **Codon frequency table module** --------

To complement the codon usage analysis, we have developed an independent module to generate the **codon frequency per gene matrix** for reference organisms. This matrix can be obtained using Ensembl (through **biomaRt**) to access the reference sequence of an organism or by directly providing a reference genome or sequence files of interest. This flexibility allows users to retrieve codon freqeuncy per gene matrices for any organim available in Ensembl, enabling broad applicability of the computations. However, it is worth mentioning that accessing Ensembl may occasionally result in connection errors due to platform-related issues unrelated to **tTEscanR**. To facilitate a smoother analysis, we have included pre-computed codon frequency per gene matrices for human ("hg38") and mouse ("mm39") as default options. 

**`ObtainCodonComposition()`**. Based on a Ensembl dataset name retrieves the reference genome of the organism and computes the codon frequency per gene matrix. This function can also be used in a targeted mode by using the "transcripts" parameter to give a vector of gene ids to retrieve. 

**`ExtractCodonComposition()`**. Computes the codon composition of a given set of DNA sequences, analyzing the frequency of codons within each sequence.  

<hr>

-------- **Differential expression analysis module** --------

The differential expression analysis module provides a comprehensive ser of functions for identifying and analyzing differentially expressed features across multiple biological levels, including gene expression, codon-anticodon interactions, and amino acid usage. It supports various statistical frameworks, normalization techniques, and visualization tools to ensure accurate and comprehensive analysis.

**`DESeq2runner()`**. Performs differential expression analysis using the DESeq2 package and provides several output visualizations (heatmap, PCA plot and/or volcano plot), as well as normalizes the input data and gives the corrected data object. 

## 4. Usage

```{r}
library(tTEscanR)
data(mRNA_data)
data(tRNA_data)

# Adding the mRNA and tRNA datasets to the object
tTEobject <- Create_tTEscanR_Object(counts = c(mRNA_data, tRNA_data), assay = c("mRNA", "tRNA)) 

# Adding metadata to the object
matching_celltypes <- intersect(colnames(mRNA_data), colnames(tRNA_data)) 
tTEobject <- Update_tTEscanR_Object(object = tTEobject, 
                                    meta.data = list(matching_celltypes), 
                                    meta.data.ids = list("matching_celltypes"), 
                                    overwrite.metadata = TRUE)

# Compute codon usage with default hg38 canonical codon frequency per gene table.
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
                         name_sep = "-")
```

## 5. References

Gao W, Gallardo-Dodd CJ, Kutter C. *Cell type-specific analysis by single-cell profiling identifies a stable mammalian tRNA-mRNA interface and increased translation efficiency in neurons.* Genome Res. 2022;32(1):97-110. [doi:10.1101/gr.275944.121](https://doi.org/10.1101/gr.275944.121)
