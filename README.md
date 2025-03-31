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

The **input** data for **tTEscanR** should be pre-processed according to the steps in [Gao et al., 2022](#5-references) or any custom pipeline to obtain gene expression count matrices. Both mRNA and tRNA data must be provided as matrices with features (mRNA transcripts or tRNA genes) as rows and conditions (e.g., model andreplicate in bulk datasets or tissue and cell type in single-cell datasets) as columns. Some functions are integrated into the [pre-processing module](#32-helper-functions) for convenience.

## 3. Workflow functionalities

The **tTEscanR** package provides a strucutred framework for analyzing codon-anticodon pools and transaltion efficiency.

![tTEscanR_workflow](https://github.com/user-attachments/assets/2cf26e1c-5e48-40d6-b81b-300e6293c98a)

### 3.1. Core functions

<table>
  <thead>
    <tr>
      <th>Category</th>
      <th>Function</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td rowspan="2"><b>Defining the tTEscanR object</b></td>
      <td><code>Create_tTEscanR_Object()</code></td>
      <td>Initializes a <b>tTEscanR</b> object to store analysis data.</td>
    </tr>
    <tr>
      <td><code>Update_tTEscanR_Object()</code></td>
      <td>Modifies or extends an existing <b>tTEscanR</b> object.</td>
    </tr>
    <tr>
      <td rowspan="2"><b>Codon-anticodon usage assessment</b></td>
      <td><code>ComputeCodonUsage()</code></td>
      <td>Calculates codon usage by matrix multiplication of mRNA expression data with codon frequency tables.</td>
    </tr>
    <tr>
      <td><code>ComputeAnticodonUsage()</code></td>
      <td>Determines anticodon usage by aggregating tRNA expression data at the anticodon level.</td>
    </tr>
    <tr>
      <td><b>Amino acid level assessment</b></td>
      <td><code>ComputeAAUsage()</code></td>
      <td>Computes amino acid demand (from codon usage) and supply (from anticodon usage), either separately or together.</td>
    </tr>
    <tr>
      <td><b>tTE computation</b></td>
      <td><code>Compute_tTE()</code></td>
      <td>Calculates <b>Theoretical Translation Efficiency (tTE)</b> by correlating amino acid demand and supply across matched conditions.</td>
    </tr>
  </tbody>
</table>

### 3.2. Helper functions

The **tTEscanR** package includes helper functions to support specific steps of the analysis. These functions are organized into modules and can be used independently or as part of the core functions described above.

<table>
  <thead>
    <tr>
      <th>Module</th>
      <th>Function</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td rowspan="3"><b>Pre-processing</b></td>
      <td><code>Get_tRNAexpressionMatrix()</code></td>
      <td>Extracts and structures tRNA expression data from various sources.</td>
    </tr>
    <tr>
      <td><code>tRNACutsFilter()</code></td>
      <td>Filters out conditions with low tRNA expression to ensure data quality.</td>
    </tr>
    <tr>
      <td><code>TranslateGeneName()</code></td>
      <td>Converts gene annotations between Ensembl IDs and gene names for consistency.</td>
    </tr>
    <tr>
      <td rowspan="2"><b>Codon frequency table</b></td>
      <td><code>ObtainCodonComposition()</code></td>
      <td>Retrieves reference genomes and computes codon frequency per gene matrices.</td>
    </tr>
    <tr>
      <td><code>ExtractCodonComposition()</code></td>
      <td>Analyzes codon composition from a given DNA sequence set.</td>
    </tr>
    <tr>
      <td><b>Differential expression analysis</b></td>
      <td><code>DESeq2runner()</code></td>
      <td>Conducts differential expression analysis using DESeq2.</td>
    </tr>
  </tbody>
</table>

## 4. Usage

```{r}
library(tTEscanR)
data(mRNA_data)
data(tRNA_data)

# Adding the mRNA and tRNA datasets to the object
tTEobject <- Create_tTEscanR_Object(counts = c(mRNA_data, tRNA_data),
                                    assay = c("mRNA", "tRNA)) 

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
