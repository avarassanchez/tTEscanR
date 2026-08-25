# tTEscanR <a href="https://your-link-here"><img src="https://github.com/user-attachments/assets/8cca530a-f0a5-4284-bf2e-cf030a8193fa" alt="logo_tTEscanR" align="right" width="120"></a>
<img src="https://img.shields.io/badge/Language-R-blue.svg" style="zoom:100%;" />

## 1. What is tTEscanR?

**tTEscanR**, is a powerful and user-friendly R-based package designed to quantify translation efficiency from bulk and single-cell sequencing data. **tTEscanR** offers a comprehensive approach to integrate gene expression data and chromatin accessibility data. The modular design of **tTEscanR** ensures flexibility, allowing users to either run independent components or a complete pipeline based on their research needs. **tTEscanR** has a user-friendly R-based interface that simplifies the analysis of complex data, even for researchers with minimal computational experience. 

Additionally, **tTEscanR** includes an advanced visualization module that generates high-quality plots, aiding in interpretation of results and enhancing the ability to communicate findings effectively. Applicable to both bulk and single-cell sequencing data, **tTEscanR** provides a versatile tool for a wide range of experimental setups and biological contexts, enabling deeper exploration of translation efficiency and its role in cellular processes, disease mechanisms, and therapeutic development. 

**Key features:**
<li>User-friendly R-based interface.</li>
<li>Compatible with bulk and single-cell data.</li>
<li>Modular and customizable pipeline.</li>
<li>High-quality visualizations.</li>

## 2. Setup

The **tTEscanR** source code is entirely written in R. The latest full distribution release can be downloaded from GitHub:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("tTEscanR")
```

```{r}
library(tTEscanR)
```

### 2.1. Data specifications

The **input** data for **tTEscanR** should be preprocessed according to the steps described in [Gao et al., 2022](#5-references) or any custom pipeline to obtain feature expression count matrices. Both mRNA and tRNA data must be provided as matrices with features (mRNA transcripts or tRNA genes) as rows and conditions (e.g., model and replicate in bulk datasets or tissue and cell type in single-cell datasets) as columns. Additionally, a metadata table should be provided for proper integration of the data. Some functions are integrated into the [preprocessing module](#32-helper-functions) for convenience. 

## 3. Workflow functionalities

The **tTEscanR** package provides a structured framework for analyzing codon-anticodon pools in the context of translation efficiency.
<img width="1950" height="2084" alt="schema_tTEscanR_figure1" src="https://github.com/user-attachments/assets/45d40d3e-57be-4f23-b5da-4a6c810d18c1" />

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
      <td rowspan="2"><b>tTEscanR object definition</b></td>
      <td><code>createObject()</code></td>
      <td>Initializes a <b>tTEscanR</b> object to store analysis data.</td>
    </tr>
    <tr>
      <td><code>updateObject()</code></td>
      <td>Modifies or extends an existing <b>tTEscanR</b> object.</td>
    </tr>
    <tr>
      <td rowspan="2"><b>Codon-anticodon usage assessment</b></td>
      <td><code>computeCodonUsage()</code></td>
      <td>Calculates codon usage by matrix multiplication of mRNA expression data with codon frequency tables.</td>
    </tr>
    <tr>
      <td><code>computeAnticodonUsage()</code></td>
      <td>Determines anticodon usage by aggregating tRNA expression data at the anticodon level.</td>
    </tr>
    <tr>
      <td><b>Amino acid level assessment</b></td>
      <td><code>computeAAUsage()</code></td>
      <td>Computes amino acid demand (from codon usage) and supply (from anticodon usage), either separately or together.</td>
    </tr>
    <tr>
      <td><b>tTE computation</b></td>
      <td><code>computeTheoreticalTE()</code></td>
      <td>Calculates <b>Theoretical Translation Efficiency (tTE)</b> by correlating amino acid demand and supply across matched conditions.</td>
    </tr>
    <tr>
      <td><b>Single execution</b></td>
      <td><code>runPipeline()</code></td>
      <td>Uses all the functions listed above to compute the (i) codon-anticodon usage, (ii) amino acid supply-demand ratios, and (iii) theoretical translation efficiency.</td>
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
      <td rowspan="3"><b>General preprocessing</b></td>
      <td><code>mergeMatrices()</code></td>
      <td>Effectively combines matrices.</td>
    </tr>
    <tr>
      <td><code>groupConditions()</code></td>
      <td>Aggregates the individual columns into conditions based on a reference metadata.</td>
    </tr>
    <tr>
      <td><code>featuresToAA()</code></td>
      <td>Translates based on a specified genetic codes codons and/or anticodons to amino acids.</td>
    </tr>
    <tr>
      <td rowspan="4"><b>tRNA preprocessing</b></td>
      <td><code>tRNGetAMatrix()</code></td>
      <td>Extracts and structures tRNA expression data from various sources.</td>
    </tr>
    <tr>
      <td><code>tRNASetGenes()</code></td>
      <td>Translates the tRNA cut coordinates to their corresponding tRNA gene names.</td>
    </tr>
    <tr>
      <td><code>tRNASetCutoff()</code></td>
      <td>Identifies a dataset-specific cutoff to select the low tRNA abundance cutoff.</td>
    </tr>
    <tr>
      <td><code>tRNAFilterCuts()</code></td>
      <td>Filters out conditions with low tRNA abundance to ensure data quality.</td>
    </tr>
    <tr>
      <td rowspan="2"><b>Codon frequency reference</b></td>
      <td><code>getCodonFreq()</code></td>
      <td>Retrieves reference genomes and computes codon frequency per gene matrices.</td>
    </tr>
    <tr>
      <td><code>extractCodons()</code></td>
      <td>Analyzes codon composition from a given DNA sequence set.</td>
    </tr>
    <tr>
      <td rowspan="4"><b>Additional metrics</b></td>
      <td><code>computeMeanUsage()</code></td>
      <td>Computes the feature's usage across conditions.</td>
    </tr>
    <tr>
      <td><code>computeExonicBackground()</code></td>
      <td>Computes the baseline codon/anticodon usage regardless of the expression levels.</td>
    </tr>
    <tr>
      <td><code>getCorrelationBackground()</code></td>
      <td>Correlates the mean usage against the background.</td>
    </tr>
    <tr>
      <td><code>showPoolContribution()</code></td>
      <td>Examines the contribution of the top features to the overall codon usage.</td>
    </tr>
  </tbody>
</table>

## 4. Usage

```{r}
library(tTEscanR)
data(mRNA_data, tRNA_data, metadata)

# Adding the mRNA and tRNA datasets and the metadata to the object
tTEscanR_obj <- createObject(counts = list(mRNA_data, tRNA_data),
                    assay = list("mRNA", "tRNA"),
                    meta.data = list(metadata, "tissue"),
                    meta.data.ids = list("ConditionsLabels", "CorrectionFactor")) 

# Adding extra information to the object
matching_celltypes <- union(colnames(mRNA_data), colnames(tRNA_data)) 
tTEscanR_obj <- updateObject(object = tTEscanR_obj, 
                    meta.data = list(matching_celltypes), 
                    meta.data.ids = list("matching_celltypes"), 
                    overwrite = TRUE)

# Compute codon usage with default hg38 canonical codon frequency per gene table.
tTEscanR_obj <- computeCodonUsage(object = tTEscanR_obj, 
                    codon_freq = NULL, 
                    species = "hg38", 
                    additional_metrics = TRUE)

tTEscanR_obj <- computeAnticodonUsage(object = tTEscanR_obj)

# The both parameter allows to perform at the same time the AA demand and supply assessment
tTEscanR_obj <- computeAAUsage(object = tTEscanR_obj, level = "both")

tTEscanR_obj <- computeTheoreticalTE(object = tTEscanR_obj, level = "both")
```

## 5. AI statement

During the development of this package, generative AI tools (specifically ChatGPT and Gemini) were utilized to optimize specific code implementations and to assist in designing robust test cases for core functions. All AI-assisted code and unit tests were thoroughly reviewed, manually verified, and validated by the authors to ensure technical accuracy and compliance with Bioconductor standards.

## 5. References

Gao W, Gallardo-Dodd CJ, Kutter C. *Cell type-specific analysis by single-cell profiling identifies a stable mammalian tRNA-mRNA interface and increased translation efficiency in neurons.* Genome Res. 2022;32(1):97-110. [doi:10.1101/gr.275944.121](https://doi.org/10.1101/gr.275944.121)
