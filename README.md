# tTEscanR <a href="https://your-link-here"><img src="https://github.com/user-attachments/assets/8cca530a-f0a5-4284-bf2e-cf030a8193fa" alt="logo_tTEscanR" align="right" width="120"></a>
<img src="https://img.shields.io/badge/Language-R-blue.svg" style="zoom:100%;" />

## 1. What is tTEscanR?

**tTEscanR**, is a powerful and user-friendly R-based package designed to quantify translation efficiency from bulk and single-cell sequencing data. **tTEscanR** offers a comprehensive approach to quantify translation efficiency by integrating gene expression data and chromatin accessibility data. The modular design of **tTEscanR** ensures flexibility, allowing users to either run independent components or a complete pipeline based on their research needs. Its user-friendly R-based interface simplifies the analysis of complex data, even for researchers with minimal computational experience. 

Additionally, **tTEscanR** includes an advanced visualization module that generates high-quality plots, aiding in interpretation of results and enhancing the ability to communicate findings effectively. Applicable to both bulk and single-cell sequencing data, **tTEscanR** provides a versatile tool for a wide range of experimental setups and biological contexts, enabling deeper exploration of translation efficiency and its role in cellular processes, disease mechanisms, and therapeutic development. 

**Key features:**
<li>Compatible with bulk and single-cell data.</li>
<li>Modular and customizable pipeline.</li>
<li>High-quality visualizations.</li>
<li>User-friendly R-based interface.</li>

## 2. Setup

The **tTEscanR** source code is entirely written in R. The latest full distribution release can be downloaded from GitHub:

```{r}
install.packages("/avarassanchez/tTEscanR")
library(tTEscanR)
```

### 2.1. Data specifications

The **input** data for **tTEscanR** should be preprocessed according to the steps described in [Gao et al., 2022](#5-references) or any custom pipeline to obtain feature expression count matrices. Both mRNA and tRNA data must be provided as matrices with features (mRNA transcripts or tRNA genes) as rows and conditions (e.g., model and replicate in bulk datasets or tissue and cell type in single-cell datasets) as columns. Additionally, a metadata table should be provided for proper integration of the data. Some functions are integrated into the [preprocessing module](#32-helper-functions) for convenience. 

## 3. Workflow functionalities

The **tTEscanR** package provides a strucutred framework for analyzing codon-anticodon pools and transaltion efficiency.
<img width="2024" height="1312" alt="tTEscanR_fig1_schema_GH" src="https://github.com/user-attachments/assets/9efc89dc-6acd-41ec-929d-b0cc76445067" />

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
    <tr>
      <td><b>Single execution</b></td>
      <td><code>Run_tTEscanR_pipeline()</code></td>
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
      <td><code>MergeMatrices()</code></td>
      <td>Effectively combines matrices.</td>
    </tr>
    <tr>
      <td><code>MergeConditions()</code></td>
      <td>Aggregates the individual columns into conditions based on a reference metadata.</td>
    </tr>
    <tr>
      <td><code>TransformFormat()</code></td>
      <td>Modifies the data to fit the requirements of the visualization module.</td>
    </tr>
    <tr>
      <td rowspan="4"><b>tRNA preprocessing</b></td>
      <td><code>Get_tRNAMatrix()</code></td>
      <td>Extracts and structures tRNA expression data from various sources.</td>
    </tr>
    <tr>
      <td><code>Set_tRNAgenes()</code></td>
      <td>Translates the tRNA cut coordinates to their corresponding tRNA gene names.</td>
    </tr>
    <tr>
      <td><code>Set_tRNACutoff()</code></td>
      <td>Identifies a dataset-specific cutoff to select the low tRNA abundance cutoff.</td>
    </tr>
    <tr>
      <td><code>Filter_tRNACuts()</code></td>
      <td>Filters out conditions with low tRNA abundance to ensure data quality.</td>
    </tr>
    <tr>
      <td rowspan="2"><b>Codon frequency reference</b></td>
      <td><code>GetCodonFreq()</code></td>
      <td>Retrieves reference genomes and computes codon frequency per gene matrices.</td>
    </tr>
    <tr>
      <td><code>ExtractCodons()</code></td>
      <td>Analyzes codon composition from a given DNA sequence set.</td>
    </tr>
    <tr>
      <td rowspan="4"><b>Additional metrics</b></td>
      <td><code>ComputeMeanUsage()</code></td>
      <td>Computes the feature's usage across conditions.</td>
    </tr>
    <tr>
      <td><code>ComputeExonicBackground()</code></td>
      <td>Computes the baseline codon/anticodon usage regardless of the expression levels.</td>
    </tr>
    <tr>
      <td><code>GetCorrelationBackground()</code></td>
      <td>Correlates the mean usage against the background.</td>
    </tr>
    <tr>
      <td><code>ShowPoolContribution()</code></td>
      <td>Examines the contribution of the top features to the overall codon usage.</td>
    </tr>
  </tbody>
</table>

## 4. Usage

```{r}
library(tTEscanR)
data(mRNA_data, tRNA_data, metadata)

# Adding the mRNA and tRNA datasets and the metadata to the object
tTEscanR_obj <- Create_tTEscanR_Object(counts = list(mRNA_data, tRNA_data),
                                       assay = list("mRNA", "tRNA"),
                                       meta.data = list(metadata, "tissue"),
                                       meta.data.ids = list("ConditionsLabels", "CorrectionFactor")) 

# Adding extra information to the object
matching_celltypes <- intersect(colnames(mRNA_data), colnames(tRNA_data)) 
tTEscanR_obj <- Update_tTEscanR_Object(object = tTEscanR_obj, 
                                       meta.data = list(matching_celltypes), 
                                       meta.data.ids = list("matching_celltypes"), 
                                       overwrite.metadata = TRUE)

# Compute codon usage with default hg38 canonical codon frequency per gene table.
tTEscanR_obj <- ComputeCodonUsage(object = tTEscanR_obj, 
                                  codon_freq = NULL, 
                                  species = "hg38", 
                                  filter = "canonical",
                                  additional.metrics = TRUE)

tTEscanR_obj <- ComputeAnticodonUsage(object = tTEscanR_obj)

# The both parameter allows to perform at the same time the AA demand and supply assessment
tTEscanR_obj <- ComputeAAUsage(object = tTEscanR_obj, level = "both")

tTEscanR_obj <- Compute_tTE(object = tTEscanR_obj, level = "both")
```

## 5. References

Gao W, Gallardo-Dodd CJ, Kutter C. *Cell type-specific analysis by single-cell profiling identifies a stable mammalian tRNA-mRNA interface and increased translation efficiency in neurons.* Genome Res. 2022;32(1):97-110. [doi:10.1101/gr.275944.121](https://doi.org/10.1101/gr.275944.121)
