# tTEscanR

### 1. What is tTEscanR?
<img width="365" alt="tTEscanR_logo" src="https://github.com/user-attachments/assets/213de690-12cd-4fa8-862a-cadbd8872bf4" />


**tTEscanR** is a versatile and user-friendly **R package** designed to quantify and analyze the relationship between codon usage in mRNA and the availability of corresponding anticodons in tRNA. The package computes a **theoretical translation efficiency (tTE)** score as a proxy of translation efficiency. As a result, both gene expression and chromatin accessibility profile data are required.

The analysis can be carried out across three hierarchical layers of informaiton: gene expression, codon and anticodon pool, and amino acid levle. This multi-layered approach provides a comprehensive view of translation efficiency.

With its modular structure, the **tTEscanR** package allows to run specific compoennts independently or as part of a comprehensive pipeline, offering flexibility to enhance and complement the analysis of **codon-anticodon dynamics** across various biological contexts.

### 2. Setup

The tTEscanR source code is entirely written in R. The latest full distribution release can be downloaded from GitHub:

**Data specifications:** The accepted mRNA and tRNA inputs for tTEscanR consist of pre-processed gene expression count matrices, with features as rows and conditions as columns. The package is optimized for both bulk and single-cell datasets. In **bulk** datasets, conditions typically represent a combination of model and replicate, whereas in **single-cell** datasets, conditions correspond to specific tissue and cell type designations. For the features tTEscanR considers mRNA genes and confidently predicted tRNA genes.

Further details can be found in the article referenced.
