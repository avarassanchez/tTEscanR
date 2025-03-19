# tTEscanR <p align="right"><img width="120" alt="tTEscanR_logo" src="https://github.com/user-attachments/assets/213de690-12cd-4fa8-862a-cadbd8872bf4"></p> 

## 1. What is tTEscanR?

**tTEscanR** is a versatile and user-friendly **R package** designed to quantify and analyze the relationship between codon usage in mRNA and the availability of corresponding anticodons in tRNA. The package computes a **theoretical translation efficiency (tTE)** score as a proxy of translation efficiency. As a result, both gene expression and chromatin accessibility profile data are required.

The analysis can be carried out across three hierarchical layers of informaiton: gene expression, codon and anticodon pool, and amino acid levle. This multi-layered approach provides a comprehensive view of translation efficiency. The package is optimized for both bulk and single-cell datasets.

With its modular structure, the **tTEscanR** package allows to run specific compoennts independently or as part of a comprehensive pipeline, offering flexibility to enhance and complement the analysis of **codon-anticodon dynamics** across various biological contexts.

## 2. Setup

The tTEscanR source code is entirely written in R. The latest full distribution release can be downloaded from GitHub:

```{r}
install.packages("(avarassanchez/tTEscanR")
library(tTEscanR)
```

**Data specifications:** The **input** data is expected to be **already pre-processed** following the steps outlined in [Add reference] or other desired user-custom pipelines to obtain gene expression count matrices. Both mRNA and tRNA inputs should be given as gene expression count matrices with features as rows and conditions as columns. In this context, features represent mRNA transcripts or tRNA genes. The conditions will differ based on the data source.  In **bulk** datasets, conditions typically represent a combination of model and replicate, whereas in **single-cell** datasets, conditions correspond to specific tissue and cell type designations. For the features tTEscanR considers mRNA genes and confidently predicted tRNA genes.

Further details can be found in the article referenced.

## 3. List of functionalities

### 3.1. Difinging the tTEscanR object

| Functions | Description | Parameters                                                                                                      | 
| --------- | --------- | --------------------------------------------------------------------------------------------------------------- |
| `Create_tTEscanObject()` | Creates a **tTEscanR** object that will contain an assays and meta.data (optional) slot and individual slots as specified in the parameters. | `counts`: Count matrix (or a list of matrices). <br> `assay`: Label (or list of labels) to identify the `counts`. <br> `meta.data`: (Optional) Additional data (single variable or list of variable) to include in the object. <br> `meta.data.ids`: (Optional) Labels to identify the metadata. <br> `verbose`: (Optional) Controls whether messages regarding function execution are displayed. |
| `Update_tTEscanObject()` | Updates an existing **tTEscanR** object by overwriting existing sections or adding new ones, as specified in the parameters. | `object`: The existing tTEscanR object to be updated. <br> `counts`: Count matrix (or a list of matrices). <br> `assay`: Label (or list of labels) to identify the `counts`. <br> `meta.data`: (Optional) Additional data (single variable or list of variable) to include in the object. <br> `meta.data.ids`: (Optional) Labels to identify the metadata. <br> `overwrite.assay`: Logical parameter indicating whether to overwrite the existing assays if the `assay` label already exists in the **tTEscanR** object. <br> `overwrite.metadata`: Logical parameter indicating whether to overwrite existing metadata if the `meta.data.ids` label already exists in the **tTEscanR** object. <br> `verbose`: (Optional) Logical value to controls whether messages regarding the function's execution are displayed. | 


