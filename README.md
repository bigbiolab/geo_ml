# 🧬 Differential Expression Analysis and ML Dataset Preparation 
This repository provides a **reproducible R-based pipeline** for analyzing **gene expression datasets** from the **NCBI GEO (Gene Expression Omnibus)** database. It performs differential expression analysis using `DESeq2`, annotates genes using Bioconductor annotation packages, and prepares a **machine learning-ready dataset** based on the top differentially expressed genes.

> ✅ This framework is dataset-agnostic and can be adapted to any GEO dataset with expression data (e.g., Affymetrix microarrays).  
> 📌 The included example uses [GSE21942](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21942), a dataset related to Multiple Sclerosis.

---

## 📦 Features

- Download GEO expression datasets programmatically
- Extract and preprocess expression and metadata
- Perform differential expression analysis with `DESeq2`
- Annotate probe IDs using Bioconductor databases (e.g., `hgu133plus2.db`)
- Select top upregulated and downregulated genes
- Export machine learning-ready dataset for classification or feature selection

---
## 🛠️ Dependencies

Ensure the following R packages are installed:

```r
install.packages("tidyverse")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2",
  "GEOquery",
  "org.Hs.eg.db",      # Human gene annotation
  "hgu133plus2.db"     # Platform-specific probe annotation (change as needed)
))
```

---

## 📬 Contact

For questions, suggestions, or contributions, feel free to open an issue or submit a pull request.

---

## 📄 Citation & Attribution

If you use this pipeline in your research, publication, or analysis workflow, please cite or include a reference to this GitHub repository:https://github.com/bigbiolab/geo_ml/



