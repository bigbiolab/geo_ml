# Load required libraries
library(tidyverse)
library(DESeq2)         # For differential expression analysis
library(GEOquery)       # For fetching GEO datasets
library(org.Hs.eg.db)   # Annotation database for Homo sapiens
library(hgu133plus2.db) # Affymetrix probe annotation for hgu133plus2 arrays
library(dplyr)          # For data manipulation

# Step 1: Download GEO Dataset
# Fetch GEO dataset (example accession: "GSE21942")
geo_data <- getGEO("GSE21942", GSEMatrix = TRUE)
geo_data

# Step 2: Extract Expression Data
# Extract count data (genes x samples)
counts_data <- exprs(geo_data[[1]])

# Extract sample metadata
col_data <- pData(geo_data[[1]])
col_data$condition <- c(rep("Control", 15), rep("MS", 14))
col_data$condition <- as.factor(col_data$condition)

# Verify that sample names match between count data and metadata
all(colnames(counts_data) == rownames(col_data))

# Step 3: DESeq2 Differential Expression Analysis
# Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_data),
  colData = col_data,
  design = ~ condition
)

# Run DESeq2 pipeline
dds <- DESeq(dds)

# Extract results with adjusted p-value threshold (alpha = 0.05)
res <- results(dds, alpha = 0.05)
res_df <- as.data.frame(res)

# Step 4: Gene Annotation
# Check available keytypes
keytypes(hgu133plus2.db)

res_df$symbol <- mapIds(
  hgu133plus2.db,
  keys = rownames(res_df),
  keytype = "PROBEID",
  column = "SYMBOL"
)

res_df_clean <- drop_na(res_df)

# Step 5: Filter Significant Genes
# Filter genes with adjusted p-value < 0.05 and absolute log2FC > 1
sig_genes_df <- res_df %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1)

# Check number of significant genes
nrow(sig_genes_df)

# Step 6: Prepare ML Dataset
# Extract significant probe IDs and their corresponding gene symbols
sig_gene_probes <- rownames(sig_genes_df)
sig_gene_symbols <- sig_genes_df$symbol

# Subset count data to only significant genes
exprs_filtered <- counts_data[sig_gene_probes, , drop = FALSE]

# Replace probe IDs with gene symbols for the dataset column names
rownames(exprs_filtered) <- sig_gene_symbols

exprs_filtered <- exprs_filtered[complete.cases(exprs_filtered), ]


# Transpose expression data: samples as rows, genes as columns
ml_dataset <- as.data.frame(t(exprs_filtered))

# Add class labels (conditions) for ML
ml_dataset$Outcome <- col_data$condition



names(ml_dataset)
