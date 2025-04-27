# Load required libraries
library(tidyverse)
library(DESeq2)         # For differential expression analysis
library(GEOquery)       # For fetching GEO datasets
library(org.Hs.eg.db)   # Annotation database for Homo sapiens
library(hgu133plus2.db) # Affymetrix probe annotation for hgu133plus2 arrays
library(tidyverse)      # For data manipulation

# Step 1: Download GEO Dataset
# Fetch GEO dataset (example accession: "GSE21942")
geo_data <- getGEO("GSE21942", GSEMatrix = TRUE, destdir = "./data")
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

# Pre-filtering: remove rows with low gene counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[rowSums(counts(dds)) >= 10]

# Set the factor level
dds$condition <- relevel(dds$condition, ref = "Control")

# Run DESeq2 pipeline
dds <- DESeq(dds)

# Extract results
res <- results(dds)
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

# Step 5: Filter Significant Genes
# Filter genes with adjusted p-value < 0.05 and absolute log2FC > 1
top25_downregulated <- res_df  |>
  dplyr::filter(!is.na(padj), !is.na(symbol), padj < 0.05)  |>
  dplyr::arrange(log2FoldChange)  |>
  dplyr::slice(1:25)

top25_upregulated <- res_df  |>
  dplyr::filter(!is.na(padj), !is.na(symbol), padj < 0.05)  |>
  dplyr::arrange(desc(log2FoldChange))  |>
  dplyr::slice(1:25)

# Combined datasets
combined_df <- bind_rows(top25_upregulated, top25_downregulated)

# Check number of significant genes
nrow(combined_df)

# Step 6: Prepare ML Dataset
# Extract  probe IDs and their corresponding gene symbols
probes_id <- rownames(combined_df)
gene_symbols <- combined_df$symbol

# Subset count data
exprs_filtered <- counts_data[probes_id, , drop = FALSE]


# Replace probe IDs with gene symbols for the dataset column names
rownames(exprs_filtered) <- gene_symbols

# Transpose expression data: samples as rows, genes as columns
ml_dataset <- as.data.frame(t(exprs_filtered))

# Add class labels (conditions) for ML
ml_dataset$Outcome <- col_data$condition
names(ml_dataset)

write.csv(ml_dataset, "data/MS_GSE21942.csv", row.names = FALSE)
