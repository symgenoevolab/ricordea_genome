# Ryum_visualization.R
# This R script loads the Ricordea yuma single-cell RNA-seq dataset into Seurat, integrates cell-type annotations and BLASTP best-hit annotations against the human UniProt proteome, and performs standard preprocessing and UMAP dimensionality reduction. Gene identifiers are supplemented with human gene symbols in lowercase to facilitate gene identification. A custom function, plot_gene(), allows users to visualize the expression of genes of interest by searching for their human UniProt best-hit symbols.

library(Seurat)

# Load cell-type annotations
meta <- read.delim(
  "Ryum.cell_to_cts.csv.gz",
  row.names = 1
)

# Load gene expression count matrix
counts <- ReadMtx(
  mtx = "Ryum.counts.mtx.gz",
  cells = "Ryum.cell_to_cts.csv.gz",
  cell.column = 1,
  skip.cell = 1,
  features = "Ryum.genes.txt.gz",
  feature.column = 1
)

# Load BLASTP annotations against the human UniProt proteome
blast <- read.delim(
  "Ryum_to_Hsap_UniProt_blastp.tsv",
  header = FALSE
)

# Extract human gene symbols and convert to lowercase
blast$symbol <- sub(
  ".* GN=([^ ]+).*",
  "\\1",
  blast$V7
)

blast$symbol[!grepl(" GN=", blast$V7, fixed = TRUE)] <- NA
blast$symbol <- tolower(blast$symbol)

# Retain the first BLASTP hit for each gene
blast <- blast[!duplicated(blast$V1), ]

# Create gene ID to human UniProt symbol mapping
gene_map <- setNames(blast$symbol, blast$V1)

# Append human gene symbols to original gene IDs
genes <- rownames(counts)
symbols <- gene_map[genes]

new_genes <- ifelse(
  !is.na(symbols),
  paste0(genes, "-", symbols),
  genes
)

rownames(counts) <- make.unique(new_genes)

# Align metadata with cell barcodes
stopifnot(all(colnames(counts) %in% rownames(meta)))
meta <- meta[colnames(counts), , drop = FALSE]

# Create Seurat object
ryum <- CreateSeuratObject(
  counts = counts,
  meta.data = meta,
  project = "Ryum"
)

# Inspect cell-type distribution
table(ryum$cell_type)

# Normalize data and identify variable features
ryum <- NormalizeData(ryum)
ryum <- FindVariableFeatures(ryum)

# Scale data and perform PCA
ryum <- ScaleData(ryum)
ryum <- RunPCA(ryum, npcs = 50)

# Construct nearest-neighbor graph and compute UMAP
ryum <- FindNeighbors(ryum, dims = 1:30)
ryum <- RunUMAP(ryum, dims = 1:30)

# Define function to plot gene expression by gene symbol
plot_gene <- function(symbol) {
  genes <- grep(
    paste0("-", symbol, "$"),
    rownames(ryum),
    value = TRUE,
    ignore.case = TRUE
  )
  
  if (length(genes) == 0) {
    stop("Gene not found: ", symbol)
  }
  
  FeaturePlot(ryum, features = genes)
}

# Visualize UMAP by annotated cell type
DimPlot(
  ryum,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
)

# Visualize expression of an annotated gene
FeaturePlot(
  ryum,
  features = "g1348-bmp2"
)

# Visualize expression of targeted orthologous genes
plot_gene("bmp2")
plot_gene("wnt7b")

# Save Seurat object
saveRDS(ryum, "Ryum_seurat.rds")