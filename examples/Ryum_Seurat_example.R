# Ryum_Seurat_example.R
# This R script provides a simple example of loading the Ricordea yuma single-cell RNA-seq dataset into Seurat. It demonstrates how to import the gene expression count matrix and cell-type annotations, perform standard preprocessing and dimensionality reduction, and visualise the annotated cell types using UMAP.

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

# Visualize UMAP by annotated cell type
DimPlot(
  ryum,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
)