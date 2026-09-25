###############################################################################
# 02_seurat_analysis.R
# Seurat v5 processing for G. fascicularis and R. yuma + loading of public species
#
# Author : Victor Pinon-Gonzalez
# Date   : May 2026
#
# Input:  CR_mapping_output/Gast/   (CellRanger output)
#         CR_mapping_output/Ryum/   (CellRanger output)
# Output: seurat_gfas_ryum/Gast/   (annotated RDS, UMAP, metadata, markers)
#         seurat_gfas_ryum/Ryum/
#
# Run order:
#   Rscript code/02_Amur_seurat_obj_reconstruction.R   # produces Amur RDS
#   Rscript code/03_seurat_analysis.R
#
# Overview
# --------
# 1. SoupX ambient RNA correction (Gast and Ryum)
# 2. QC filtering, SCTransform v2, PCA / UMAP / clustering
# 3. Manual cell type assignment (guided by step 1 heatmap in 04_calicoblast_cells_analysis.R)
# 4. UMAP plots
# 5. [OPTIONAL — uncomment to run] Expressed gene extraction for Gast CellRanger re-mapping
# 6. Save annotated Seurat objects
# 7. Per-cell metadata tables
# 8. Top 50 marker genes per cell type (FindAllMarkers, Wilcoxon)
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(SoupX)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# ---------------------------------------------------------------------------
# Project root: works with Rscript, RStudio Source button, or RStudio console.
# ---------------------------------------------------------------------------
BASE_DIR <- local({
  f <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) return(dirname(dirname(normalizePath(sub("--file=", "", f)))))
  p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
  if (isTRUE(nzchar(p))) return(dirname(dirname(normalizePath(p))))
  wd <- normalizePath(getwd())
  if (basename(wd) == "code") dirname(wd) else wd
})
cat("Project root:", BASE_DIR, "\n")

# CellRanger output directories (symlinks under CR_mapping_output/)
GFAS_RAW  <- file.path(BASE_DIR, "CR_mapping_output", "Gfas", "raw_feature_bc_matrix")
GFAS_FILT <- file.path(BASE_DIR, "CR_mapping_output", "Gfas", "filtered_feature_bc_matrix")
RYUM_RAW  <- file.path(BASE_DIR, "CR_mapping_output", "Ryum", "raw_feature_bc_matrix")
RYUM_FILT <- file.path(BASE_DIR, "CR_mapping_output", "Ryum", "filtered_feature_bc_matrix")

# Output
OUT_DIR <- file.path(BASE_DIR, "seurat_gfas_ryum")
dir.create(file.path(OUT_DIR, "Gfas"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "Ryum"), showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Cell type colour palette (shared across all plots)
# ---------------------------------------------------------------------------
cell_type_colors <- c(
  "Gastrodermis"        = "#D99A7C",
  "Gastrodermis_1"      = "#D99A7C",
  "Gastrodermis_2"      = "#E3AA8F",
  "Epidermis"           = "#6F8FBF",
  "Epidermis_1"         = "#6F8FBF",
  "Epidermis_2"         = "#89A6CF",
  "Calicoblast"         = "#9A88B5",
  "Neurons"             = "#7FB69A",
  "Neurons_1"           = "#7FB69A",
  "Neurons_2"           = "#97C6AE",
  "Neurons-like"        = "#CAE8D8",
  "Cnidocytes"          = "#D6B26E",
  "Algae_hosting_cells" = "#E4CF85",
  "Gland cells-like"    = "#B79C7A",
  "Immune cells"        = "#C48CB3",
  "Unknown"             = "#9E9E9E",
  "unknown"             = "#9E9E9E"
)

# =============================================================================
# Helper — SoupX ambient RNA correction
# CellRanger output uses uncompressed TSV files (barcodes.tsv / features.tsv /
# matrix.mtx), so Read10X() is not used.
# SoupX requires preliminary cluster labels, so a first-pass SCTransform →
# PCA → clustering is run on the filtered matrix before building SoupChannel.
# Dimensions and resolution must match the original analysis.
# =============================================================================
run_soupx <- function(raw_dir, filt_dir, species_label,
                      prelim_dims, prelim_resolution) {
  cat("  Loading matrices for", species_label, "...\n")

  tod <- readMM(file.path(raw_dir, "matrix.mtx"))
  rownames(tod) <- read.delim(file.path(raw_dir, "features.tsv"),
                               header = FALSE)$V1
  colnames(tod) <- read.delim(file.path(raw_dir, "barcodes.tsv"),
                               header = FALSE)$V1

  toc <- readMM(file.path(filt_dir, "matrix.mtx"))
  rownames(toc) <- read.delim(file.path(filt_dir, "features.tsv"),
                               header = FALSE)$V1
  colnames(toc) <- read.delim(file.path(filt_dir, "barcodes.tsv"),
                               header = FALSE)$V1

  cat("  Preliminary clustering for SoupX (", species_label, ")...\n")
  pre_obj <- CreateSeuratObject(counts = toc, min.cells = 0, min.features = 0)
  pre_obj <- SCTransform(pre_obj, vst.flavor = "v2", verbose = FALSE)
  pre_obj <- RunPCA(pre_obj, npcs = 50, verbose = FALSE)
  pre_obj <- RunUMAP(pre_obj, dims = 1:prelim_dims, verbose = FALSE)
  pre_obj <- FindNeighbors(pre_obj, dims = 1:prelim_dims, verbose = FALSE)
  pre_obj <- FindClusters(pre_obj, resolution = prelim_resolution, verbose = FALSE)
  cluster_labels <- Idents(pre_obj)

  cat("  Running SoupX for", species_label, "...\n")
  common_barcodes <- intersect(colnames(toc), colnames(pre_obj))
  cat("    Common barcodes:", length(common_barcodes), "\n")
  toc_sub <- toc[, common_barcodes]
  sc <- SoupChannel(tod = tod, toc = toc_sub)
  sc <- setClusters(sc, clusters = cluster_labels)
  sc <- autoEstCont(sc)
  adjustCounts(sc)
}

# =============================================================================
# Helper — SCTransform + PCA + clustering pipeline (post-SoupX)
#   CreateSeuratObject(min.cells=3, min.features=100) → PercentageFeatureSet →
#   VlnPlot → subset (strict > / <) → SCTransform(glmGamPoi, regress percent.mt) →
#   RunPCA → ElbowPlot → FindNeighbors(n_dims) → FindClusters → RunUMAP
# Reproducibility relies on Seurat's built-in seed.use defaults (no global seed).
# =============================================================================
run_seurat_pipeline <- function(counts_mat, species, mt_pattern,
                                n_dims, resolution,
                                nFeature_min, nFeature_max,
                                nCount_min,   nCount_max,
                                pct_mt_max,   out_dir) {
  cat("  Creating Seurat object for", species, "...\n")
  obj <- CreateSeuratObject(counts = counts_mat, project = species,
                            min.cells = 3, min.features = 100)
  obj[["species"]] <- species

  if (length(mt_pattern) > 1) {
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, features = mt_pattern)
  } else {
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
  }

  cat(sprintf("  Pre-filter: %d cells\n", ncol(obj)))

  p_vln <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                   pt.size = 0) &
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
  ggsave(file.path(out_dir, "QC_violin_before_filter.pdf"), p_vln,
         width = 10, height = 4)

  obj <- subset(obj,
                nFeature_RNA > nFeature_min & nFeature_RNA < nFeature_max &
                nCount_RNA   > nCount_min   & nCount_RNA   < nCount_max   &
                percent.mt   < pct_mt_max)
  cat(sprintf("  Post-filter: %d cells\n", ncol(obj)))

  obj <- SCTransform(obj, method = "glmGamPoi", vars.to.regress = "percent.mt",
                     verbose = FALSE)
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)

  p_elbow <- ElbowPlot(obj, ndims = 50)
  ggsave(file.path(out_dir, "PCA_elbow_plot.pdf"), p_elbow, width = 8, height = 5)

  # FindClusters before RunUMAP — matches original order
  obj <- FindNeighbors(obj, dims = 1:n_dims,     verbose = FALSE)
  obj <- FindClusters(obj,  resolution = resolution, verbose = FALSE)
  obj <- RunUMAP(obj,       dims = 1:n_dims,     verbose = FALSE)

  cat(sprintf("  Clusters at res %.2f: %d\n", resolution,
              length(levels(Idents(obj)))))
  obj
}

# =============================================================================
# ID conversion functions — OrthoFinder protein ID → Seurat gene ID
# =============================================================================
sc_gfas <- function(x) paste0("g.", gsub("_", "-", x))
sc_amur <- function(x) {
  x <- sub("^evm\\.model\\.", "", x)
  x <- sub("^scaffold_", "scaffold-", x)
  paste0("EVM%20prediction%20", x)
}
sc_spis <- function(x) gsub("_", "-", x)
sc_amil <- function(x) { x <- sub("-R[A-Z]$", "", x); gsub("_", "-", x) }
sc_opat <- function(x) { x <- sub("\\.t[0-9]+$", "", x); gsub("_", "-", x) }
sc_ryum <- function(x) { x <- sub("GeneExt~", "", x); sub("\\.t[0-9]+.*$", "", x) }


###############################################################################
# 1. G. fascicularis
###############################################################################
cat("\n=== G. fascicularis processing ===\n")

gfas_corrected <- run_soupx(GFAS_RAW, GFAS_FILT, "G. fascicularis",
                            prelim_dims = 30, prelim_resolution = 0.1)

# Mitochondrial genes identified from Trinity assembly (locus TRINITY-GG-119)
GFAS_MT_PATTERN <- "TRINITY-GG-119-c0-g1-i(12|14|16)"

gfas_seurat <- run_seurat_pipeline(
  counts_mat   = gfas_corrected,
  species      = "G. fascicularis",
  mt_pattern   = GFAS_MT_PATTERN,
  n_dims       = 30,   resolution   = 0.1,
  nFeature_min = 192,  nFeature_max = 1612,
  nCount_min   = 657,  nCount_max   = 5342,
  pct_mt_max   = 10,
  out_dir      = file.path(OUT_DIR, "Gfas")
)

###############################################################################
# 2. R. yuma
###############################################################################
cat("\n=== R. yuma processing ===\n")

ryum_corrected <- run_soupx(RYUM_RAW, RYUM_FILT, "R. yuma",
                            prelim_dims = 30, prelim_resolution = 0.1)

# Mitochondrial genes from MitoFinder annotation (exact names — includes ND4L)
RYUM_MT_GENES <- c("COX1", "COX2", "COX3",
                   "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
                   "ATP6", "ATP8", "CYTB")

ryum_seurat <- run_seurat_pipeline(
  counts_mat   = ryum_corrected,
  species      = "R. yuma",
  mt_pattern   = RYUM_MT_GENES,
  n_dims       = 25,   resolution   = 0.1,
  nFeature_min = 223,  nFeature_max = 1658,
  nCount_min   = 326,  nCount_max   = 2347,
  pct_mt_max   = 0.51,
  out_dir      = file.path(OUT_DIR, "Ryum")
)

###############################################################################
# 3. Manual cell type assignment
# Cluster labels were determined by cross-species Spearman correlation
# (04_calicoblast_cells_analysis.R step 1 heatmap) + canonical marker inspection.
###############################################################################
cat("\n=== Note: cluster annotation guided by step 1 heatmap (see 04_calicoblast_cells_analysis.R) ===\n")

gfas_cell_types <- c(
  "0"  = "Gastrodermis",
  "1"  = "Epidermis_1",
  "2"  = "Neurons",
  "3"  = "Calicoblast",
  "4"  = "Epidermis_2",
  "5"  = "Algae_hosting_cells",
  "6"  = "Gland cells-like",
  "7"  = "Immune cells",
  "8"  = "Cnidocytes",
  "9"  = "unknown"
)

ryum_cell_types <- c(
  "0"  = "Algae_hosting_cells",
  "1"  = "Gastrodermis_1",
  "2"  = "Epidermis",
  "3"  = "Neurons_2",
  "4"  = "Neurons-like",
  "5"  = "Gastrodermis_2",
  "6"  = "Unknown",
  "7"  = "Immune cells",
  "8"  = "Neurons_1",
  "9"  = "Cnidocytes"
)

Idents(gfas_seurat) <- "SCT_snn_res.0.1"
gfas_seurat <- RenameIdents(gfas_seurat, gfas_cell_types)
gfas_seurat[["cell_type"]] <- Idents(gfas_seurat)

Idents(ryum_seurat) <- "SCT_snn_res.0.1"
ryum_seurat <- RenameIdents(ryum_seurat, ryum_cell_types)
ryum_seurat[["cell_type"]] <- Idents(ryum_seurat)

###############################################################################
# 4. UMAP plots
###############################################################################
plot_umap <- function(obj, title, out_path) {
  p <- DimPlot(obj, reduction = "umap", group.by = "cell_type",
               label = TRUE, repel = TRUE, pt.size = 0.4) +
    scale_color_manual(values = cell_type_colors) +
    ggtitle(title) +
    theme_classic(base_size = 12) +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "italic"))
  ggsave(out_path, p, width = 9, height = 7)
  invisible(p)
}

cat("\n=== Saving UMAP plots ===\n")
plot_umap(gfas_seurat, "G. fascicularis", file.path(OUT_DIR, "Gfas", "UMAP_GFAS_cell_types.pdf"))
plot_umap(ryum_seurat, "R. yuma",     file.path(OUT_DIR, "Ryum", "UMAP_RYUM_cell_types.pdf"))

###############################################################################
# 5. [OPTIONAL] G. fascicularis expressed gene extraction for CellRanger re-mapping
#
# Trinity produces a large redundant transcriptome (~90,000 transcripts).
# This step extracts transcripts detected in ≥1 cell so the reference can be
# filtered to ~30,000 sequences and CellRanger re-run for better alignment.
# Remove the comment marks (#) below to run this step.
###############################################################################
# cat("\n=== Extracting expressed Gast transcript IDs ===\n")
# gfas_counts    <- GetAssayData(gfas_seurat, assay = "RNA", layer = "counts")
# gfas_expressed <- rownames(gfas_counts)[rowSums(gfas_counts) > 0]
# writeLines(gfas_expressed,
#            file.path(OUT_DIR, "Gfas", "gfas_expressed_transcript_ids.txt"))
# cat("  Expressed transcripts:", length(gfas_expressed), "\n")

###############################################################################
# 6. Save annotated Seurat objects
###############################################################################
cat("\n=== Saving annotated Seurat objects ===\n")
saveRDS(gfas_seurat, file.path(OUT_DIR, "Gfas", "gfas_seurat_annotated.rds"))
saveRDS(ryum_seurat, file.path(OUT_DIR, "Ryum", "ryuma_seurat_annotated.rds"))

cat("\nCell type breakdown — G. fascicularis:\n")
print(sort(table(gfas_seurat$cell_type), decreasing = TRUE))
cat("\nCell type breakdown — R. yuma:\n")
print(sort(table(ryum_seurat$cell_type), decreasing = TRUE))

###############################################################################
# 7. Per-cell metadata tables
###############################################################################
cat("\n=== Saving per-cell metadata tables ===\n")

export_metadata <- function(obj, out_path) {
  meta <- obj@meta.data
  umap <- Embeddings(obj, "umap")
  out  <- cbind(
    data.frame(barcode = rownames(meta), stringsAsFactors = FALSE),
    meta[, intersect(c("cell_type", "seurat_clusters", "nCount_RNA",
                       "nFeature_RNA", "percent.mt", "species"),
                     names(meta)), drop = FALSE],
    UMAP_1 = umap[, 1],
    UMAP_2 = umap[, 2]
  )
  write.table(out, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  Saved:", basename(out_path), "|", nrow(out), "cells\n")
}

export_metadata(gfas_seurat, file.path(OUT_DIR, "Gfas", "gfas_cell_metadata.tsv"))
export_metadata(ryum_seurat, file.path(OUT_DIR, "Ryum", "ryuma_cell_metadata.tsv"))

###############################################################################
# 8. Top 50 marker genes per cell type (Wilcoxon, FindAllMarkers)
###############################################################################
cat("\n=== Finding top markers per cell type ===\n")
options(future.globals.maxSize = Inf)

find_top_markers <- function(obj, out_path, top_n = 50) {
  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.10,
                            logfc.threshold = 0.25, test.use = "wilcox",
                            verbose = FALSE)
  if (nrow(markers) == 0) { warning("No markers found"); return(invisible(NULL)) }
  out <- markers %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    select(sc_gene_id      = gene,
           cell_type        = cluster,
           avg_log2FC,
           pct_expressing   = pct.1,
           pct_other_cells  = pct.2,
           p_val_adj) %>%
    mutate(across(where(is.numeric), ~ round(., 4))) %>%
    arrange(cell_type, desc(avg_log2FC))
  write.table(out, out_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  cat("  Saved:", basename(out_path), "|", nrow(out), "rows\n")
}

find_top_markers(gfas_seurat, file.path(OUT_DIR, "Gfas", "Gfas_markers_per_celltype.tsv"))
find_top_markers(ryum_seurat, file.path(OUT_DIR, "Ryum", "Ryum_markers_per_celltype.tsv"))

###############################################################################
# 9. Overlapped UMAPs — G. fascicularis + R. yuma
#
#  Two versions produced:
#
#  A) SAMap joint UMAP  (file: UMAP_Gfas_Ryum_samap.pdf/png)
#     Joint embedding from SAMap pairwise (gfas × ryum). Cross-species homology
#     (DIAMOND) bridges the gene namespaces, so equivalent cell types cluster
#     together across species. Run SAMap/pairwise_comps/gfas_ryum/run_samap_gfas_ryum.py
#     first to generate umap_coords.csv.
#
#  B) Seurat merged UMAP  (file: UMAP_Gfas_Ryum_seurat.pdf/png)
#     Simple merge of both objects on common genes + re-normalisation. Species
#     tend to occupy distinct regions (gene namespace differences), showing the
#     overall species-level structure.
#
#  Both versions use the same 3-color scheme:
#    G. fascicularis cells        → muted gray-red  (#DCBAB6)
#    R. yuma cells            → muted gray-blue (#B6CADC)
#    G. fascicularis Calicoblast  → #9A88B5 (highlighted, on top)
###############################################################################

overlay_colors <- c(
  "G. fascicularis"               = "#DCBAB6",
  "R. yuma"                   = "#B6CADC",
  "Calicoblast (G. fascicularis)" = cell_type_colors[["Calicoblast"]]
)

add_plot_group <- function(df, species_col, celltype_col) {
  df$plot_group <- ifelse(
    df[[species_col]] %in% c("gfas", "G. fascicularis") &
      df[[celltype_col]] == "Calicoblast",
    "Calicoblast (G. fascicularis)",
    ifelse(df[[species_col]] %in% c("gfas", "G. fascicularis"),
           "G. fascicularis", "R. yuma")
  )
  df$plot_group <- factor(df$plot_group,
    levels = c("R. yuma", "G. fascicularis", "Calicoblast (G. fascicularis)"))
  df[order(df$plot_group), ]
}

make_overlay_plot <- function(df, title) {
  ggplot(df, aes(x = UMAP1, y = UMAP2, color = plot_group)) +
    geom_point(size = 0.3, alpha = 0.7) +
    scale_color_manual(values = overlay_colors,
                       guide  = guide_legend(override.aes = list(size = 3))) +
    ggtitle(title) +
    theme_classic(base_size = 12) +
    theme(legend.position = "right",
          axis.title      = element_text(size = 10),
          plot.title      = element_text(hjust = 0.5, face = "italic"))
}

# ── A) SAMap joint UMAP ───────────────────────────────────────────────────────
cat("\n=== 9A. SAMap joint UMAP (Gast + Ryum) ===\n")

SAMAP_UMAP <- file.path(BASE_DIR, "SAMap", "pairwise_comps",
                        "gfas_ryum", "figures", "umap_coords.csv")

if (file.exists(SAMAP_UMAP)) {
  umap_samap <- add_plot_group(
    read.csv(SAMAP_UMAP, stringsAsFactors = FALSE),
    species_col  = "species",
    celltype_col = "cell_type"
  )
  p_samap <- make_overlay_plot(umap_samap,
    "G. fascicularis + R. yuma - SAMap joint UMAP\n(equivalent cell types co-cluster)")
  ggsave(file.path(OUT_DIR, "UMAP_Gfas_Ryum_samap.pdf"), p_samap, width = 9, height = 7)
  ggsave(file.path(OUT_DIR, "UMAP_Gfas_Ryum_samap.png"), p_samap, width = 9, height = 7, dpi = 150)
  cat("  Saved: UMAP_Gfas_Ryum_samap.pdf + .png\n")
} else {
  cat("  SKIPPED: umap_coords.csv not found.",
      "Run SAMap/pairwise_comps/gfas_ryum/run_samap_gfas_ryum.py first.\n")
}

# ── B) Seurat merged UMAP ─────────────────────────────────────────────────────
cat("\n=== 9B. Seurat merged UMAP (Gast + Ryum) ===\n")

merged <- merge(gfas_seurat, y = ryum_seurat, add.cell.ids = c("gfas", "ryum"))
DefaultAssay(merged) <- "RNA"
merged <- JoinLayers(merged)
merged <- NormalizeData(merged, verbose = FALSE)
merged <- FindVariableFeatures(merged, nfeatures = 3000, verbose = FALSE)
merged <- ScaleData(merged, verbose = FALSE)
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)
merged <- RunUMAP(merged, dims = 1:25, verbose = FALSE)

umap_seurat        <- as.data.frame(Embeddings(merged, "umap"))
colnames(umap_seurat) <- c("UMAP1", "UMAP2")
umap_seurat$species   <- merged$species
umap_seurat$cell_type <- as.character(merged$cell_type)
umap_seurat <- add_plot_group(umap_seurat,
                              species_col  = "species",
                              celltype_col = "cell_type")

p_seurat <- make_overlay_plot(umap_seurat,
  "G. fascicularis + R. yuma - Seurat merged UMAP\n(common genes, species-level structure)")
ggsave(file.path(OUT_DIR, "UMAP_Gfas_Ryum_seurat.pdf"), p_seurat, width = 9, height = 7)
ggsave(file.path(OUT_DIR, "UMAP_Gfas_Ryum_seurat.png"), p_seurat, width = 9, height = 7, dpi = 150)
cat("  Saved: UMAP_Gfas_Ryum_seurat.pdf + .png\n")
