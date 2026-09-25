###############################################################################
# 02_Amur_seurat_obj_reconstruction.R
# Reconstruct A. muricata annotated Seurat object from Han et al. 2025
# (Communications Biology 8:652)
#
# Author : Victor Pinon-Gonzalez
# Date   : May 2026
#
# Input:  sc_RNA_public_data/Amur/control_folder/ 
#         sc_RNA_public_data/Amur/day21_folder/
# Output: sc_RNA_public_data/Amur/Amur.seurat_ctrl_annotated.rds
#
# Run before 03_seurat_analysis.R:
#   Rscript code/02_Amur_seurat_obj_reconstruction.R
#
# Workflow:
#   1. Load ctrl + day21 with ReadMtx (feature.column = 2)
#   2. QC: nFeature_RNA 200–2500 (no mt filter — data already clean)
#   3. SCTransform integration (ctrl + day21, 3000 features)
#   4. PCA → FindNeighbors (dims 1:30) → FindClusters (res 0.1) → UMAP
#   5. Annotate with cluster_to_celltype; harmonise to cross-species palette
#   6. Subset to control cells, used for this analysis (orig.ident == "control")
#   8. Save RDS
#
# Cluster labels (Han et al. 2025 TSNE annotation):
#   0  → NA1 (Originally: Unassigned cell cluster1)
#   1  → Gastrodermis
#   2  → NA2 (Originally: Unassigned cell cluster2)
#   3  → Neurons
#   4  → Calicoblast
#   5  → Alga-hosting cells
#   6  → Epidermis
#   7  → Cnidocytes
#   8  → Gland cells
#   9  → Progenitors
#   10 → Immune_cells
#   11 → Muscle cells
#   12 → Digestive filaments
###############################################################################

suppressPackageStartupMessages(library(Seurat))

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

AMUR_CTRL_DIR  <- file.path(BASE_DIR, "sc_RNA_public_data", "Amur", "control_folder")
AMUR_DAY21_DIR <- file.path(BASE_DIR, "sc_RNA_public_data", "Amur", "day21_folder")
OUT_RDS        <- file.path(BASE_DIR, "sc_RNA_public_data", "Amur", "Amur.seurat_ctrl_annotated.rds")

# =============================================================================
# 1. Load data
# =============================================================================
cat("Loading A. muricata matrices...\n")

control_data <- ReadMtx(
  mtx      = file.path(AMUR_CTRL_DIR,  "matrix.mtx.gz"),
  features = file.path(AMUR_CTRL_DIR,  "features.tsv.gz"),
  cells    = file.path(AMUR_CTRL_DIR,  "barcodes.tsv.gz"),
  feature.column = 2
)
control <- CreateSeuratObject(counts = control_data, project = "control", min.cells = 3)
control <- RenameCells(control, add.cell.id = "Control")
control[["species"]] <- "A. muricata"

day21_data <- ReadMtx(
  mtx      = file.path(AMUR_DAY21_DIR, "matrix.mtx.gz"),
  features = file.path(AMUR_DAY21_DIR, "features.tsv.gz"),
  cells    = file.path(AMUR_DAY21_DIR, "barcodes.tsv.gz"),
  feature.column = 2
)
day21 <- CreateSeuratObject(counts = day21_data, project = "day21", min.cells = 3)
day21 <- RenameCells(day21, add.cell.id = "Day21")

cat(sprintf("  Control cells before QC: %d\n", ncol(control)))
cat(sprintf("  Day21  cells before QC:  %d\n", ncol(day21)))

# =============================================================================
# 2. QC
# =============================================================================
control <- subset(control, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
day21   <- subset(day21,   subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
cat(sprintf("  Control cells after QC:  %d\n", ncol(control)))
cat(sprintf("  Day21  cells after QC:   %d\n", ncol(day21)))

# =============================================================================
# 3. SCTransform + integration (ctrl + day21, 3000 features)
# =============================================================================
cat("SCTransform + integration...\n")
list1    <- lapply(list(control, day21), SCTransform, vst.flavor = "v2", verbose = FALSE)
features <- SelectIntegrationFeatures(object.list = list1, nfeatures = 3000)
list2    <- PrepSCTIntegration(object.list = list1, anchor.features = features)
anchors  <- FindIntegrationAnchors(object.list = list2, normalization.method = "SCT",
                                   anchor.features = features)
object   <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# =============================================================================
# 4. PCA + clustering + UMAP  (dims 1:30, resolution 0.1 — same as Han et al.)
# =============================================================================
cat("PCA / clustering / UMAP...\n")
DefaultAssay(object) <- "integrated"
object <- RunPCA(object,      verbose = FALSE)
object <- FindNeighbors(object, dims = 1:30, verbose = FALSE)
object <- FindClusters(object,  resolution = 0.1, verbose = FALSE)
object <- RunUMAP(object,     dims = 1:30, verbose = FALSE)

cat("\nCluster sizes (integrated, res = 0.1):\n")
print(table(Idents(object)))

# =============================================================================
# 5. Cell type annotation (Han et al. 2025; see header for full mapping)
# =============================================================================
cluster_to_celltype <- c(
  "0"  = "NA1",
  "1"  = "Gastrodermis",
  "2"  = "NA2",
  "3"  = "Neurons",
  "4"  = "Calicoblast",
  "5"  = "Algae_hosting_cells-like",
  "6"  = "Epidermis",
  "7"  = "Cnidocytes",
  "8"  = "Gland_cells",
  "9"  = "Progenitors",
  "10" = "Immune_cells",
  "11" = "Muscle_cells",
  "12" = "Digestive_filaments"
)

cluster_vec  <- object$seurat_clusters
celltype_vec <- cluster_to_celltype[as.character(cluster_vec)]
names(celltype_vec) <- names(cluster_vec)
object <- AddMetaData(object, metadata = celltype_vec, col.name = "cell_type")

# =============================================================================
# 6. Subset to control cells only
# =============================================================================
cat("\nSubsetting to control cells (orig.ident == 'control')...\n")
amur_ctrl <- subset(object, subset = orig.ident == "control")
cat(sprintf("  Control cells retained: %d\n", ncol(amur_ctrl)))

# =============================================================================
# 8. Save
# =============================================================================
saveRDS(amur_ctrl, OUT_RDS)
cat("\nSaved:", OUT_RDS, "\n")
