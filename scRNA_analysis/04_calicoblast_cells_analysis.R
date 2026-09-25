###############################################################################
# 04_calicoblast_cells_analysis.R
# Cross-species cell type similarity and calicoblast marker extraction
#
# Author : Victor Pinon-Gonzalez
# Date   : May 2026
#
# Input:  seurat_gfas_ryum/Gast/gfas_seurat_annotated.rds    (from 03_seurat_analysis.R)
#         seurat_gfas_ryum/Ryum/ryuma_seurat_annotated.rds   (from 03_seurat_analysis.R)
#         sc_RNA_public_data/Amur/Amur.seurat_ctrl_annotated.rds   (02_Amur_seurat_obj_reconstruction.R)
#         sc_RNA_public_data/Spis/Spis.seurat_final.rds  (public source)
#         sc_RNA_public_data/Amil/Amil.seurat_final.rds  (public source)
#         sc_RNA_public_data/Opat/Ocupat.seurat_final.rds (public source)
#         data/somp_rbh_output_for_scRNA.tsv (from reciprocal best hit analysis using 6 species proteomes see paper method section)
# Output: results_calicoblast/heatmap_cell_type_similarity.pdf
#         results_calicoblast/pairwise_similarity.tsv
#         results_calicoblast/calicoblast_markers.tsv
#
# Run: Rscript code/04_calicoblast_cells_analysis.R
#
# Overview
# --------
# 1. Load all 6 Seurat objects + ortholog table (to integrate scRNA objects)
# 2. Build combined expression matrix
# 3. All-vs-all Spearman correlation → heatmap (Ward.D2)
# 4. FindMarkers (calicoblast vs all, Wilcoxon) for 5 coral species
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(pheatmap)
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

ORTHOLOG_TABLE <- file.path(BASE_DIR, "data", "somp_rbh_output_for_scRNA.tsv")
OUT_DIR        <- file.path(BASE_DIR, "results_calicoblast")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Load Seurat objects
# ---------------------------------------------------------------------------
cat("Loading Seurat objects...\n")
gfas_seurat <- readRDS(file.path(BASE_DIR, "seurat_gfas_ryum", "Gfas", "gfas_seurat_annotated.rds"))
cat("  G. fascicularis loaded\n")
ryum_seurat <- readRDS(file.path(BASE_DIR, "seurat_gfas_ryum", "Ryum", "ryuma_seurat_annotated.rds"))
cat("  R. yuma loaded\n")
amur_seurat <- readRDS(file.path(BASE_DIR, "sc_RNA_public_data", "Amur", "Amur.seurat_ctrl_annotated.rds"))
cat("  A. muricata loaded\n")
spis_rds    <- readRDS(file.path(BASE_DIR, "sc_RNA_public_data", "Spis", "Spis.seurat_final.rds"))
cat("  S. pistillata loaded\n")
amil_rds    <- readRDS(file.path(BASE_DIR, "sc_RNA_public_data", "Amil", "Amil.seurat_final.rds"))
cat("  A. millepora loaded\n")
opat_rds    <- readRDS(file.path(BASE_DIR, "sc_RNA_public_data", "Opat", "Ocupat.seurat_final.rds"))
cat("  O. patagonica loaded\n")

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

get_ct_col <- function(obj) {
  candidates <- c("cell_type", "celltype", "CellType", "cell_type_broad",
                  "seurat_clusters", "ident")
  found <- candidates[candidates %in% names(obj@meta.data)]
  if (length(found)) return(found[1])
  obj[["._ct"]] <- as.character(Idents(obj)); "._ct"
}

###############################################################################
# 1. Load ortholog table
###############################################################################
cat("Loading 6-way ortholog table...\n")
ortho <- read.delim(ORTHOLOG_TABLE, sep = "\t", stringsAsFactors = FALSE,
                    check.names = FALSE)
ortho <- ortho[ortho$confidence %in% c("HIGH", "MEDIUM"), ]
cat(sprintf("  %d ortholog groups (HIGH + MEDIUM confidence)\n", nrow(ortho)))

id_converters <- list(gfas = sc_gfas, amur = sc_amur, spis = sc_spis,
                      amil = sc_amil, opat = sc_opat, ryum = sc_ryum)
for (sp in names(id_converters)) {
  raw <- sub("\\.p[0-9]+$", "", ortho[[sp]])
  ortho[[paste0("sc_", sp)]] <- id_converters[[sp]](raw)
}

sc_cols       <- paste0("sc_", names(id_converters))
complete_mask <- apply(ortho[, sc_cols], 1, function(r) all(nzchar(r) & r != "NA"))
ortho_complete <- ortho[complete_mask, ]
cat(sprintf("  Complete ortholog groups (all 6 species): %d\n", nrow(ortho_complete)))

###############################################################################
# 2. Build combined expression matrix
###############################################################################

# Average expression per cell type (log2)
calculate_cell_type_avg <- function(seurat_obj, species_name, assay = "SCT") {
  if (!assay %in% Assays(seurat_obj)) {
    assay <- DefaultAssay(seurat_obj)
    cat("  Using assay:", assay, "for", species_name, "\n")
  }
  avg_expr <- AverageExpression(seurat_obj, group.by = "cell_type",
                                assays = assay)[[assay]]
  avg_log <- log2(as.matrix(avg_expr) + 1)
  colnames(avg_log) <- paste0(species_name, "_", colnames(avg_log))
  cat(sprintf("  %s: %d cell types | %d genes\n",
              species_name, ncol(avg_log), nrow(avg_log)))
  avg_log
}

cat("\nCalculating average expression per cell type...\n")
gfas_avg <- calculate_cell_type_avg(gfas_seurat, "GFAS", "SCT")
amur_avg <- calculate_cell_type_avg(amur_seurat, "AMUR", "SCT")
spis_avg <- calculate_cell_type_avg(spis_rds,    "SPIS", "SCT")
amil_avg <- calculate_cell_type_avg(amil_rds,    "AMIL", "SCT")
opat_avg <- calculate_cell_type_avg(opat_rds,    "OPAT", "SCT")
ryum_avg <- calculate_cell_type_avg(ryum_seurat, "RYUM", "SCT")

# Combine species matrices via ortholog table
create_combined_matrix <- function(ortho_complete,
                                   gfas_avg, amur_avg, spis_avg,
                                   ryum_avg, amil_avg, opat_avg) {
  gfas_expr <- gfas_avg[match(ortho_complete$sc_gfas, rownames(gfas_avg)), , drop = FALSE]
  amur_expr <- amur_avg[match(ortho_complete$sc_amur, rownames(amur_avg)), , drop = FALSE]
  spis_expr <- spis_avg[match(ortho_complete$sc_spis, rownames(spis_avg)), , drop = FALSE]
  ryum_expr <- ryum_avg[match(ortho_complete$sc_ryum, rownames(ryum_avg)), , drop = FALSE]
  amil_expr <- amil_avg[match(ortho_complete$sc_amil, rownames(amil_avg)), , drop = FALSE]
  opat_expr <- opat_avg[match(ortho_complete$sc_opat, rownames(opat_avg)), , drop = FALSE]

  rownames(gfas_expr) <- ortho_complete$coral_og_id
  rownames(amur_expr) <- ortho_complete$coral_og_id
  rownames(spis_expr) <- ortho_complete$coral_og_id
  rownames(ryum_expr) <- ortho_complete$coral_og_id
  rownames(amil_expr) <- ortho_complete$coral_og_id
  rownames(opat_expr) <- ortho_complete$coral_og_id

  combined <- cbind(gfas_expr, amur_expr, spis_expr, ryum_expr, amil_expr, opat_expr)
  cat(sprintf("  Combined matrix: %d genes x %d cell types\n",
              nrow(combined), ncol(combined)))
  combined
}

combined_expr <- create_combined_matrix(
  ortho_complete, gfas_avg, amur_avg, spis_avg, ryum_avg, amil_avg, opat_avg
)

# Z-score normalization (batch correction)
normalize_per_species <- function(combined_matrix) {
  cat("\nZ-score normalization (batch correction)...\n")
  species_cols <- list(
    GFAS = grep("^GFAS_", colnames(combined_matrix)),
    AMUR = grep("^AMUR_", colnames(combined_matrix)),
    SPIS = grep("^SPIS_", colnames(combined_matrix)),
    RYUM = grep("^RYUM_", colnames(combined_matrix)),
    AMIL = grep("^AMIL_", colnames(combined_matrix)),
    OPAT = grep("^OPAT_", colnames(combined_matrix))
  )

  normalized <- combined_matrix
  for (sp in names(species_cols)) {
    cols <- species_cols[[sp]]
    if (length(cols) > 1) {
      sp_data   <- combined_matrix[, cols]
      row_means <- rowMeans(sp_data, na.rm = TRUE)
      row_sds   <- apply(sp_data, 1, sd, na.rm = TRUE)
      row_sds[row_sds == 0] <- 1
      normalized[, cols] <- (sp_data - row_means) / row_sds
      cat(sprintf("  %s: normalized %d cell types\n", sp, length(cols)))
    }
  }

  valid_rows <- complete.cases(normalized) &
                !apply(normalized, 1, function(x) any(is.infinite(x)))
  cat(sprintf("\n  Genes after filtering: %d of %d\n",
              sum(valid_rows), nrow(normalized)))
  normalized[valid_rows, ]
}

combined_normalized <- normalize_per_species(combined_expr)

cat("\nCell types per species:\n")
for (sp in c("GFAS", "AMUR", "SPIS", "RYUM", "AMIL", "OPAT")) {
  cat(sprintf("  %s: %d\n", sp,
              sum(grepl(paste0("^", sp, "_"), colnames(combined_normalized)))))
}
cat(sprintf("Total cell types: %d | Orthologous genes: %d\n",
            ncol(combined_normalized), nrow(combined_normalized)))

###############################################################################
# 3. All-vs-all Spearman correlation
###############################################################################
cat("Computing all-vs-all Spearman correlation...\n")

calculate_pairwise_similarity <- function(expr_matrix) {
  cor_matrix <- cor(expr_matrix, method = "spearman")
  cat("  Correlation matrix:", nrow(cor_matrix), "x", ncol(cor_matrix), "\n")

  cor_df <- as.data.frame(as.table(cor_matrix))
  colnames(cor_df) <- c("celltype1", "celltype2", "similarity")
  cor_df <- cor_df %>%
    filter(as.character(celltype1) != as.character(celltype2)) %>%
    rowwise() %>%
    mutate(pair = paste(sort(c(as.character(celltype1),
                               as.character(celltype2))), collapse = "___")) %>%
    ungroup() %>%
    distinct(pair, .keep_all = TRUE) %>%
    select(-pair) %>%
    mutate(
      species1        = sub("_.*", "", celltype1),
      species2        = sub("_.*", "", celltype2),
      type1           = sub("^[^_]+_", "", celltype1),
      type2           = sub("^[^_]+_", "", celltype2),
      same_species    = species1 == species2,
      comparison_type = ifelse(same_species, "within_species", "between_species")
    ) %>%
    arrange(desc(similarity))

  cat("  Total unique pairs:", nrow(cor_df), "\n")
  cat("  Within-species pairs:", sum(cor_df$same_species), "\n")
  cat("  Between-species pairs:", sum(!cor_df$same_species), "\n")
  list(matrix = cor_matrix, pairs = cor_df)
}

similarity_results <- calculate_pairwise_similarity(combined_normalized)
cor_matrix <- similarity_results$matrix
all_pairs  <- similarity_results$pairs

calico_pairs <- all_pairs %>%
  filter(!same_species,
         grepl("calicoblast", tolower(type1)),
         grepl("calicoblast", tolower(type2)))
cat(sprintf("\n  Inter-calicoblast pairs: %d | mean similarity: %.4f\n",
            nrow(calico_pairs), mean(calico_pairs$similarity)))
cat(sprintf("  Mean overall between-species similarity: %.4f\n",
            mean(all_pairs$similarity[!all_pairs$same_species])))

###############################################################################
# 4. Heatmap
###############################################################################
create_calicoblast_heatmap <- function(cor_matrix, gamma = 0.5) {
  # Clustering uses raw Spearman ρ; display uses signed power transform (gamma)
  clustering_matrix <- cor_matrix
  heatmap_matrix    <- cor_matrix
  diag(heatmap_matrix) <- NA
  scaled_matrix <- sign(heatmap_matrix) * abs(heatmap_matrix)^gamma
  off_diag <- scaled_matrix[!is.na(scaled_matrix)]
  max_abs  <- max(abs(off_diag))

  cat(sprintf("  Raw range: %.3f to %.3f\n",
              min(cor_matrix), max(cor_matrix)))
  cat(sprintf("  Visual gamma: %.2f | Scaled range: %.3f to %.3f\n",
              gamma, -max_abs, max_abs))

  annotation_df <- data.frame(
    Species   = sub("_.*", "", colnames(cor_matrix)),
    row.names = colnames(cor_matrix)
  )
  annotation_df$Highlight <- dplyr::case_when(
    grepl("Calicoblast", rownames(annotation_df), ignore.case = TRUE) ~ "Calicoblast",
    grepl("Epidermis",   rownames(annotation_df), ignore.case = TRUE) ~ "Epidermis",
    TRUE ~ "NA"
  )

  annotation_colors <- list(
    Species = c(
      "GFAS" = "#d53e4f", "AMUR" = "#fc8d59", "SPIS" = "#fee08b",
      "AMIL" = "#e6f598", "OPAT" = "#99d594", "RYUM" = "#3288bd"
    ),
    Highlight = c(
      "Calicoblast" = "#9A88B5", "Epidermis" = "#6F8FBF", "NA" = "white"
    )
  )

  row_dist <- as.dist(1 - clustering_matrix)
  col_dist <- as.dist(1 - clustering_matrix)
  stopifnot(!anyNA(row_dist), !anyNA(col_dist))

  pheatmap::pheatmap(
    scaled_matrix,
    annotation_row    = annotation_df,
    annotation_col    = annotation_df,
    annotation_colors = annotation_colors,
    color  = colorRampPalette(c("#4E79A7", "white", "darkred"))(100),
    breaks = seq(-max_abs, max_abs, length.out = 101),
    main   = paste0("Cell Type Similarity\n",
                    "Signed power scaling (gamma = ", gamma,
                    "), diagonal excluded"),
    fontsize_row = 6, fontsize_col = 6,
    na_col = "grey95",
    border_color = NA,
    clustering_distance_rows = row_dist,
    clustering_distance_cols = col_dist,
    clustering_method = "ward.D2",
    silent = TRUE
  )
}

cat("\nBuilding heatmap...\n")
p_heatmap <- create_calicoblast_heatmap(cor_matrix)

pdf(file.path(OUT_DIR, "heatmap_cell_type_similarity.pdf"), width = 14, height = 13)
print(p_heatmap)
dev.off()
cat("  Saved: heatmap_cell_type_similarity.pdf\n")

write.table(all_pairs, file.path(OUT_DIR, "pairwise_similarity.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved: pairwise_similarity.tsv\n")

###############################################################################
# 5. Calicoblast marker genes (Wilcoxon, p_val_adj < 0.05) — 5 coral species
# R. yuma excluded: no calicoblast population detected.
# Note: gfas/amur use "Calicoblast" (capital); spis/amil/opat use "calicoblast".
###############################################################################
cat("\nExtracting calicoblast markers...\n")

coral_objects <- list(
  list(sp = "gfas", obj = gfas_seurat, ident = "Calicoblast"),
  list(sp = "amur", obj = amur_seurat, ident = "Calicoblast"),
  list(sp = "spis", obj = spis_rds,    ident = "calicoblast"),
  list(sp = "amil", obj = amil_rds,    ident = "calicoblast"),
  list(sp = "opat", obj = opat_rds,    ident = "calicoblast")
)

all_markers <- list()
for (x in coral_objects) {
  DefaultAssay(x$obj) <- "SCT"
  Idents(x$obj)       <- "cell_type"
  x$obj <- PrepSCTFindMarkers(x$obj)

  m <- FindMarkers(x$obj, ident.1 = x$ident, min.pct = 0.10,
                   logfc.threshold = 0.25, only.pos = TRUE, test.use = "wilcox")
  m <- m[m$p_val_adj < 0.05, ]
  m$gene    <- rownames(m)
  m$species <- x$sp
  all_markers[[x$sp]] <- m
  cat(sprintf("  %s: %d markers\n", x$sp, nrow(m)))
}

calicoblast_markers <- do.call(rbind, all_markers)
calicoblast_markers <- calicoblast_markers[, c("species", "gene",
                                               "p_val", "avg_log2FC",
                                               "pct.1", "pct.2", "p_val_adj")]
rownames(calicoblast_markers) <- NULL

write.table(calicoblast_markers,
            file.path(OUT_DIR, "calicoblast_markers.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Saved: calicoblast_markers.tsv (%d rows)\n",
            nrow(calicoblast_markers)))

cat("\nDone. Results in:", OUT_DIR, "\n")
