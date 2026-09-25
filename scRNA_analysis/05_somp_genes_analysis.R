###############################################################################
# 05_somp_analysis.R
# SOMP (skeletal organic matrix protein) expression analysis across cell types
#
# Author : Victor Pinon-Gonzalez
# Date   : May 2026
#
# Input:  seurat_gfas_ryum/Gast/gfas_seurat_annotated.rds    (from 03_seurat_analysis.R)
#         seurat_gfas_ryum/Ryum/ryuma_seurat_annotated.rds   (from 03_seurat_analysis.R)
#         sc_RNA_public_data/Amur/Amur.seurat_ctrl_annotated.rds
#         sc_RNA_public_data/Spis/Spis.seurat_final.rds
#         sc_RNA_public_data/Amil/Amil.seurat_final.rds
#         sc_RNA_public_data/Opat/Ocupat.seurat_final.rds
#         data/somp_gene_id_keys.tsv  (from SOMPs annotation see method section)
# Output: results_calicoblast/somp_heatmap_all_species.pdf
#         results_calicoblast/somp_expression_long.tsv
#
# Run: Rscript code/05_somp_analysis.R
#
# Overview
# --------
# 1. Load SOMP gene table + all 6 Seurat objects
# 2. Extract mean SCT expression and % expressing cells per gene × cell type
# 3. Aggregate to broad cell types (median)
# 4. Normalize per gene: log2(mean_expr / median_nonzero + 1); flag top cell type
# 5. ggplot heatmap (faceted by species, normalized expression 0-4)
# 6. Save long-format expression table (somp_name | species | sc_gene_id | top_cell_type | mean_expression | pct | log2FC1)
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
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

SOMP_MAP <- file.path(BASE_DIR, "data", "somp_gene_id_keys.tsv")
OUT_DIR  <- file.path(BASE_DIR, "results_calicoblast")

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
# ID conversion functions (same as scripts 02 and 03)
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

# =============================================================================
# Broad cell type mapping
# =============================================================================
map_cell_type <- function(ct) {
  ct_lower <- tolower(ct)
  dplyr::case_when(
    grepl("calicoblast",                  ct_lower) ~ "Calicoblast",
    grepl("epiderm",                      ct_lower) ~ "Epidermis",
    grepl("alga.*host|symbio|zoox",       ct_lower) ~ "Algae_hosting_cells",
    grepl("gastroderm",                   ct_lower) ~ "Gastrodermis",
    grepl("gland",                        ct_lower) ~ "Gland cells",
    grepl("neuron|neurosecret",           ct_lower) ~ "Neurons",
    grepl("cnidocyte|cnidocyt",           ct_lower) ~ "Cnidocytes",
    grepl("immune|macrophage",            ct_lower) ~ "Immune cells",
    grepl("digest",                       ct_lower) ~ "Digestive filaments",
    grepl("muscle",                       ct_lower) ~ "Muscle cells",
    grepl("progenitor",                   ct_lower) ~ "Progenitors",
    grepl("germline|oocyte",              ct_lower) ~ "Germline",
    grepl("unknown|^na[0-9]",             ct_lower) ~ "Unknown/NA",
    TRUE                                             ~ ct
  )
}

# =============================================================================
# Extract mean expression and % expressing cells per gene × cell type
# =============================================================================
extract_gene_celltype_expr <- function(obj, genes, sp_label) {
  assay         <- if ("SCT" %in% names(obj@assays)) "SCT" else "RNA"
  genes_present <- intersect(genes, rownames(obj[[assay]]))
  if (!length(genes_present)) {
    warning(sp_label, ": no matching gene IDs"); return(NULL)
  }
  mat  <- GetAssayData(obj, assay = assay, layer = "data")[genes_present, , drop = FALSE]
  cts  <- obj@meta.data[[get_ct_col(obj)]]
  bind_rows(lapply(unique(cts), function(ct) {
    idx <- which(cts == ct)
    sub <- mat[, idx, drop = FALSE]
    data.frame(species   = sp_label,
               gene      = genes_present,
               cell_type = ct,
               n_cells   = length(idx),
               mean_expr = rowMeans(sub),
               pct_expr  = rowMeans(sub > 0) * 100,
               stringsAsFactors = FALSE)
  }))
}

###############################################################################
# 1. Load SOMP table
###############################################################################
cat("Loading SOMP gene table...\n")
somp_map <- read.delim(SOMP_MAP, sep = "\t", stringsAsFactors = FALSE,
                       check.names = FALSE)
cat(sprintf("  %d SOMP genes\n", nrow(somp_map)))

SP_SC_COLS <- c(gfas = "sc_gfas", amur = "sc_amur", spis = "sc_spis",
                amil = "sc_amil", opat = "sc_opat", ryum = "sc_ryum")
SP_FULL    <- c(gfas = "G. fascicularis", amur = "A. muricata", spis = "S. pistillata",
                amil = "A. millepora", opat = "O. patagonica", ryum = "R. yuma")

seurat_list <- list(gfas = gfas_seurat, amur = amur_seurat, spis = spis_rds,
                    amil = amil_rds,    opat = opat_rds,    ryum = ryum_seurat)

gene_name_map <- bind_rows(lapply(names(SP_SC_COLS), function(sp) {
  sc_col <- SP_SC_COLS[sp]
  if (!sc_col %in% names(somp_map)) return(NULL)
  data.frame(sc_id = somp_map[[sc_col]], gene_name = somp_map$gene_name,
             stringsAsFactors = FALSE)
})) %>% filter(!is.na(sc_id), sc_id != "", sc_id != "NA") %>%
  distinct(sc_id, gene_name)

###############################################################################
# 2. Extract expression per species
###############################################################################
cat("Extracting per-gene, per-cell-type mean expression...\n")
all_expr <- list()
for (sp in names(seurat_list)) {
  sc_col <- SP_SC_COLS[sp]
  if (!sc_col %in% names(somp_map)) next
  genes <- unique(somp_map[[sc_col]])
  genes <- genes[!is.na(genes) & genes != "" & genes != "NA"]
  expr  <- extract_gene_celltype_expr(seurat_list[[sp]], genes, SP_FULL[sp])
  if (!is.null(expr)) all_expr[[sp]] <- expr
}
expr_df <- bind_rows(all_expr) %>%
  left_join(gene_name_map, by = c("gene" = "sc_id")) %>%
  mutate(gene_name = coalesce(gene_name, gene))

###############################################################################
# 3. Broad cell types, aggregate with median
###############################################################################
expr_df <- expr_df %>%
  mutate(cell_type_broad = map_cell_type(cell_type)) %>%
  filter(cell_type_broad != "Unknown/NA")

heatmap_data <- expr_df %>%
  group_by(gene, gene_name, species, cell_type_broad) %>%
  summarise(mean_expr = median(mean_expr, na.rm = TRUE),
            pct_expr  = mean(pct_expr,   na.rm = TRUE),
            n_cells   = sum(n_cells), .groups = "drop")

###############################################################################
# 4. Normalize expression and flag top-expressing cell type
# Per-gene normalization: divide by the median of NONZERO groups (across all
# species x cell types), then log2(fc + 1). This baseline = "typical level when
# detected" — avoids white-heatmap collapse for sparsely expressed SOMP genes.
###############################################################################
heatmap_data <- heatmap_data %>%
  group_by(gene_name) %>%
  mutate(median_nonzero = {
           nz <- mean_expr[mean_expr > 0]
           if (length(nz) > 0) median(nz, na.rm = TRUE) else NA_real_
         },
         fc      = ifelse(!is.na(median_nonzero) & median_nonzero > 0,
                          mean_expr / median_nonzero, 0),
         expr_fc = ifelse(mean_expr == 0 | is.na(mean_expr), 0,
                          log2(fc + 1))) %>%
  ungroup() %>%
  group_by(gene_name, species) %>%
  mutate(expr_rank = rank(-mean_expr, ties.method = "first"),
         is_top1   = (expr_rank == 1),
         label     = ifelse(is_top1 & mean_expr > 0, "*", "")) %>%
  ungroup()

###############################################################################
# 5. Heatmap
###############################################################################
cat("Building SOMP expression heatmap...\n")

CT_ORDER <- c("Calicoblast", "Epidermis", "Gastrodermis", "Algae_hosting_cells", "Neurons")
SP_ORDER <- c("A. millepora", "A. muricata", "G. fascicularis",
              "S. pistillata", "O. patagonica", "R. yuma")
GENE_ORDER <- c("GXN2", "USOM6", "CPP1/2", "GXN", "USOM5", "USOM1", "CA", "USOM2", "CTXL", "USOM3", 
                "ASOMP", "SAAR1", "SAAR2", "SLC4y", "CDP", "MLP", "CADN", "CARP3", "USOM7", "USOM8",
                "MLRP1/2", "ECT", "ZPP", "HEPHL", "PK1L", "PCDL", "CARP1",  "COA", "ELP", "FP", "SLC4B" )

heatmap_data <- heatmap_data %>%
  filter(cell_type_broad %in% CT_ORDER) %>%
  mutate(cell_type_broad = factor(cell_type_broad, levels = CT_ORDER),
         species         = factor(species,         levels = SP_ORDER),
         gene_name       = factor(gene_name,       levels = GENE_ORDER))

p_heatmap <- ggplot(heatmap_data,
                    aes(x = cell_type_broad, y = gene_name, fill = expr_fc)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_grid(. ~ species, scales = "fixed", space = "fixed") +
  coord_fixed() +
  scale_fill_gradient(
    low    = "white",
    high   = "darkred",
    name   = "Normalized\nexpression",
    limits = c(0, 4),
    oob    = scales::squish,
    na.value = "grey90"
  ) +
  theme_classic(base_size = 10) +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y   = element_text(size = 8),
        strip.text    = element_text(face = "italic", size = 9),
        panel.spacing = unit(0.5, "lines")) +
  labs(x = NULL, y = NULL,
       title   = "SOMP gene expression across cell types and species")

ggsave(file.path(OUT_DIR, "somp_heatmap_all_species.pdf"),
       p_heatmap, width = 18, height = 10)
cat("  Saved: somp_heatmap_all_species.pdf\n")

###############################################################################
# 6. Long-format SOMP expression table (all genes x cell types x species)
###############################################################################
cat("Building long-format SOMP expression table...\n")

somp_long <- heatmap_data %>%
  filter(is_top1) %>%
  select(
    somp_name        = gene_name,
    species,
    sc_gene_id       = gene,
    top_cell_type    = cell_type_broad,
    mean_expression  = mean_expr,
    pct              = pct_expr,
    log2FC1          = expr_fc
  ) %>%
  mutate(across(where(is.numeric), ~ round(., 4))) %>%
  arrange(species, somp_name)

write.table(somp_long,
            file.path(OUT_DIR, "somp_expression_long.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
cat("  Saved: somp_expression_long.tsv |", nrow(somp_long), "rows\n")

# Full heatmap data: all genes x all cell types x all species (used to build somp_heatmap_all_species.pdf)
somp_heatmap_long <- heatmap_data %>%
  select(
    somp_name       = gene_name,
    species,
    cell_type       = cell_type_broad,
    mean_expression = mean_expr,
    pct             = pct_expr,
    normalized_expr = expr_fc,
    is_top_cell_type = is_top1
  ) %>%
  mutate(across(where(is.numeric), ~ round(., 4))) %>%
  arrange(species, somp_name, cell_type)

write.table(somp_heatmap_long,
            file.path(OUT_DIR, "somp_heatmap_data_long.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
cat("  Saved: somp_heatmap_data_long.tsv |", nrow(somp_heatmap_long), "rows\n")
