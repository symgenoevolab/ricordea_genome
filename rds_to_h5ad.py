#!/usr/bin/env python3
"""Convert Seurat v4/v5 RDS file to h5ad for TranscriptFormer inference.

Usage:
    python rds_to_h5ad.py <input.rds> <output.h5ad> [--species gast|ryuma]

Species-specific gene name handling
------------------------------------
Gastreata (--species gast):
    RDS row names: g.TRINITY-GG-24-c4-g2-i3
    FASTA-derived: TRINITY_GG_24_c4_g2_i3  (prefix "g." stripped, "-" -> "_")
    Transformation applied automatically.

Ryuma (--species ryuma):
    RDS row names: g4, g5, g6, ...  already match FASTA-derived IDs.
    No transformation needed.

Count source
------------
SCT@counts (raw integer counts, subset to SCT HVGs) is preferred and used by
default when available.  Falls back to RNA@counts when SCT assay is absent.
Both Corallimorpharia RDS files store SCT Pearson residuals in RNA@counts,
so --use-rna-counts is provided only for debugging.
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io


R_SCRIPT = """\
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly=TRUE)
rds_path      <- args[1]
mtx_path      <- args[2]
genes_path    <- args[3]
barcodes_path <- args[4]
meta_path     <- args[5]
use_rna       <- as.logical(args[6])

cat("Reading RDS:", rds_path, "\\n")
obj <- readRDS(rds_path)
cat("Seurat version:", as.character(slot(obj, "version")), "\\n")
cat("Available assays:", paste(names(slot(obj, "assays")), collapse=", "), "\\n")

use_sct <- !use_rna && ("SCT" %in% names(slot(obj, "assays")))

if (use_sct) {
    cat("Using SCT@counts (raw integer counts)\\n")
    sct <- slot(obj, "assays")[["SCT"]]
    counts_mat <- slot(sct, "counts")
    genes <- rownames(counts_mat)
} else {
    cat("Using RNA@counts\\n")
    rna <- slot(obj, "assays")[["RNA"]]
    rna_class <- class(rna)[1]
    cat("RNA assay class:", rna_class, "\\n")
    if (rna_class == "Assay5") {
        counts_mat <- slot(rna, "layers")[["counts"]]
        if (is.null(counts_mat)) {
            idx <- grep("^counts", names(slot(rna, "layers")))[1]
            counts_mat <- slot(rna, "layers")[[idx]]
        }
        genes <- rownames(slot(rna, "features"))
    } else {
        counts_mat <- slot(rna, "counts")
        genes <- rownames(counts_mat)
    }
}

if (is.null(counts_mat) || length(counts_mat) == 0) {
    stop("No count matrix found.")
}

meta <- slot(obj, "meta.data")
barcodes <- rownames(meta)

cat("Dims:", nrow(counts_mat), "genes x", ncol(counts_mat), "cells\\n")
writeMM(counts_mat, mtx_path)
write.table(genes,    genes_path,    row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(barcodes, barcodes_path, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.csv(meta, meta_path)
cat("Done.\\n")
"""


def transform_gast(gene: str) -> str:
    """g.TRINITY-GG-24-c4-g2-i3  ->  TRINITY_GG_24_c4_g2_i3"""
    if gene.startswith("g."):
        gene = gene[2:]
    return gene.replace("-", "_")


TRANSFORMS = {
    "gast": transform_gast,
    "ryuma": None,
}


def run_r_export(rds_path: Path, tmpdir: Path, use_rna: bool) -> tuple[Path, Path, Path, Path]:
    r_script_path = tmpdir / "export.R"
    mtx_path      = tmpdir / "counts.mtx"
    genes_path    = tmpdir / "genes.txt"
    barcodes_path = tmpdir / "barcodes.txt"
    meta_path     = tmpdir / "metadata.csv"

    r_script_path.write_text(R_SCRIPT)

    cmd = [
        "Rscript", "--vanilla", str(r_script_path),
        str(rds_path), str(mtx_path), str(genes_path), str(barcodes_path), str(meta_path),
        "TRUE" if use_rna else "FALSE",
    ]

    print("Running R export ...")
    result = subprocess.run(cmd, capture_output=False, text=True)
    if result.returncode != 0:
        print("\nR export failed. Make sure SeuratObject and Matrix are installed:")
        print('  Rscript -e \'install.packages(c("SeuratObject","Matrix"))\'')
        sys.exit(1)

    return mtx_path, genes_path, barcodes_path, meta_path


def build_h5ad(
    mtx_path: Path,
    genes_path: Path,
    barcodes_path: Path,
    meta_path: Path,
    output_path: Path,
    gene_transform=None,
    assay: str = "10x 3' v3",
) -> None:
    print("Reading exported files ...")

    # Matrix from R is genes x cells; transpose to cells x genes
    X = scipy.io.mmread(mtx_path).T.tocsr()

    genes = pd.read_csv(genes_path, header=None)[0].values
    barcodes = pd.read_csv(barcodes_path, header=None)[0].values
    meta = pd.read_csv(meta_path, index_col=0)

    print(f"  {X.shape[0]} cells x {X.shape[1]} genes")

    if X.shape[0] != len(barcodes):
        raise ValueError(f"Cell count mismatch: {X.shape[0]} rows vs {len(barcodes)} barcodes.")
    if X.shape[1] != len(genes):
        raise ValueError(f"Gene count mismatch: {X.shape[1]} cols vs {len(genes)} genes.")

    if gene_transform is not None:
        genes_transformed = np.array([gene_transform(g) for g in genes])
        n_changed = int(np.sum(genes_transformed != genes))
        print(f"  Gene name transform applied: {n_changed}/{len(genes)} names changed")
        genes = genes_transformed

    sample_data = X.data[:5000] if X.nnz > 5000 else X.data
    frac_integer = np.mean(np.abs(sample_data - np.round(sample_data)) < 1e-6)
    if frac_integer < 1.0:
        print(
            f"WARNING: {(1-frac_integer)*100:.1f}% of non-zero values are non-integer. "
            "TranscriptFormer works best with raw integer counts. "
            "Consider using --use-rna-counts=False (default) which reads SCT@counts."
        )
    else:
        print(f"  Counts verified as raw integers (checked {len(sample_data)} non-zero values).")

    var = pd.DataFrame({"ensembl_id": genes}, index=genes)
    obs = meta.reindex(barcodes)

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs_names = barcodes
    adata.var_names = genes
    adata.obs["assay"] = assay

    print(f"Writing {output_path} ...")
    adata.write_h5ad(output_path)
    print(f"Done. {X.shape[0]} cells x {X.shape[1]} genes -> {output_path}")
    print("  .var['ensembl_id'] contains gene IDs for TranscriptFormer matching.")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("rds",  help="Input Seurat RDS file")
    parser.add_argument("h5ad", help="Output h5ad file")
    parser.add_argument(
        "--species",
        choices=list(TRANSFORMS.keys()),
        default=None,
        help="Apply species-specific gene name transformation (gast or ryuma)",
    )
    parser.add_argument(
        "--assay",
        default="10x 3' v3",
        help="Sequencing assay label for TranscriptFormer (default: '10x 3' v3'). "
             "Must match a key in the checkpoint's assay_vocab.json.",
    )
    parser.add_argument(
        "--use-rna-counts",
        action="store_true",
        default=False,
        help="Force use of RNA@counts instead of SCT@counts (not recommended for these datasets).",
    )
    args = parser.parse_args()

    rds_path = Path(args.rds).resolve()
    output_path = Path(args.h5ad).resolve()

    if not rds_path.exists():
        print(f"Error: RDS file not found: {rds_path}")
        sys.exit(1)

    gene_transform = TRANSFORMS.get(args.species) if args.species else None
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="rds2h5ad_") as tmpdir:
        tmpdir = Path(tmpdir)
        mtx, genes, barcodes, meta = run_r_export(rds_path, tmpdir, use_rna=args.use_rna_counts)
        build_h5ad(mtx, genes, barcodes, meta, output_path, gene_transform=gene_transform, assay=args.assay)


if __name__ == "__main__":
    main()
