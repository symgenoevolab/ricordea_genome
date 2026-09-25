#!/bin/bash
# =============================================================================
# 01_cellranger_mapping.sh
# CellRanger reference building and read alignment for G. fascicularis and R. yuma
#
# Author : Victor Pinon-Gonzalez
# Date   : May 2026
#
# Overview
# --------
# G. fascicularis (Gast): de novo Trinity genome-guided transcriptome used as
#   "genome" for CellRanger. A synthetic GTF is generated from the FASTA
#   headers so that every transcript is treated as a single-exon gene.
#   Run inside a Singularity container that provides cellranger v10.
#
# R. yuma (Ryum): chromosomal assembly with BRAKER2 gene models. The BRAKER
#   GTF is extended at 3' ends with GeneExt (to capture 10x read coverage in
#   UTR regions), mitochondrial gene models from MitoFinder are merged in, and
#   the combined annotation is used for mkref.
#
# Sequencing: 10x Chromium v3, paired-end.
#   Gast library: XCSC25-SK01   Ryum library: XCSC25-SK02
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# User configuration — set these paths before running
# ---------------------------------------------------------------------------
CELLRANGER_SIF="path/to/cellranger_v10.sif"   # Singularity image (Gast only)

# G. fascicularis inputs
GFAS_TRANSCRIPTOME="path/to/filtered_transcripts.fasta"  # Trinity okayset (from EvidentialGene), expression-filtered
GFAS_FASTQ="path/to/XCSC25-SK01"                         # directory with R1/R2 FASTQs
GFAS_OUTDIR="path/to/gfas_cellranger"                    # output directory

# R. yuma inputs
RYUM_GENOME="path/to/Ryuma_genome.fna"                   # unmasked genome assembly
RYUM_GENEEXT_GTF="path/to/Ryum_geneext.gtf"              # GeneExt-extended annotation (see below)
RYUM_FASTQ="path/to/XCSC25-SK02"
RYUM_OUTDIR="path/to/ryum_cellranger"

THREADS=16

# =============================================================================
# G. fascicularis — Trinity transcriptome reference
# =============================================================================
echo "=== G. fascicularis: building CellRanger reference from Trinity transcriptome ==="

GFAS_REF="${GFAS_OUTDIR}/cellranger_ref"
mkdir -p "${GFAS_REF}"

# Step 1: Generate a synthetic GTF — one single-exon gene per transcript.
# CellRanger requires genome FASTA + GTF; for a transcriptome we treat each
# transcript as a chromosome and its full length as a single exon.
echo "  Generating synthetic GTF from transcript FASTA headers..."
python3 - <<PYEOF
fa  = "${GFAS_TRANSCRIPTOME}"
gtf = "${GFAS_REF}/trinity_transcripts.gtf"

lengths = {}
current = None
with open(fa) as fh:
    for line in fh:
        line = line.rstrip()
        if line.startswith(">"):
            current = line[1:].split()[0]
            lengths[current] = 0
        else:
            lengths[current] += len(line)

with open(gtf, "w") as out:
    for tid, length in lengths.items():
        for feat in ("gene", "transcript", "exon"):
            out.write(
                f'{tid}\tTrinity\t{feat}\t1\t{length}\t.\t+\t.\t'
                f'gene_id "{tid}"; transcript_id "{tid}";\n'
            )

print(f"Written {len(lengths)} transcript entries to {gtf}")
PYEOF

# Step 2: mkref
echo "  Running cellranger mkref for G. fascicularis..."
singularity exec "${CELLRANGER_SIF}" cellranger mkref \
    --genome=Gfascicularis_trinity_ref \
    --fasta="${GFAS_TRANSCRIPTOME}" \
    --genes="${GFAS_REF}/trinity_transcripts.gtf" \
    --nthreads="${THREADS}" \
    --output-dir="${GFAS_REF}"

# Step 3: count
echo "  Running cellranger count for G. fascicularis (XCSC25-SK01)..."
mkdir -p "${GFAS_OUTDIR}/cellranger_count"
cd "${GFAS_OUTDIR}/cellranger_count"

singularity exec "${CELLRANGER_SIF}" cellranger count \
    --id="gfas_count" \
    --transcriptome="${GFAS_REF}" \
    --fastqs="${GFAS_FASTQ}" \
    --sample="XCSC25-SK01" \
    --localcores="${THREADS}" \
    --localmem=128 \
    --chemistry=10xv3

echo "  G. fascicularis count complete."

# =============================================================================
# R. yuma — chromosomal assembly reference with GeneExt-extended annotation
# =============================================================================
echo ""
echo "=== R. yuma: building CellRanger reference ==="
mkdir -p "${RYUM_OUTDIR}"

# The GTF used for mkref (Ryum_geneext.gtf) was produced by:
#
#   1. fix_mitofinder_gff.sh
#      Converts MitoFinder GFF3 output to GTF format compatible with BRAKER.
#      Adds gene / transcript / exon features; sets source to "MitoFinder".
#
#   2. Merge BRAKER + MitoFinder annotations:
#         cat braker.gtf mitofinder_genes.gtf > merged_Ryuma.gtf
#
#   3. GeneExt 3'-UTR extension:
#         geneext merged_Ryuma.gtf \
#             --bam <10x_BAM> \
#             --genome Ryuma_genome.fna \
#             --output Ryum_geneext.gtf
#
#   The final GTF contains BRAKER + MitoFinder genes with extended 3' UTRs
#   to maximise read capture during CellRanger count.

RYUM_REF="${RYUM_OUTDIR}/cellranger_ref"
echo "  Running cellranger mkref for R. yuma..."
cellranger mkref \
    --genome=Ryuma_genome_ref \
    --fasta="${RYUM_GENOME}" \
    --genes="${RYUM_GENEEXT_GTF}" \
    --nthreads="${THREADS}" \
    --output-dir="${RYUM_REF}"

echo "  Running cellranger count for R. yuma (XCSC25-SK02)..."
mkdir -p "${RYUM_OUTDIR}/cellranger_count"
cd "${RYUM_OUTDIR}/cellranger_count"

cellranger count \
    --id="ryum_count" \
    --transcriptome="${RYUM_REF}" \
    --fastqs="${RYUM_FASTQ}" \
    --sample="XCSC25-SK02" \
    --localcores="${THREADS}" \
    --localmem=128 \
    --chemistry=10xv3

echo "  R. yuma count complete."
echo ""
echo "Count output:"
echo "  Gfas: ${GFAS_OUTDIR}/cellranger_count/gfas_count/outs/"
echo "  Ryum: ${RYUM_OUTDIR}/cellranger_count/ryum_count/outs/"
