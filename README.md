This repository contains scripts used to generate figures for our Corallimorpharia (_Ricordea yuma_) paper.

Python file:
- rds_to_h5ad.py

  Converts a Seurat v4/v5 RDS file to h5ad for TranscriptFormer inference.

R Markdown files:
- calicoblast_subsampling.Rmd

  Analysis and visualisation of calicoblast marker gene detectability using a subsampling approach.
- topGO_calicoblast_epidermis.Rmd

  Analysis and visualisation of GO enrichment using topGO.
- transcriptformer_calicoblast_centroid.Rmd

  Analysis and visualisation of cross-species comparisons using TranscriptFormer.

R file:
- hexacoral_microsynteny.R

  Analysis and visualisation of microsynteny across hexacoral genomes.
  
Reference:

Lewin TD, Sakagami T, Piñon-Gonzalez VM, Yoshioka Y, Kao LJ, Chiu YL, Sasaki K, Chen YH, Li JY, Tin KX, Lu MYJ, Miller DJ, Shikina S, Lin MF, Luo YJ (2026) Corallimorpharian genome supports the monophyly of Scleractinia and illuminates the cellular evolution of coral calcification. bioRxiv 2026.09.25.754439. https://doi.org/10.64898/2026.09.25.754439
