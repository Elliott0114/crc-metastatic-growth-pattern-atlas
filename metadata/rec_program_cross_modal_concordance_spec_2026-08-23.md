# Cross-modal concordance of the frozen REC program

Date: 2026-08-23  
Status: analysis specification written before inspecting gene-level HGP effects

## Purpose

Reduce the 50-gene REC program to an interpretable cross-modal result without replacing the frozen validation score. The analysis asks which REC genes show the same rHGP direction in reconstructed epithelial single-cell pseudobulks and in near-interface Visium pseudobulks.

## Inputs

- Ogden REC-versus-other epithelial effect and program rank: `analysis_results/ogden_anchor_program/selected_anchor_program_top50.tsv`.
- E-MTAB-12022 epithelial patient pseudobulks: `analysis_results/rec_program_hgp_scrna_projection/patient_gene_pseudobulk.tsv`.
- E-MTAB-12043 spatial patient-region pseudobulks: `analysis_results/rec_program_hgp_spatial_projection/patient_region_gene_pseudobulk.tsv.gz`.

The single-cell and spatial layers originate from the same source study and share some patients. Agreement is therefore cross-modal consistency, not replication in two independent cohorts.

## Fixed comparisons

- Single-cell effect: mean patient epithelial log2 CPM in rHGP minus dHGP.
- Spatial effect: mean patient log2 CPM in rHGP minus dHGP at the primary 200-detected-gene threshold, all tumour-side spots in the 0–500 µm band.
- Use the intersection of frozen program genes available in both HGP modalities.
- Report the number and fraction positive in each modality, positive in both, and opposite in the two modalities.
- Report Spearman correlation between single-cell and spatial gene effects. The gene-level correlation is descriptive because genes are correlated and were selected in the Ogden source analysis.

## Interpretability subsets

Summarize the same four pre-defined subsets used in the HGP single-cell projection: Hypoxia MP6 overlap, regenerative overlap, experimental AP-1 targets, and NF-κB regulon overlap. For each subset, report median effects and the fraction of genes positive in both HGP modalities. No gene-level or subset-level P values are used for mechanism selection.

## Concordant core

A gene is labelled `cross_modal_positive` when both HGP effects are greater than zero. This label is a descriptive core for interpretation and visualization. It must not replace the original 50-gene program in validation analyses, and no claim is based on any single gene.
