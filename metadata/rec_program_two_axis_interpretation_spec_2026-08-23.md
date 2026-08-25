# REC-program two-axis interpretation

Date: 2026-08-23  
Status: descriptive integration specified after completion of the HGP and micro-to-macro projections

## Purpose

Place each measurable gene from the frozen 50-gene REC program on two directly interpretable biological axes:

1. **HGP axis:** whether the gene is higher in rHGP in both epithelial single-cell pseudobulks and near-interface FFPE spatial pseudobulks.
2. **Outgrowth axis:** whether the gene is higher in macrometastatic than micrometastatic tumour regions in at least 8 of 11 paired patients.

This is a descriptive synthesis, not a new validation test. It is intended to separate an HGP-associated epithelial component from a lesion-outgrowth component without fitting a latent or black-box model.

## Inputs

- HGP gene effects: `analysis_results/rec_program_cross_modal_concordance/gene_level_concordance.tsv`.
- Paired micro-to-macro effects: `analysis_results/gse294385_rec_program_extension/gene_level_paired_effects.tsv`.
- Frozen program definition: `analysis_results/ogden_anchor_program/selected_anchor_program_top50.tsv`.

The analysis uses only genes measurable in both HGP modalities and the GSE294385 liver series.

## Fixed descriptive rules

- `hgp_supported`: rHGP-minus-dHGP effect is positive in both HGP modalities.
- `outgrowth_supported`: mean macro-minus-micro effect is positive and at least 8 of 11 patients have a positive paired difference.
- `shared_hgp_outgrowth`: both labels are true.
- `hgp_biased`: HGP support only.
- `outgrowth_biased`: outgrowth support only.
- `unresolved_or_opposed`: neither rule is met.

The 8-of-11 threshold denotes a clear patient majority and is not treated as a hypothesis-test cutoff. Exact paired P values and false-discovery rates are carried forward only as descriptive annotations.

## Mechanistic summaries

For the four previously defined subsets (experimental AP-1 targets, hypoxia MP6, regenerative and NF-kB regulon overlap), report category counts and median effects on both axes. No enrichment P values are calculated for the four-way classification because the gene universe is small and source-selected.

## Claim boundary

The two-axis categories are an explanatory visualization and gene-prioritisation device. They do not imply temporal lineage, causality, or independence between the datasets. GSE294385 has no HGP label, while the HGP single-cell and spatial modalities partially share patients.
