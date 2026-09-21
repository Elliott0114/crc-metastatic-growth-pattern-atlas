# Fixed REC-program summary in the Latacz tumour–liver interface cohort

Date: 2026-08-23  
Status: specified before inspecting REC-program effects in the Latacz workbook

## Question

Do the genes in the frozen REC program show a coherent replacement-HGP direction in a published tumour–liver interface RNA-sequencing cohort?

## Inputs

- Published Supplementary Table 1 from Latacz et al.: `data_sources/Latacz_2024_targeted_HGP_interface/SupplTablesLataczforupload.xlsx`.
- Frozen 50-gene REC program: `analysis_results/ogden_anchor_program/selected_anchor_program_top50.tsv`.
- Descriptive 15-gene cross-modal HGP core: `analysis_results/rec_program_cross_modal_concordance/top15_interpretable_core.tsv`.

The source table contains gene-level log2 replacement/desmoplastic effects from an uncorrected model and models adjusted for systemic treatment, consensus molecular subtype, and tumour–stroma ratio. It represents 57 pure-HGP interface specimens from 51 patients (33 replacement and 24 desmoplastic specimens).

## Fixed summaries

For each model and each fixed gene set, report:

- number of measurable genes;
- mean and median log2 replacement/desmoplastic effect;
- number and fraction of genes with a positive effect;
- interquartile range of gene effects.

Also report how many program genes are positive in all four models. The uncorrected model is the primary descriptive view; adjusted models are direction checks.

## Claim boundary

This source provides published gene-level summary statistics, not patient-level expression. Correlated genes cannot be treated as independent replicates, so no gene-level sign test or combined program P value is calculated. The cohort is corroborating study-level evidence; GSE151165 supplies the independent patient-level bulk test.
