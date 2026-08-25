# Fixed REC-program projection into GSE151165 bulk CRLM

Date: 2026-08-23  
Status: specified before inspecting REC-program scores in GSE151165

## Question

Is the frozen Ogden-derived REC epithelial program higher in replacement- than desmoplastic-HGP tumour samples in an independent bulk CRLM cohort?

## Inputs

- Raw counts: `data_sources/GSE151165/GSE151165_RNA_seq.raw_read_count.xlsx`.
- Sample contract: 9 desmoplastic tumour samples and 6 replacement tumour samples, with unique public KR identifiers.
- Frozen program: `analysis_results/ogden_anchor_program/selected_anchor_program_top50.tsv`.

Normal-adjacent samples are not used in the primary comparison. The 15 tumour samples are treated as independent patient-level observations because the public columns have distinct patient identifiers.

## Primary summary

1. Apply TMM library-size normalisation.
2. Convert each available program gene to log2 CPM.
3. Standardise each gene across the 15 tumour samples and average the gene z scores within each sample (`rec_equal_gene`).
4. Estimate the mean replacement-minus-desmoplastic difference, Hedges' g, a group-stratified patient bootstrap 95% interval (10,000 replicates; seed 42), and an exact label-permutation P value over all 5,005 allocations of six replacement labels.

## Transparent companion summaries

- Aggregate raw counts for the matched REC genes, divided by the TMM effective library size and expressed as pooled-program log2 CPM.
- The pre-existing 15-gene cross-modal HGP core and four previously defined mechanistic subsets, calculated by the same equal-gene rule.
- Gene-level replacement-minus-desmoplastic mean differences and patient-direction summaries.
- Leave-one-patient-out direction of the primary REC-program effect.

Only the frozen 50-gene program is the primary test. The compact core and mechanistic subsets are explanatory summaries.

## Claim boundary

Bulk tissue cannot identify the cellular source of a program and can contain tumour, stromal and liver-derived RNA. Epithelial single-cell pseudobulks provide the source check. GSE151165 is therefore an independent regional HGP test, not proof of a tumour-cell-autonomous mechanism.
