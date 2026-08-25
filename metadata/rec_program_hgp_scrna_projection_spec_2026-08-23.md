# REC program projection into HGP single-cell RNA sequencing

Date: 2026-08-23  
Status: prospective analysis specification written before inspecting the REC-program projection in E-MTAB-12022

## Purpose

Determine whether the REC program observed in E-MTAB-12043 Visium data is carried by reconstructed CRC/epithelial cells rather than being explained only by mixed spatial spots. E-MTAB-12022 and E-MTAB-12043 arise from the same source study and share some patients; they are complementary modalities, not independent patient cohorts.

## Cohort and statistical unit

- Input: `analysis_results/e_mtab_12022_primary_reconstruction.h5ad`, reconstructed from public tumour FASTQs using the project's fixed cell-calling and broad-compartment annotation pipeline.
- Primary compartment: cells annotated as `CRC/epithelial` by the reference-marker method.
- Patient eligibility: at least 20 reconstructed CRC/epithelial cells.
- Statistical unit: patient. Counts are summed within patient before HGP comparison; cells are never treated as independent replicates.
- Expected eligible set from the already-audited reconstruction: five patients (three rHGP and two dHGP). PT52 has fewer than 20 epithelial cells and is excluded by the existing threshold.

## Frozen program and score

- Full program: the 50 REC-up genes frozen in `analysis_results/ogden_anchor_program/selected_anchor_program_top50.tsv`.
- Use every gene present in the reconstructed matrix; require at least 40 genes.
- For each patient, calculate gene-level log2 CPM from epithelial pseudobulk counts, standardize each gene across eligible patients, and average genes equally.
- Also report pooled-program log2 CPM after summing program counts. This is a direct abundance-weighted summary and is secondary to the equal-gene score.

## Mechanistic decomposition

Four small, source-defined subsets are projected descriptively. They were fixed from the Ogden signature/regulon enrichment before this HGP single-cell projection:

- Hypoxia MP6 overlap: ADM, ANGPTL4, ANKRD37, EGLN3, LDHA, NDRG1, NDUFA4L2, SLC16A3.
- Regenerative overlap: ADM, APOL1, DUOXA2, ISG15, OAS1, SLC16A3.
- Experimental AP-1-target overlap: CYP3A5, FHL2, GSN, KRT80, PLAUR.
- NF-κB-regulon overlap: ABCG1, BIRC3, DUSP5, MXD1, PLAUR, SDCBP2.

These subset scores are interpretability readouts, not independent validation endpoints or evidence of transcription-factor causality.

## Inference and robustness

- Report every patient value, rHGP-minus-dHGP mean difference, Hedges' g, stratified patient bootstrap 95% interval, and complete exact label allocation preserving three rHGP and two dHGP labels.
- Primary endpoint: equal-gene full REC-program score.
- Secondary endpoint: pooled full-program log2 CPM.
- Cell-recovery sensitivity: sample 38 epithelial cells without replacement from each eligible patient in 1,000 iterations (seed 42), rebuild patient pseudobulks, and report the fraction of iterations with a positive rHGP-minus-dHGP effect.
- QC-threshold sensitivity: rebuild epithelial patient pseudobulks across minimum-feature thresholds of 150, 200, 300 and 500 and maximum mitochondrial fractions of 20%, 30%, 50% and 70%. Report cohort retention and direction; settings that change the eligible patient set are descriptive and do not replace the source-study reconstruction rule.
- Leave-one-patient-out estimates are directional influence checks only.

## Interpretation

A stable positive patient-level program effect supports an epithelial source for the HGP-associated regional signal. It does not distinguish increased frequency of a discrete REC-like state from higher expression along a continuous program, and it does not prove lineage identity with Ogden REC cells.
