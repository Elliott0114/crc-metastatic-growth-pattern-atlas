# GSE294385 micro-to-macrometastasis extension of the REC program

Date: 2026-08-23  
Status: prospective analysis specification written before inspecting REC-program expression in GSE294385

## Purpose

Test where the HGP-associated REC/AP-1 epithelial program lies along liver metastatic outgrowth. GSE294385 contains source pathology labels for liver micrometastasis tumour, micrometastasis stroma, macrometastasis tumour and macrometastasis stroma in 11 patients.

## Data and unit of analysis

- Expression: public filtered 10x Visium matrices for all 24 liver sections in GSE294385.
- Regions: source-provided `visium_meta_after_qc.rds` labels from the authors' code repository.
- Primary unit: patient. Raw counts are summed over all labelled spots within each patient-region and normalized after aggregation.
- Primary contrast: paired liver macrometastasis tumour minus liver micrometastasis tumour across 11 patients.

## Fixed programs

- Primary: frozen 50-gene Ogden REC-up program.
- Secondary abundance-weighted summary: pooled counts over the same 50 genes.
- Compact secondary program: the 15 highest cross-modal positive HGP genes, frozen in `analysis_results/rec_program_cross_modal_concordance/top15_interpretable_core.tsv` before GSE294385 expression is inspected.
- Mechanistic subsets: Hypoxia MP6 overlap, regenerative overlap, experimental AP-1 targets and NF-κB regulon overlap, with the same definitions as prior analyses.
- Source positive control: the published micrometastasis six-gene signature RNF40, AEN, WEE1, BCL7B, YME1L1 and COX17. It is expected to be higher in micrometastasis tumour and verifies region extraction and direction.

No GSE294385 result will be used to add, remove or reweight genes in these programs.

## Scoring and paired inference

- For each comparison, calculate patient-region gene log2 CPM, standardize each gene across the paired patient-region pseudobulks, and average genes equally.
- Report paired macro-minus-micro differences, positive-patient count, mean and median paired difference, paired standardized effect (mean difference divided by the standard deviation of paired differences), patient-pair bootstrap 95% interval, and exact sign-flip P value.
- The full equal-gene REC score is the primary endpoint. Pooled REC CPM and the compact HGP core are secondary. Mechanistic subsets are interpretability analyses.

## Robustness and compartment checks

- Spot-count balance: in 1,000 iterations (seed 42), downsample each patient's macrometastasis tumour spots to that patient's micrometastasis tumour spot count, while retaining all micrometastasis spots.
- Micro-focus size: repeat paired summaries among patients with at least 20 micrometastasis tumour spots.
- Section control: repeat the analysis using only sections containing both micro- and macrometastasis tumour labels; average section-level differences within patient. PAT5, whose micro and macro lesions occur in separate sections, is excluded from this check.
- Compartment specificity: compare tumour with matched stroma separately for micro- and macrometastases at the patient level.

## Interpretation

- A positive macro-minus-micro REC effect places the program with metastatic outgrowth/expansion.
- A negative effect places it with early micrometastatic persistence.
- A null effect suggests the HGP-associated program is not primarily ordered by lesion size.

This extension does not provide HGP labels and cannot establish that micro- or macrometastases use a specific histopathological growth pattern.
