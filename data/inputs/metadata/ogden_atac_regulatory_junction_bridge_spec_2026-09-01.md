# Ogden ATAC-informed regulatory-to-junction bridge specification

Date: 2026-09-01  
Status: frozen before inspecting the new regulon-to-junction score outcomes  
Scope: source-constrained mechanistic prioritisation in the Ogden multiome cohort

## Scientific question and value

The manuscript currently shows that repair/plasticity and epithelial-junction
transcription coexist in REC cells. That observation is biologically useful but
does not identify a candidate regulatory system capable of coordinating the two
features. Ogden et al. inferred transcription-factor regulons by integrating RNA
expression with chromatin accessibility and identified AP-1 and NF-κB as major
regulators of regenerative malignant states.

This analysis asks: **does an ATAC-informed regulator, particularly the AP-1
axis, connect REC identity to junctional transcription more consistently than
NF-κB, YAP/TEAD, WNT/stemness or differentiated-intestinal regulator controls?**
The intended output is a short, experimentally actionable regulator-to-target
bridge, not another catalogue of enriched pathways.

Regulon membership is a source-predicted regulatory link. Expression covariance
can support prioritisation but cannot prove TF binding, regulator activity,
assembled junctions or causality.

## Fixed inputs

- Ogden liver-metastasis RNA counts and published epithelial-state labels.
- Source Supplementary Table `mmc8.xlsx`, including 18 SCENIC+ regulons and the
  published combined AP-1 and NF-κB target lists.
- Frozen `REACTOME_TIGHT_JUNCTION_INTERACTIONS` genes used in the existing
  within-cell analysis.
- Existing frozen 50-gene REC programme and recurrent REC-versus-Hypoxia/UPR/
  iREC differential-expression results; neither will be reselected.
- Primary state: REC. Context states: Hypoxia, UPR and iREC.
- Patient is the inferential unit; cells are measurement units.
- Random seed 42; patient bootstrap 10,000 iterations.

## Regulator families

All 18 source SCENIC+ regulons will be retained in the complete output. The
reader-facing candidates and controls are fixed as follows:

- **Primary AP-1 bridge:** combined AP-1 regulon, JUNB regulon and FOSL2 regulon.
  JUNB and FOSL2 are prioritised because the already-completed, pre-existing
  frozen-REC-program ORA identified them as the two strongest individual
  SCENIC+ overlaps; this prior result is recorded rather than rediscovered.
- **Inflammatory comparator:** combined NF-κB regulon and RELB regulon.
- **Mechanotransduction comparators:** TEAD1 and TEAD4 regulons.
- **Lineage/WNT comparators:** ASCL2, LEF1, CDX2 and HNF4A regulons.
- The remaining source regulons are reported as exploratory context.

## Analysis A: regulon-to-junction target map

1. Export every source regulon without changing membership.
2. Use the genes tested in the existing Ogden REC-versus-other epithelial
   limma-voom analysis as the ORA background; do not use the whole genome.
3. Test over-representation of the frozen tight-junction set within each
   regulon by one-sided Fisher exact test, with BH adjustment across all source
   regulons.
4. Separately enumerate genes that satisfy all three conditions: junction-set
   member, source regulon target and recurrently REC-up versus Hypoxia, UPR and
   iREC. These genes are the concrete experimental target bridge.
5. Report odds ratios, overlap counts, complete gene lists and the background
   size. An overlap of one or two genes is not treated as a mechanistic module
   regardless of its P value.

## Analysis B: patient-paired state activity

1. Reuse existing patient-state pseudobulk counts and eligibility rules.
2. Transform counts to log1p CPM and standardise each tested gene across
   pseudobulks before averaging the regulon targets, so large regulons do not
   dominate by absolute expression.
3. For each eligible patient and comparator, calculate REC-minus-comparator
   regulon-score differences.
4. Report the equal-patient mean, median, positive-patient count, patient
   bootstrap 95% CI and exhaustive paired sign-flip P value.
5. BH adjustment is applied across all tested regulons within each comparator;
   candidate-family summaries retain raw effects and uncertainty.

## Analysis C: within-cell regulon-to-junction covariance

1. For each regulon-junction pair, remove every shared target gene from both
   sides before scoring. This prevents a gene from correlating with itself.
2. Require at least five detected genes on each side.
3. Construct disjoint expression-bin-matched control sets from epithelial genes
   outside all source target sets. The primary score is mean target expression
   minus mean matched-control expression on log1p counts per 10,000.
4. Within every patient-state group, residualise both scores on log1p library
   size and log1p detected-gene count, then calculate Spearman correlation.
5. Primary summary: REC groups with at least 20 cells. Sensitivities use 10 and
   30 cells, unadjusted matched scores and gene-balanced z scores.
6. Combine patient correlations on the Fisher-z scale, report patient-equal
   correlation, median, range, sign consistency, bootstrap 95% CI and exhaustive
   sign-flip P value.
7. Compare REC correlations descriptively and by paired patient differences
   with Hypoxia, UPR and iREC where the same patients are evaluable.
8. BH adjustment is applied across all tested regulons within state, method and
   cell-threshold families.

## Interpretation without a hard gate

The conclusion will follow the whole evidence pattern rather than a single
threshold:

- An AP-1-centred bridge would be strongest if AP-1/JUNB/FOSL2 regulons contain
  multiple recurrent REC-up junction genes, rise in REC versus the three stress/
  inflammatory comparators, and covary with disjoint junction scores in most
  REC patients.
- A signal shared equally by AP-1, NF-κB and TEAD controls would support broad
  epithelial activation rather than one named regulatory axis.
- State elevation without within-cell covariance would indicate that the
  regulator marks REC but does not explain its junction component.
- Within-cell covariance without state elevation would indicate a continuous
  epithelial response shared across states.
- Sparse or inconsistent evidence will be retained as a candidate for
  perturbation, not promoted to the manuscript's causal conclusion.

Positive findings may justify naming an ATAC-informed candidate regulator in
the mechanistic model. Negative findings will keep the paper centred on the
observable repair-with-junction phenotype rather than an unsupported upstream
driver.
