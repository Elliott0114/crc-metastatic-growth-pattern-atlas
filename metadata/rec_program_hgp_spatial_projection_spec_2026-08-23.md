# REC program projection into HGP spatial transcriptomics

Date: 2026-08-23  
Status: prospective analysis specification written after the Ogden REC program was frozen and before inspecting its E-MTAB-12043 projection

Correction (2026-08-24): an overlap audit identified KRT20 in both the original
epithelial marker set and the frozen REC programme. KRT20 was removed from mask
construction, and the spatial analysis and all dependent summaries were rerun.
The remaining epithelial and hepatocyte markers are disjoint from the programme.

Closed-submission addendum (2026-08-24): the five-hop region remains the primary
boundary definition. A focused three-/five-/seven-hop sensitivity analysis was
added for reader interpretation. None of the three definitions supported stable
near-interface enrichment or an HGP-specific localisation effect; the complete
patient-level results are reported in Supplementary Table 2.

## Purpose

Test whether the fixed full-transcriptome program of the Ogden REC malignant state is visible in an independent, HGP-labelled CRLM spatial cohort and whether it is concentrated on the tumour side of the tumour–liver boundary. This analysis is a spatial projection of a frozen program, not a new signature-discovery step.

## Fixed program

- Source: the 50 highest-ranked genes in `analysis_results/ogden_anchor_program/selected_anchor_program_top50.tsv`.
- Direction: all genes were higher in REC than in other epithelial states in the paired Ogden analysis.
- No gene will be added, removed or reweighted according to E-MTAB-12043 results.
- Genes absent from all six spatial expression matrices will be reported and omitted uniformly. At least 40 of 50 genes must be available for the projection to be considered adequately covered.

## Spatial cohort and unit of analysis

- Dataset: E-MTAB-12043, six CRLM Visium sections from six patients, comprising three rHGP and three dHGP lesions.
- Statistical unit: patient/section. Spots are aggregated within patient and region before between-HGP comparisons.
- Primary spot QC: at least 200 detected genes. Thresholds of 100 and 500 detected genes are descriptive sensitivity analyses.

## Region definition

The sample-internal epithelial-versus-hepatocyte contrast defines tumour-side and liver-side spots without using REC-programme genes. The corrected epithelial set is EPCAM, KRT8, KRT18, KRT19 and CEACAM5; the hepatocyte set is ALB, APOA1, APOA2, ASGR1, CPS1 and TTR. Graph distance on the Visium hexagonal lattice is measured from liver-side spots.

- Near-interface tumour: tumour-side spots one to five graph hops from liver-side spots (approximately 0–500 µm).
- Deep tumour: tumour-side spots at least six graph hops from liver-side spots (>500 µm).
- Tumour side: all tumour-side spots, used as a secondary regional comparison.

These are computationally defined spatial bands. They are not manual pathology annotations and are not described as direct cell contact.

## Program score

For each patient-region pseudobulk, calculate gene-level log2 counts per million. Standardize each program gene across the patient-region values entering a given comparison, then average the standardized values with equal gene weights. The fixed direction means a larger score represents greater REC-program expression.

For interpretability, also report pooled program log2 CPM obtained by summing counts over all available program genes. It is secondary because abundant genes receive more weight.

## Primary questions

1. Is the REC-program score higher in rHGP than dHGP in the near-interface tumour band?
2. Is the within-patient near-interface-minus-deep-tumour program difference larger in rHGP than dHGP?

Report all six patient values, the rHGP-minus-dHGP mean difference, Hedges' g, and the exact 3-versus-3 label-permutation P value. Because only 20 HGP allocations exist, effect size and patient-level consistency take priority over an arbitrary significance threshold.

## Simple sensitivity checks

- Repeat both questions at the 100- and 500-detected-gene QC thresholds.
- Repeat the near-interface and localisation summaries with three, five and seven
  graph hops assigned to the boundary band; five hops remains the primary analysis.
- Repeat the primary analysis after retaining the more epithelial half of tumour-side spots within each sample. This is a composition check, not a replacement main analysis.
- Report leave-one-patient-out direction for the two primary effects.

## Interpretation rule

- Strong spatial support: near-interface rHGP effect and localization interaction are both positive at the primary threshold, keep their direction at both QC thresholds, and at least five of six leave-one-patient-out estimates remain positive for each endpoint.
- Regional support: the near-interface rHGP effect is stable, but preferential near-versus-deep localization is not.
- No independent HGP spatial support: the near-interface effect is not directionally stable.

The result may support spatial compatibility of the REC program with rHGP. It cannot prove that Visium spots are pure malignant cells, that REC and NC3 are identical lineages, or that AP-1/NF-κB activity causes replacement growth.
