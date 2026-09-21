# External projection of the REC regulatory-junction model

Date: 2026-09-01  
Status: frozen before inspecting any new external regulatory-module outcome  
Random seed: 42

## Biological question

The current manuscript separates three candidate functions within the REC-like
state enriched in replacement-pattern colorectal liver metastasis (rHGP CRLM):

1. a JUNB target bridge to epithelial-junction genes;
2. a FOSL2-associated motility/plasticity arm; and
3. a RELB/non-canonical NF-kB-associated survival arm.

The source Ogden multiome cohort nominates these components but does not provide
independent HGP validation of them. The primary question is therefore:

**Do source-frozen JUNB, FOSL2 and RELB/NF-kB programmes show an rHGP direction
outside the Ogden cohort after removing every frozen REC-program and tight-
junction gene from the regulatory scores?**

The disjoint scoring is essential. It asks whether the candidate regulatory
arms carry information beyond the observable REC and junction phenotypes rather
than recovering the same genes under a new label.

## Evidence hierarchy and non-independence

- **Discovery and mechanism nomination:** Ogden RNA-ATAC multiome data.
- **Primary independent patient-level test:** GSE151165 bulk CRLM tumours (9
  dHGP, 6 rHGP).
- **Independent interface-level corroboration:** Latacz et al. published
  gene-level rHGP/dHGP coefficients from 57 pure-HGP interface specimens in 51
  patients. These correlated gene summaries do not support a patient-level
  programme P value.
- **Epithelial-compartment and spatial localization:** E-MTAB-12022 scRNA-seq
  and E-MTAB-12043 FFPE Visium. Four patients (PT44, PT54, PT55 and PT61) occur
  in both modalities. They are one multimodal source, not two independent
  replications, and will never be meta-analysed as independent cohorts.

All data are public and de-identified. Patient, not cell, spot or gene, is the
inferential unit wherever patient-level measurements are available.

## Frozen gene sets

### Observable anchors

- `FROZEN_REC_TOP50`: the existing 50-gene REC programme.
- `REACTOME_TIGHT_JUNCTION_INTERACTIONS`: the existing 30-gene Reactome set.
- `JUNB_TIGHT_JUNCTION_TARGETS_5`: CLDN3, CLDN4, CLDN7, CRB3 and F11R, defined
  by the source JUNB-regulon/tight-junction overlap.

### Candidate regulatory arms

- JUNB regulon.
- FOSL2 regulon.
- RELB regulon.
- source combined NF-kB target set.
- source combined AP-1 target set.

### Comparator programmes

- HNF4A and CDX2 regulons: differentiated-intestinal lineage controls that
  overlap junction targets in the source data but move oppositely to REC.
- ASCL2 regulon: WNT/stemness control.
- TEAD1 regulon: mechanotransduction control.

For every candidate and comparator regulon, all genes in either
`FROZEN_REC_TOP50` or `REACTOME_TIGHT_JUNCTION_INTERACTIONS` are removed before
external scoring. Gene membership is otherwise unchanged. The freeze script
will record source size, removed genes, retained size and a SHA-256 digest for
each set.

### Literature-anchored named genes

JUNB, FOSL2, RELB, ITPR3 and NFKB2 will be reported as targeted upstream genes;
CLDN3, CLDN4, CLDN7, CRB3 and F11R will be reported as the direct junction
targets. ITPR3 is included because an in-vivo functional screen placed
ITPR3/calcium upstream of RELB in CRC liver colonization. Named-gene results are
descriptive secondary evidence and will not replace the regulon analyses.

## Common scoring and statistics

1. Counts are normalized within each dataset using its existing audited
   approach: TMM log2-CPM for GSE151165 and library-size log expression for the
   reconstructed E-MTAB data.
2. Each measurable gene is standardized across the patient-level samples being
   compared; a programme score is the equal-gene mean of those z scores.
3. Coverage, zero-variance exclusions and exact gene membership are reported.
   No outcome-based gene deletion or reweighting is allowed.
4. The effect is mean score in rHGP minus mean score in dHGP, accompanied by
   Hedges' g, a patient-bootstrap 95% interval (10,000 resamples within HGP),
   and an exhaustive two-sided label-permutation P value.
5. BH adjustment is applied across all frozen module tests within each dataset
   and region. The ten named genes form a separate BH family.
6. Effects, intervals and patient directions are interpreted continuously; no
   self-created biological pass/fail threshold will be used.

## Analysis A: GSE151165 independent bulk projection

- Tumour samples only; the existing fixed sample order and HGP labels are
  reused.
- The unadjusted HGP difference is primary.
- One-covariate sensitivity models separately adjust for the existing frozen
  pan-epithelial and hepatocyte-context scores. Exact label permutation is
  repeated for each adjusted HGP coefficient. These surrogates do not constitute
  direct tumour-purity measurement.
- Leave-one-patient-out effects are reported for every module.
- The most informative pattern would be a positive disjoint FOSL2 or RELB/NF-kB
  effect that is not mirrored uniformly by lineage and TEAD controls. A broad
  positive shift across all controls would instead indicate nonspecific tissue
  composition or epithelial activation.

## Analysis B: Latacz interface-coefficient projection

- For the uncorrected, treatment-adjusted, CMS-adjusted and tumour-stroma-ratio-
  adjusted source models, report measurable genes, mean and median published
  log2(rHGP/dHGP) coefficient, interquartile range, and positive-gene fraction.
- Report how many genes retain a positive direction in all four source models.
- Do not calculate a sign test, enrichment P value or combined programme P
  value from genes, because genes are correlated and patient-level expression is
  unavailable.

## Analysis C: E-MTAB epithelial pseudobulk

- Use reconstructed cells annotated `CRC/epithelial` under the fixed primary
  cell-calling rules.
- Include a patient only when at least 20 epithelial cells are available. This
  retains PT44, PT54, PT55, PT59 and PT61 and excludes PT52 (14 epithelial
  cells).
- Aggregate raw counts by patient and score each patient equally.
- Because the eligible patients have strongly unequal cell counts, repeat the
  scoring after sampling 38 epithelial cells per patient without replacement in
  500 fixed-seed iterations. Report the median, 2.5th-97.5th percentile and
  positive fraction of the rHGP-dHGP effect. This is a measurement-depth
  sensitivity, not a confidence interval and not a second inferential test.
- With only ten possible HGP allocations, P values are resolution-limited;
  direction and effect uncertainty are emphasized.

## Analysis D: E-MTAB FFPE spatial projection and co-state analysis

### Patient-level regional projection

- Reuse the fixed QC rule of at least 200 detected genes per in-tissue spot and
  the existing epithelial-versus-hepatocyte computational masks.
- Primary region: all `tumour_side` spots, matching the manuscript's broad
  malignant-compartment model.
- Sensitivity region: computational `interface` spots within two Visium hops of
  liver-side spots.
- Aggregate raw counts within patient and region, then apply the common
  patient-level module scoring and exact HGP comparison.

### Within-region disjoint covariance

- In tumour-side spots, and separately in interface spots, calculate gene-
  balanced programme scores within each patient for the five candidate
  regulatory arms and four comparator programmes.
- Pair every regulatory programme with the full tight-junction score; also pair
  JUNB with the five-gene direct-target score. Regulatory and junction genes are
  disjoint by construction.
- Residualize both scores within each patient on log1p UMI count, log1p detected
  genes and the epithelial-minus-hepatocyte marker score, then calculate
  Spearman correlation.
- Report each patient, patient-equal Fisher-z mean correlation with a 10,000-
  patient bootstrap interval, exact paired sign-flip P value for the overall
  correlation, rHGP and dHGP summaries, and an exhaustive HGP-label test of the
  correlation difference.
- Spot-level covariance is a localization measurement. It cannot establish
  tumour-cell identity for each spot, TF binding, assembled junctions or causal
  regulation.

## Interpretation boundary

External concordance would strengthen the candidate model but would not make it
causal. The most defensible mechanistic conclusion will follow the joint pattern:

- JUNB has the strongest source target-level bridge if its five direct targets
  retain an rHGP direction, even if the broader disjoint JUNB regulon does not;
- FOSL2 is strengthened if its disjoint programme is rHGP-directed in an
  independent cohort and spatially covaries with junction expression;
- RELB is strengthened if its epithelial disjoint programme is supported, which
  helps distinguish tumour-intrinsic activity from immune admixture in bulk;
- comparator behaviour determines whether the result is candidate-specific or
  part of a broad epithelial/lineage shift.

Discordant or null findings will be retained. They will narrow the manuscript to
the observed low-cycling repair-with-junction phenotype and keep unsupported
regulators out of the title, abstract and central conclusion.
