# Focused evidence bridge, 2026-09-16

Status: exploratory extension specified after the junction-context and measurement
results were known, before running the analyses below. This is not preregistration.

## Questions and decision rules

1. Do paired raw-count models support an overall junction decrease from micro- to
   macrometastatic tumour regions, after accounting for library size and patient?
2. In the six HGP-labelled E-MTAB-12043 spatial sections, do the fixed junction
   components covary with HRC in epithelial-enriched spots, differ by HGP, and
   localize preferentially to the inferred tumour–liver boundary?
3. Are the same junction members responsible for the bulk HGP, epithelial HRC,
   regional and spatial observations? Report disagreements, not a selected core.

Do not equate association with regulation, boundary proximity with cell contact,
or different assays from overlapping patients with independent replication.
Failure to support localization or HGP specificity leaves that connection open.

## Fixed genes

Use `rec_junction_context_common.definitions()`: Claudin 21, Polarity 8,
Junction 30 (also including F11R), and HRC 98 after removal of junction genes.
For E-MTAB-12043 only, exclude CEACAM5 from HRC because the fixed epithelial mask
uses that gene. No junction gene overlaps the mask. Retain every measured member
in coverage tables, including low-abundance members excluded from count tests.
No outcome-driven gene selection or new signature.

## Raw-count models

GSE294385: 22 patient-region profiles, 11 paired patients; design `~patient+macro`.
Use the common assayed feature intersection already audited, and full original
library sizes (not sums of the gene intersection). GSE151165: 15 tumour profiles,
6 rHGP and 9 dHGP; design `~rHGP`; the count table supplies the library totals.
The two studies have different endpoints and are not pooled.

Primary: edgeR negative-binomial quasi-likelihood fits, robust dispersion,
`prior.count=0` for reported coefficients, with library offsets and TMM factors.
Sensitivity: raw library offsets without TMM, keeping the primary gene eligibility
fixed. Eligibility is CPM >=1 in at least half the profiles (rounded up), plus
>=15 total counts, with no outcome labels used. TMM is estimated from eligible
common genes while retaining full library sizes. Report factor ranges and design
rank. These models do not turn regional expression into per-cell expression.

For each fixed set, edgeR's `fry.DGEList` tests directional and mixed changes
using negative-binomial transformed residuals and accounts for gene correlation;
these are approximate set tests, not QL tests of a single programme coefficient.
BH correction across the four predefined sets for each contrast/normalization;
genome-wide BH and fixed-junction-family BH for gene-level QL tests. Median gene
logFC and number of positive coefficients are descriptive, not inferential set
effect estimates. No cross-gene independence assumption. Do not select between
normalizations based on the smaller P value.

## HGP spatial bridge

Use local reconstructed E-MTAB-12043 Space Ranger outputs, six patients, one
section each, three rHGP and three dHGP. HGP comes from source SDRF/manifest.
No source spot-level pathology or malignant-cell annotation is available here.
The same sample-internal epithelial/hepatocyte tertile contrast and hex-grid graph
distance as the frozen REC projection define epithelial-enriched tumour-side
spots and near (1–5 hops) versus deep (>=6 hops) bands. Never call these masks
pathology labels or pure tumour cells. In-tissue spots with >=200 detected genes;
all tumour-side spots primary, upper epithelial half sensitivity. Additional
3-/7-hop boundaries are fixed sensitivity checks, not candidate cutoffs.

Spot scores: mean log1p(CP10k) of common assayed members, no gene z-weighting.
HRC excludes junction and mask genes. Within each patient, partial Spearman HRC
versus Claudin/Polarity controlling epithelial and liver proxy scores, log library
size and log detected genes. Report unadjusted and adjusted correlations, then
equal-patient Fisher-z summaries with patient bootstrap intervals and exhaustive
sign flips. These within-section descriptive associations do not remove spatial
autocorrelation or establish cause; no spot-level inferential P values.

Regions: sum raw counts and full libraries within each patient-band; >=20 spots
required for an evaluable region. Report equal-patient near-minus-deep programme
effects on library-scaled logCPM mean expression and mean gene-z expression.
HGP contrasts and HGP-by-localization differences use exhaustive 3-versus-3
patient allocations (20 assignments; minimum two-sided P is 0.10), bootstrap
intervals and leave-one-patient-out ranges. No multivariable between-patient
regression at n=6. Correction across the four sets per endpoint; all outcomes
and all sensitivity results are retained.

For the primary five-hop/all-tumour region counts only, also fit edgeR models:
`~rHGP` for tumour-side totals and near-band counts, and
`~patient+near+near:rHGP` for the paired HGP-by-localization interaction.
These exploratory model P values may be smaller than exact patient permutations;
report both honestly and do not treat modelling assumptions as new replication.

## Verification and deliverables

Check unique identifiers, count integrality, barcode-coordinate alignment,
source HGP labels, common feature coverage, disjoint masks/targets, full-library
preservation, and design rank. Independently reproduce key spatial summaries
in R and verify primary count fits against a separately constructed design.
Save source hashes, scripts, coverage, patient values, complete result tables,
figures, Chinese evidence-chain report and candidate manuscript paragraphs.
Leave the submitted manuscript intact.
