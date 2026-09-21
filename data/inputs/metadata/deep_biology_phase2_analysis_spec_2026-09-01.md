# Deep-biology Phase 2 analysis specification

Date: 2026-09-01  
Status: timestamped before inspecting the new Phase 2 endpoint outputs  
Scope: mechanistic specificity analyses for the frozen REC programme  

## Scientific question

Do replacement-pattern colorectal liver metastases enrich a regenerative-like
malignant epithelial state that is distinguishable from generic hypoxia,
unfolded-protein stress, proliferation, lesion outgrowth and the recently
reported fatty-acid/MYC/proline/collagen mechanism?

The analysis is observational and discriminates among biological explanations.
It does not test whether the REC state causes replacement growth or whether
damaged hepatocytes induce REC.

## Fixed objects

- The existing 50-gene REC programme remains unchanged.
- The NC3-to-Ogden state mapping remains unchanged.
- New external gene sets are frozen from source publications or MSigDB
  2025.1.Hs before the corresponding Phase 2 scores are inspected.
- GSE151165 HGP labels and GSE294385 micro/macro labels are never used to select
  genes.
- Patient or pseudonymous donor is the inferential unit. Cells and spots are
  measurement units only.

## Operational definitions

| Construct | Operational definition | Role |
|---|---|---|
| Replacement growth | Source HGP label `rHGP`/`RHGP`; desmoplastic/encapsulated comparator `dHGP`/`EHGP` | Exposure |
| REC state | Source Ogden epithelial subtype `REC` and the already frozen top-50 paired REC-versus-other programme | Primary malignant-state construct |
| Generic stress | Ogden `Hypoxia`, `UPR` and `iREC` states; frozen Hallmark/author hypoxia, UPR and interferon-related sets | Rival explanation |
| Proliferation/outgrowth | Ogden TA/intermediate states, Hallmark E2F/G2M/MYC programmes, and paired GSE294385 macro-minus-micro tumour contrast | Rival explanation |
| Metabolic co-option | Peng-Winkler palmitate perturbation programme, Hallmark fatty-acid metabolism/glycolysis/WNT, GO proline and Reactome collagen programmes | Rival explanation |
| Regenerative plasticity | Source `coreHRC`, `EpiHR`, fetal, YAP, HRC KRT20+/-, Nusse regenerative and selected GO repair/regeneration programmes | Biological interpretation |
| Liver-parenchymal integration | Fraction of neoplastic cells with a liver epithelial cell among 10 nearest segmented neighbours | Architectural endpoint |
| State-specific liver adjacency | NC3-minus-other-neoplastic difference in DHC or HC1 neighbour outcomes, evaluated at patient level and against within-ROI state-label exchange | Specificity endpoint |

## Analysis A: direct epithelial-state contrasts

1. Reconstruct patient-state pseudobulk counts from the Ogden epithelial cells.
2. A direct REC-versus-comparator contrast is eligible when at least five
   patients contribute at least 20 cells to both states.
3. Fit paired limma-voom models (`~ patient + state`) separately for REC versus
   Hypoxia, UPR, iREC, TA1, Colonocyte, Stem NOTUM, Intermediate and Goblet when
   eligible.
4. Rank the full tested transcriptome by the moderated t statistic.
5. Run targeted preranked GSEA against the frozen programme registry. Hallmark
   and selected Reactome/GO results are adjusted within their declared family.
6. Report NES, FDR, leading-edge genes, paired patient consistency and the
   direction of the frozen REC programme score. Direct REC-versus-Hypoxia,
   REC-versus-UPR and REC-versus-iREC contrasts are the principal specificity
   comparisons; the remaining contrasts define differentiation context.

## Analysis B: external rival-program controls

1. In GSE151165 tumour samples, score the frozen REC programme and every frozen
   rival programme from gene-wise standardized log-CPM values.
2. Estimate the unadjusted rHGP-minus-dHGP REC effect and, one rival at a time,
   the rHGP coefficient from `REC score ~ HGP + rival score`.
3. Use exhaustive HGP-label permutation for P values and patient bootstrap for
   95% confidence intervals. Apply BH correction across rival models.
4. Repeat using the REC genes remaining after removal of overlap with each rival
   set. Record retained gene count and the change from the original effect.
5. Repeat the one-rival-at-a-time analysis on paired macro-minus-micro
   differences in GSE294385 using exhaustive sign flips. This is a lesion-
   outgrowth comparator, not HGP validation.
6. Do not fit an omnibus many-covariate model with 15 or 11 patients.

## Analysis C: spatial structure-state separation

1. Recompute 10-nearest-neighbour outcomes for every neoplastic ISS cell.
2. For NC3 and the pooled other neoplastic states, calculate DHC and HC1
   neighbour presence and neighbour-slot fraction within each ROI.
3. Aggregate ROIs equally within patient-HGP units.
4. Shuffle neoplastic state labels within ROI, region and front-status strata,
   preserving state counts and all target-cell positions. Use 1,000 iterations,
   seed 42.
5. Report patient-level NC3-minus-other contrasts, their within-patient spatial-
   null residuals, exact sign-flip intervals/tests and pure-HGP group summaries.
6. Preserve the two mixed-HGP patients as paired descriptive counterexamples.

## Evidence interpretation

No arbitrary hard pass/fail threshold is used. Results update the claim as
follows:

- Persistent REC effects after strong rival controls support a distinct
  malignant-state layer in rHGP.
- Partial attenuation supports biological overlap and requires naming the shared
  mechanism.
- Near-complete attenuation supports reframing REC as a transcriptional readout
  of the rival process.
- Broad liver adjacency without NC3-specific DHC adjacency supports a two-layer
  model: architecture-wide liver integration plus lesion-level REC enrichment,
  without a DHC-to-REC coupling claim.
- Stable NC3-specific DHC adjacency would add spatial coherence but remains
  associational.
- Failure of boundary localization is reported as absence of stable confinement,
  not evidence for a distributed field.

All positive and negative Phase 2 results enter the final evidence matrix. Gene
sets, patients or endpoints are not replaced after results are observed.

