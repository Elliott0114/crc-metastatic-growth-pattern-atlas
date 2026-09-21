# GATA6-loss functional competitor projection specification

**Frozen:** 2026-09-01, before downloading or inspecting GSE290753 expression values.  
**Status:** literature-triggered, post hoc mechanistic discrimination. It does not alter the frozen REC programme or primary HGP inference.

## Rationale and question

Goto et al. showed with serially selected mouse CRC organoids and genetic perturbation that GATA6 loss promotes liver metastasis, represses HNF4A, generates LGR5-negative cells and activates fetal-like and basal-like/squamous plasticity (Cell Stem Cell 2026;33:1304–1323.e7). This is a strong alternative explanation for any LGR5-low metastatic state.

The question is not whether GATA6 loss promotes metastasis—the source study established that experimentally. The question for this manuscript is:

> Does GATA6 loss reproduce the combined REC–junction–low-cycle phenotype, or does it generate a different LGR5-low plastic state?

This comparison can distinguish a general metastatic-plasticity signature from the more specific state associated with replacement HGP. GSE290753 contains no HGP annotation, so it cannot establish an rHGP mechanism.

## Data and units

- Dataset: GSE290753, mouse AKP CRC organoids, three GATA6-knockout and three control RNA-seq libraries.
- Each library is the inferential unit.
- GEO sample titles define condition; no condition is inferred from expression.
- The processed GEO series matrix will be used if it contains a complete gene-level expression matrix. If values are nonnegative integers, they will be analysed as counts with TMM–voom; otherwise nonnegative abundance values will be transformed as `log2(value + 0.5 * smallest positive value)` and analysed with limma.
- Mouse–human programme transfer uses case-insensitive exact gene-symbol identity, with full coverage reporting.

## Frozen programme family

1. `REC_PROGRAM_TOP50` — observed HGP-associated state.
2. `CORE_HRC` — metastatic high-relapse repair state.
3. `PUBLISHED_FETAL` — positive control for the source-reported fetal programme.
4. `PARTIAL_EMT` and `HALLMARK_EMT` — plasticity route.
5. `REACTOME_TIGHT_JUNCTION_INTERACTIONS` — broad epithelial-junction phenotype.
6. `JUNB_TIGHT_JUNCTION_TARGETS_5` — focused REC-associated junction subset.
7. `ACTIN_TURNOVER` — shape-remodelling component.
8. `HALLMARK_E2F_TARGETS` and `HALLMARK_G2M_CHECKPOINT` — cycling boundary.

Named genes are `GATA6`, `HNF4A`, `LGR5`, `MKI67`, `EMP1`, `KLF4`, `EPCAM`, `ZEB1`, `CDH17`, `TJP1`, `CLDN2`, `CLDN3`, `CLDN4`, `CLDN7`, `CRB3` and `F11R`.

## Frozen estimands and inference

For each programme:

1. calculate an equal-gene score by standardizing every measurable gene across the six libraries and averaging without abundance weighting;
2. report mean GATA6-KO minus control difference, a 95% library-bootstrap interval, Hedges' *g*, the number of KO samples above the control median and the exhaustive two-sided 3-versus-3 label-allocation *P* value;
3. report leave-one-gene-out direction ranges for REC, junction and the focused five-gene set;
4. run limma gene-level modelling and `camera` competitive gene-set tests using the same design;
5. control Benjamini–Hochberg FDR across this ten-programme family separately for module contrasts and `camera` results.

Named-gene effects are reported with limma log2 fold change, 95% interval, moderated *P* and family-wise BH FDR. Because there are only six libraries, effect direction and programme coherence carry more information than the coarse exact permutation resolution.

Base random seed is 42. Bootstrap resamples whole libraries within condition 10,000 times.

## Interpretation without hard gates

- Joint REC/HRC/fetal enrichment with retained or increased junction programmes and reduced cycling would resemble the full replacement-associated state.
- Fetal/EMT enrichment with junction loss would define a mechanistically distinct LGR5-low route and strengthen the specificity of the REC–junction combination.
- REC enrichment without junction or low-cycle agreement would show partial signature overlap and argue against calling GATA6 the upstream rHGP driver.
- No REC enrichment would place GATA6 loss as a competing metastatic mechanism rather than an explanation of the observed state.

The result will be framed as a functional competitor projection, not as independent HGP validation and not as proof that GATA6 is causally involved in replacement growth.
