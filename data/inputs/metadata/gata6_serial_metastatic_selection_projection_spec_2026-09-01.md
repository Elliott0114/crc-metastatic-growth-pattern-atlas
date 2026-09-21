# Serial metastatic-selection projection specification

**Frozen:** 2026-09-01, before downloading or inspecting GSE290752 expression values.  
**Status:** literature-triggered, post hoc functional-context analysis. It does not alter the frozen REC programme or primary HGP inference.

## Biological question

Goto et al. repeatedly transplanted mouse AKP colorectal cancer organoids derived from either primary tumours or liver metastases. This experiment selected for liver-metastatic capacity without annotating histopathological growth pattern.

The analysis asks:

> Does serial liver-metastatic selection reproduce the combined REC–junction–low-cycle phenotype, or only parts of it?

Agreement would provide functional-context convergence but not rHGP validation. Disagreement would show that the REC–junction combination is not an automatic consequence of generic liver-metastatic selection.

## Data and inferential units

- Dataset: GSE290752, mouse AKP organoid RNA-seq.
- Main analysis: `Primary2`, `Met2`, `Primary3` and `Met3`, each with three libraries. These two rounds contain both primary-tumour-derived and liver-metastasis-derived organoids.
- Context only: three `Met1` libraries, because no `Primary1` control is available.
- Each RNA-seq library is the inferential unit. GEO sample titles and source metadata define round and origin; neither is inferred from expression.
- Processed integer counts will be analysed by TMM–voom/limma. Exact case-insensitive gene-symbol identity transfers human programmes to mouse, with full coverage reporting.

## Frozen programme family

The same ten-programme family used for the GATA6-knockout competitor projection is retained:

1. `REC_PROGRAM_TOP50`;
2. `CORE_HRC`;
3. `PUBLISHED_FETAL`;
4. `PARTIAL_EMT`;
5. `HALLMARK_EMT`;
6. `REACTOME_TIGHT_JUNCTION_INTERACTIONS`;
7. `JUNB_TIGHT_JUNCTION_TARGETS_5`;
8. `ACTIN_TURNOVER`;
9. `HALLMARK_E2F_TARGETS`;
10. `HALLMARK_G2M_CHECKPOINT`.

Named genes are `GATA6`, `HNF4A`, `LGR5`, `MKI67`, `EMP1`, `KLF4`, `EPCAM`, `ZEB1`, `CDH17`, `TJP1`, `CLDN2`, `CLDN3`, `CLDN4`, `CLDN7`, `CRB3` and `F11R`.

## Frozen estimands and inference

### Primary contrast

The primary effect is the equal-weight mean of `Met2 - Primary2` and `Met3 - Primary3`, estimated with a linear model containing generation and organoid origin. For every programme, report:

1. the adjusted metastatic-derived minus primary-derived score difference;
2. a 95% bootstrap interval obtained by resampling whole libraries within each generation-by-origin cell;
3. an exhaustive two-sided permutation *P* value from all label allocations that preserve three metastatic-derived and three primary-derived labels within each generation (20 × 20 = 400 allocations);
4. generation-specific mean differences for rounds 2 and 3;
5. leave-one-gene-out direction ranges for REC, broad tight junction and the focused five-gene set;
6. a limma/`camera` competitive gene-set test for the adjusted origin coefficient.

BH FDR is controlled across the ten programmes separately for module-score and `camera` families. Gene-level effects use the same adjusted design and report limma log2 fold change, 95% interval, moderated *P*, all-gene FDR and named-gene-family FDR.

### Contextual trajectory

`Met1`, `Met2` and `Met3` group means are shown descriptively. Because each round arose after another transplantation cycle and lacks a complete round-1 primary control, these values are not treated as a longitudinal patient trajectory or an independent trend test.

## Interpretation without hard gates

- Concordant REC/HRC, junction and low-cycle directions in both controlled rounds would indicate functional-context convergence with the complete replacement-associated phenotype.
- REC/HRC enrichment without junction preservation, or junction enrichment without REC/HRC, would indicate partial convergence.
- Opposite or round-inconsistent effects would show that general liver-metastatic selection does not inevitably create the replacement-associated state.
- Regardless of direction, GSE290752 contains no HGP label and cannot establish a cause, marker or validation of rHGP.

Base random seed is 42; bootstrap resamples whole libraries 10,000 times.
