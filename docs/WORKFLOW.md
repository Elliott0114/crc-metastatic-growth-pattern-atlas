# Workflow And Provenance

The executable order is `metadata/analysis_steps.json`. Scripts run in a new,
isolated workspace; their historical internal output names are preserved so
that statistical implementations do not need a broad rewrite.

| Stage | Scope |
| --- | --- |
| Input preparation | Checksummed assay matrices, source labels, fixed programme/control definitions and annotations |
| ISS and mapping | HGP composition, liver adjacency, conditional neighbourhood nulls, 94-gene reference mapping |
| Discovery | Ogden REC50, patient-state contrasts, full original enrichment families, C8 and regulatory comparisons |
| Supporting context | Clinical, protein, within-cell and paired-state calculations retained in the current supplement |
| HGP projections | Epithelial and spatial reconstruction, bulk counts and source interface coefficients |
| Regional projection | Patient-paired micro/macro expression, source counts and retained rival definitions |
| Junction analyses | Components, all members, neighbourhoods, count eligibility and normalization |
| Reinforcement | Size/measurement matching, same-section models and bulk score/proxy sensitivities |
| External validation | Common targets/controls, five independent studies, REML/Hartung-Knapp summaries |
| Regional adjustment | Baseline and source-adjusted count models, both normalizations, eleven and ten pairs |
| Functional supplement | MRTX1133, acute Plexin B2, GATA6 and serial-selection source analyses |
| Final assembly | `materialize_sources.py`, final 21 September plots, explicit workbook source recipes |

`metadata/reproduction_map.tsv` links each of the 52 panels and 127 worksheets
to result files, producing scripts and ultimate input files. Detailed membership,
ordering and missing-value conventions are in the workbook recipes. Summary
worksheet text is a fixed presentation schema; its numerical cells are formatted
from fresh, full-precision results, never copied out of an old Excel workbook.

Some older calculations remain necessary because the current workbooks retain
them. Conversely, manuscript editing, submission packing, unused ATAC/perturbation
exploration and raw-read workflows are outside this repository's purpose.

## Historical REC44 Table

The current `Supplementary_Table_4.xlsx / REC44_cross_context` retains the spatial
column calculated before KRT20 was removed from the epithelial region markers.
All 44 retained values were reproduced with that earlier marker definition.
The contemporary corrected scripts alone do not reproduce this inherited table.

`reproduce_retained_rec44_spatial.py` therefore recomputes only this historical
column from the same matrices. It is isolated from the corrected main HGP
analysis. The old marker definition overlaps REC50 and must not be interpreted
as the corrected, target-disjoint analysis. The manuscript baseline has not been
silently edited; see the validation report before reusing this table scientifically.

## Verification

- Quick: pixels/dimensions for 11 PNGs, every panel source table, 127 worksheets,
  worksheet ordering, numerical manuscript statements and source checksums.
- Full: the same checks using newly computed results, explicit input provenance,
  module completion and retained stricter model reconstruction checks.
- Discrete identifiers, membership order, booleans and missingness are exact.
  Numerical results use `atol=1e-8`, `rtol=1e-6`; prose rounding checks have their
  own recorded display precision. PDF/ZIP timestamps are not statistical values.
- Source indices use current relative paths and runtime hashes. Panel order
  within each figure is checked; figure groups in the index follow the current
  redraw order, not the original manuscript's historical editing order.
- `verify` exits nonzero on failure and leaves its machine-readable report.
  Quick verification never implies that models were refitted.
