# ISS NC3 cohesive and dual-interface analysis specification

Date: 2026-09-01  
Status: frozen before inspection of the endpoints defined below  
Scope: post hoc, hypothesis-generating extension of the targeted ISS analysis

## Scientific question and value

Replacement growth is morphologically defined by tumour cords that remain
cohesive while inserting between hepatic plates. The current paper identifies
an rHGP-enriched NC3/REC-like malignant state and shows that repair/plasticity
and junction transcription coexist within REC cells, but transcript abundance
does not demonstrate that these cells are spatially organised as a cohesive
population.

This analysis asks a directly visible biological question: **within rHGP
tissue, do NC3 cells occur next to other NC3 cells more often than expected
from their local abundance, and are such homotypic NC3 neighbourhoods also
positioned beside damaged hepatocytes?** A positive result would connect the
transcriptional state to the defining tissue architecture. A negative result
would prevent the manuscript from turning junctional transcription into an
unsupported claim of spatial cohesion.

The analysis remains observational. “Neighbourhood” refers to segmented-cell
proximity in the targeted ISS assay, not a measured molecular junction,
direction of signalling, migration trajectory or causal interaction.

## Fixed data and labels

- Input: `prepared_cells_panel_2.h5ad` from Escriva Conde et al.
- Published `Growth pattern`, `Clusters`, anatomical region and front-status
  annotations are retained unchanged.
- NC1, NC2, NC3 and NCγ are neoplastic states.
- DHC1-DHC3 are damaged-hepatocyte states; HC1 is the steady-state hepatocyte
  specificity control.
- Ten nearest segmented neighbours within each ROI are primary, matching the
  existing paper. Five and twenty neighbours are scale sensitivities.
- The primary state is NC3. NC1, NC2 and NCγ are contextual state controls and
  are not used to redefine the primary result.
- Canonical pseudonymous patient is the inferential unit. Multiple ROIs receive
  equal weight within a patient-HGP unit. The two patients with both HGPs remain
  paired descriptive observations and do not enter pure-HGP group inference.
- Random seed: 42. State-label exchanges: 1,000. Patient bootstraps: 10,000.

## Endpoints

For every neoplastic focal cell and neighbourhood size `k`:

1. **Homotypic-neighbour presence:** at least one of the `k` neighbours has the
   same published neoplastic-state label as the focal cell.
2. **Homotypic-neighbour fraction:** fraction of all `k` neighbour slots occupied
   by the same neoplastic state.
3. **DHC dual-interface presence:** the focal cell has both at least one
   homotypic neoplastic neighbour and at least one DHC1-DHC3 neighbour.
4. **HC1 dual-interface presence:** the corresponding endpoint using HC1,
   serving as hepatocyte-state specificity control.

The prespecified reader-facing endpoints are NC3 homotypic-neighbour presence
and NC3 DHC dual-interface presence at `k=10`. The fraction endpoint tests
density rather than presence; the HC1 endpoint tests target specificity.

## Spatially constrained null

Within each ROI, neoplastic-state labels are exchanged while coordinates,
neighbourhood geometry and all non-neoplastic labels remain fixed. Exchanges
preserve the number of each neoplastic state separately within strata defined
by published anatomical `Region`, `Liver front` status and `Tumor front`
status. Thus the null retains:

- the very large rHGP-EHGP difference in NC3 abundance;
- local cell density and tissue geometry;
- coarse tumour, stromal, liver and front localisation;
- the observed positions of damaged and steady-state hepatocytes.

For each ROI-state-endpoint, the spatial residual is the observed proportion
minus the mean exchange-null proportion. Positive residuals indicate finer
homotypic or dual-interface organisation than expected from these conditioned
features.

## Aggregation and inference

1. Calculate observed proportions, exchange-null means, residuals and empirical
   two-sided exchange P values per ROI.
2. Average ROIs equally within canonical patient-HGP units, carrying the
   matched permutation draw through the same aggregation.
3. For pure-HGP patients, report HGP-specific mean residuals and rHGP-minus-EHGP
   differences with patient bootstrap 95% confidence intervals.
4. Use exhaustive two-sided sign allocation for within-HGP residuals and
   exhaustive HGP-label allocation for between-HGP differences.
5. Apply Benjamini-Hochberg adjustment within clearly labelled endpoint/state
   families, while retaining effect estimates and raw exact P values.
6. Report all `k=5,10,20` results and all four neoplastic states; no endpoint,
   state or spatial scale is selected after viewing results.

## Interpretation without an arbitrary pass/fail gate

Evidence will be interpreted as a pattern of effect size, patient consistency,
scale robustness and biological specificity rather than a self-imposed hard
threshold.

- NC3 homotypic residuals above zero would support fine-scale NC3 spatial
  cohesion beyond abundance and coarse localisation.
- An additional DHC dual-interface residual would place cohesive NC3
  neighbourhoods at the damaged-liver interface, strengthening a concrete
  two-interface model.
- Homotypic organisation without a DHC dual residual would support tumour-
  tumour cohesion but keep damaged-liver proximity as a parallel tissue feature.
- Similar DHC and HC1 residuals would argue for general hepatocyte adjacency,
  not damage specificity.
- Comparable residuals across all neoplastic states would make the finding an
  architecture-wide property rather than NC3-specific biology.
- Residuals near zero would mean that the striking raw NC3 neighbourhood pattern
  is sufficiently explained by NC3 abundance and annotated tissue position.

Whatever the direction, the result will constrain the manuscript claim; it
will not be used to alter the frozen REC programme or the NC3-to-REC mapping.
