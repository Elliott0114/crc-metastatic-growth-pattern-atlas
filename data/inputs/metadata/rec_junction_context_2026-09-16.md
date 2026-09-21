# Focused REC junction-context extension

Date: 2026-09-16. Status: exploratory extension specified before calculating the new component results. Earlier HRC, junction, REC50 and public perturbation results have already been inspected. This is not a preregistration or a claim of untouched external validation.

## Question and scope

Do the membrane-strand and polarity-complex components contained in the existing junction programme have different relationships to regenerative transcription, HGP, regional context and perturbation? Can measured regulatory evidence support a more specific explanation than a single junction total score?

The analysis does not seek a new cell type, prognostic signature, optimal gene count, HGP classifier or inferred micro-to-macro trajectory. Existing manuscripts and released analyses remain the historical baseline. All new results, including null and opposing effects, are retained together.

## Definitions fixed without outcome selection

- Start from the existing 30-gene MSigDB 2025.1.Hs Reactome tight-junction set. Use all measurable members, without a highly-variable-gene restriction.
- Claudin component: all 21 CLDN-prefixed members of this set.
- Polarity-complex component: CRB3, PALS1, PARD3, PARD6A, PARD6B, PARD6G, PATJ and PRKCI. This is a mechanistic partition based on the Reactome PAR and Crumbs complexes, not a newly discovered module. F11R is retained individually and in the broad score, not forced into either component.
- Resolve PALS1 and its previous symbol MPP5 as the same entity, keeping dataset-specific source identifiers in the coverage table. Source: NCBI Gene 64398, https://www.ncbi.nlm.nih.gov/gene/64398. Other target symbols are not guessed. A matched gene is counted once.
- Biological sources: https://reactome.org/content/detail/R-HSA-420029 and https://reactome.org/content/detail/R-HSA-420661.
- Reference programmes: existing core HRC, White RSC, iCMS2/iCMS3, E2F, hypoxia, epithelial and liver-context lists. Remove all junction genes from conditioning programmes to prevent target/predictor gene reuse. No genes are selected using HGP or intervention effects.
- Gene-level results cover the complete 30-member source set and serve to inspect heterogeneity and cell origin; no significant subset replaces the fixed component definitions.

## Human single-cell analysis

Use the deposited Ogden decontaminated expression matrix. Source Cell_type and Cell_subtype labels are preserved. Profile each junction gene across source compartments by patient to distinguish tumour-epithelial expression from liver/endothelial/stromal expression.

Normalize to 10,000 per cell and log1p transform. Primary programme scores subtract controls matched in 25 mean-expression bins, with 20 controls per occupied bin. Target and conditioning sets are excluded from the control pool. Control sets for distinct scores are disjoint. Also save uncorrected mean-expression scores for a scoring sensitivity. No imputation is used.

For every patient, estimate HRC-component partial Spearman correlations in (a) epithelial cells, requiring at least 100 cells and including source-state indicators, and (b) REC cells, requiring at least 30 cells. Technical covariates are log library size and log detected-gene count. The full conditional model additionally includes junction-disjoint RSC, iCMS3-minus-iCMS2 and E2F scores. Rank scores and continuous covariates before residualization; categorical state terms are indicators. Cells estimate within-patient relationships; patients supply independent summaries.

Report equal-patient Fisher-z summaries, 10,000 patient-bootstrap intervals, exhaustive sign-flip tests, and the paired difference between component correlations. The primary family consists of the two full-conditional all-epithelial component tests and their paired difference; REC and technical-only models are sensitivities. A second score construction checks whether a difference depends on the matched-control definition. Gene-level HRC relationships form a separate family across measurable genes. Save all patient effects and leave-one-patient summary ranges. These conditional associations do not establish causal independence or functional junction assembly.

## HGP and regional projection

Use existing complete gene-expression inputs, not only previously selected REC50 genes. Report source membership and effective coverage for each component.

GSE151165: retain the 15 tumour patients and original logCPM inputs. Mean gene-z scores use sample SD. Estimate rHGP-minus-dHGP effects, exhaustive 5,005 allocations, stratified patient-bootstrap intervals, and models adjusting for junction-disjoint HRC. Epithelial and liver-context scores are added separately as sensitivities, never interpreted as measured cell fractions. Primary component inference is a two-test family; gene-level comparisons are a separate family. No feature selection or classifier fitting occurs in these 15 patients.

GSE294385: retain the 11 patient-paired micro/macro profiles. Use the existing full-gene counts and original logCPM convention. Scores use population SD across 22 profiles. Compare both fixed components with paired sign flips and patient-bootstrap intervals. Differences adjusted for disjoint HRC and epithelial/liver proxies assess conditional regional relationships. These data do not contain HGP labels or longitudinal lineage evidence.

Source feature-table inspection confirmed different assay coverage between sections. Retain the intersection across every included source feature table, as in the historical analysis; absent features are not biological zeros. Report coverage separately from detection. The covariate-adjusted HGP and paired-change model P values use approximate HC3 t tests, whereas the unadjusted primary comparisons use exhaustive randomization tests. Flag paired adjustments evaluated outside the observed covariate-change range.

For source-annotated tumour spots, examine local liver/stromal neighbourhood context if annotation-to-coordinate matching is complete. Analyse within patient/section and region, include HRC, epithelial/liver proxies and sequencing depth, and compare simple versus context-augmented models using held-out patients. Local histological neighbourhoods are not inferred molecular communication or direct cell-contact measurements. No programme expression is imputed from the REC reference. HD data are preferable if available; failed access is recorded rather than bypassed by a proxy.

Before neighbourhood results: use the six nearest grid positions within one Visium lattice step, excluding the focal spot. Require at least three annotated measured neighbours. Liver context is source hepatic-lobule annotation; stromal context pools the three source stromal labels. Use measured log-CP10k mean programme expression, with the same gene intersection in every section. Patient-wise partial rank models include section and micro/macro labels, HRC, epithelial and liver scores, log library size and log detected genes. Hold out one patient at a time in weighted linear models (each training patient has equal total weight), comparing the same baseline predictors plus regional/platform indicators against baseline plus the two neighbourhood fractions. Train predictor scaling only on training patients. Summarize held-out RMSE improvement by patient; no spot-level inference or spatial-block claim.

## Regulatory and perturbation evidence

Retrieve deposited GSE326231 GATA6 CUT&RUN and GSE326232 H3K27ac CUT&RUN count tables, plus GSE290751 serial-selection ATAC. Preserve full SOFT metadata and hashes. GSE290751 is not paired GATA6-knockout ATAC.

Use the existing human-mouse one-to-one ortholog table for new RNA projections. Keep library-level comparisons descriptive when edited-clone or parent-lineage independence is unresolved. Analyse all component members rather than picking concordant targets after observing intervention effects.

Before joining peaks to genes, verify the assembly, peak convention and gene annotation. Restrict primary local regulatory summaries to annotated promoter-proximal peaks (TSS +/-2 kb), with +/-1 kb as a location sensitivity. A nearby or overlapping peak is evidence of local occupancy/activity, not proof of enhancer-to-gene causality. Report coverage and all target outcomes; no distal nearest-gene assignments are described as validated regulatory edges. Relative H3K27ac changes are not interpreted as absolute global acetylation gains without appropriate normalization information.

A specific regulatory explanation requires consistent measured RNA/occupancy or chromatin evidence and clear mapping of experimental units. Failure to cover target promoters is uninformative, not evidence of absent regulation. Public perturbations without HGP cannot establish HGP causality.

## Explicit exploratory follow-up checks added after the primary results

Directly compare the fixed component effects in each external context; a significant association for one component and a nonsignificant association for the other do not demonstrate their difference. Keep these two additional contrasts in a separate labelled exploratory family. Inspect leave-one-gene influence without selecting a new subset. Reuse the existing E-MTAB-12022 reconstruction and >=20-epithelial-cell patient threshold as a five-patient source sensitivity, with exact group allocations and no independence claim versus the source's overlapping spatial patients. Repeat held-out prediction in macro-only spots and spots with >=5 neighbours. These follow-up checks are not presented as untouched confirmation.

Prediction intervals summarize the 11 held-out patient losses conditional on the fitted models; resampling does not refit the overlapping training folds. Fold-summary sign-flip P values are descriptive, not exact prediction-validation inference. Use gain magnitude, direction by patient, and sensitivity analyses for interpretation.

Mapping QC found that the frozen table assigns human CLDN4 to mouse Cldn13 while MGI also lists Cldn4 as a human CLDN4 ortholog (https://www.informatics.jax.org/marker/MGI:1313314; https://www.informatics.jax.org/marker/MGI:1913102). Exclude the disputed CLDN4 mapping in this new cross-species analysis. A direct Ensembl REST confirmation timed out; no replacement is guessed. The excluded Cldn13 row was not retained by RNA expression filtering, so programme-score estimates are unchanged by the exclusion. Preserve the original mapping and manuscripts.

## Completion and interpretation

Final measurement follow-up, added after viewing the complete gene heatmap: the abundant epithelial CLDN3/4/7 members had regional directions opposite to the aggregate claudin score. This prompted an outcome-blind abundance check (>=1 raw CPM in at least half the profiles) and sensitivity to low-count handling: the original fixed count prior, a library-size-scaled 0.5-count prior, and fixed 0.1/1-CPM priors. Retain every resulting score, effect and coverage; do not choose a preferred method by significance. Also inspect mean-log-expression scoring, means of individually normalized spots, and the five-patient epithelial pseudobulk source check. These are explicitly post-primary measurement checks. Use exact original source-symbol gene coverage separately to determine the impact on the current manuscript. The regional junction direction reverses under the library-scaled prior, so the old opposing HRC/junction interpretation is not retained. Independently verify the scaled-prior calculation with edgeR. Earlier primary tables remain historical/exploratory outputs, not a guarantee that their interpretation passed these later checks.

Deliver executable scripts, frozen definitions and coverage, complete patient/gene/component tables, measured regulatory evidence, a compact figure, and a Chinese synthesis that distinguishes completed results from inaccessible or unsupported evidence. Reproduce historical broad-score effects where the original definitions are unchanged, verify sample identities and nesting, and check key statistics independently.

An improvement in scientific interpretation requires reproducible component-specific information beyond the existing broad programme. Smaller P values or a more complex method are not sufficient. If the fixed components behave similarly, or differences disappear with source/state/measurement checks, report that result and do not search further partitions for significance.
