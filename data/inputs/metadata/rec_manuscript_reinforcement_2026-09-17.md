# Fixed specification: focused manuscript reinforcement

Authorized by the user's 2026-09-17 request to implement the reviewed reinforcement plan. Existing analysis and manuscript versions remain immutable. This specification is written before calculating the new target effects.

## Component measurement checks

- Reuse the 14 eligible Ogden patients, existing cell identifiers, HRC score, conditional design, and fixed junction definitions.
- Draw 500 unique 8-member subsets uniformly without replacement within the 20 measured claudins, using seed 20260917. Compare each with all 8 polarity members.
- Expression/detection matching uses patient-equal mean log1p(CP10k) and detection fraction, transformed with natural log and logit (floor 1e-6), then standardized across the 28 measured component genes. A permissible pair differs by no more than 0.5 SD on either feature. Select the maximum-cardinality, minimum-squared-distance one-to-one match, with dummy cost 100 and forbidden cost 10000. Matching never uses HRC correlations or HGP effects.
- Score two target conventions: mean normalized target expression; and that mean minus the component's original fixed control background. HRC and conditioning covariates stay fixed in both conventions. The fixed background isolates target membership; it is not rematched to each subset. Full-set references accompany both conventions.
- Rank residuals adjust for epithelial-state indicators, log library, detected genes, RSC, iCMS and E2F exactly as the established model. Report per-patient Fisher-z differences. For the 500-subset analysis average differences within each patient before patient bootstrap/sign-flip inference. The between-subset 2.5–97.5% range is a perturbation range, not a confidence interval.
- Use 10,000 patient bootstraps. Four secondary difference tests (size and measurement matching, each under two target conventions) form one BH family. Fewer than 3 matched pairs makes the matching check descriptive/unresolved. Preserve unmatched members and full primary results regardless of outcome.

## Same-section count analysis

- Use source Layer3 annotations and the existing 24-section sample QC. Of 12 sections containing both tumour regions, choose one per patient by descending microtumour spot count, then descending macrotumour spot count, then sample identifier. This selects 10 patients without expression-effect selection.
- Sum raw counts by selected section and tumour region. Preserve full-library offsets and the genes assayed in all sections of the established regional analysis. Validate spot counts against source QC.
- Use the established CPM>=1 in at least half the profiles and total counts>=15 rule; patient+macro design; edgeR robust dispersion/QL gene tests, unshrunk logFC; fry directional and mixed tests. TMM primary and library-only sensitivity use the same genes. Four programme tests form each directional/mixed family. Compare shared eligible gene directions with the 11-patient model. Do not count multiple sections as independent patients.

## Bulk identifiers and normalization diagnostics

- Canonicalize MPP5 to PALS1 in the established bulk logCPM projection and programme lists; retain the projection's normalization to isolate symbol harmonization. Recalculate full/disjoint junction scores, all 5005 patient allocations, 10,000 within-HGP bootstraps, patient omissions, separate epithelial/liver adjustment and CAMERA. Preserve HRC, E2F/G2M and REC50 definitions.
- Record the distinct normalization used by the score projection and raw-count models as explicit analysis definitions, not iteration history. Keep current count models unchanged.
- Describe norm factors, high-abundance transcript fractions, M–A shifts and epithelial/liver proxies for bulk and HGP spatial profiles. Diagnostics do not select a preferred method by significance.

## Histology and delivery

- Audit original TIFF/Space Ranger image dimensions and coordinate correspondence; make de-identified image/spot review materials without programme scores or HGP keys. Use verified scale information only. Do not generate expert pathology labels.
- Retain all 30 junction members in an integrated evidence figure, distinguishing unassayed, ineligible and estimated entries. Different estimands use separate scales.
- Integrate completed results into a new manuscript, supplement and figure/table bundle. Preserve biological breadth, source attribution, patient units, null results, and unresolved localization/function boundaries. Existing user approval covers the specified title/heading/figure restructuring.
