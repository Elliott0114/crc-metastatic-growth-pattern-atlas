# HGP interface-ecology analysis specification

Date: 2026-08-23  
Status: post-hoc analysis specification; not a preregistration  
Dataset: Escriva Conde et al. BJC 2026 targeted ISS, panel 2

## Scientific question

Does replacement-type colorectal liver metastasis combine two directly visible features: expansion of the source-defined NC3 malignant state and increased local adjacency between malignant cells and liver epithelial cells?

## Main analysis

- Statistical unit: canonical patient.
- Primary comparison: patients contributing only one HGP class (7 EHGP/dHGP and 8 RHGP/rHGP); P03 and P05, which contribute both classes, are excluded.
- Patient aggregation: calculate each endpoint per ROI, then average ROIs equally within patient and HGP.
- Co-primary endpoint 1: NC3 fraction among source-defined neoplastic cells in the source tumour mask.
- Co-primary endpoint 2: fraction of all source-defined neoplastic cells with at least one liver epithelial cell among their 10 nearest segmented-cell neighbours in the same ROI.
- Neoplastic cells: NC1, NC2, NC3 and NCγ.
- Liver epithelial cells: HC1, DHC1, DHC2, DHC3 and CC1.
- Effect: RHGP minus EHGP difference in patient means, with a patient-stratified bootstrap 95% confidence interval.
- Inference: exact two-sided label permutation for the difference in means. Holm adjustment is applied only across the two co-primary endpoints.
- Descriptive support: group medians and interquartile ranges, probability of superiority/rank-biserial correlation, and leave-one-patient-out direction.

The local-neighbour endpoint describes tissue architecture. It is not interpreted as molecular ligand–receptor communication, preferential attraction, or a causal cell–cell interaction.

## Secondary and sensitivity analyses

- Neighbourhood sizes of 6, 15 and 30 cells.
- Liver targets excluding CC1; damaged hepatocytes only (DHC1–DHC3); and HC1 only.
- Restriction of focal neoplastic cells to the source tumour mask.
- Restriction of target liver epithelial cells to the source liver mask.
- Cell-pooled patient estimates versus the primary equal-ROI patient estimates.
- Observed-to-random local-neighbour ratio conditional on the within-ROI abundance of target cells.
- Cluster-resolved adjacency for NC1, NC2, NC3 and NCγ; this is descriptive and is not used to claim NC3-specific contact.
- Within-patient RHGP-minus-EHGP changes for mixed-HGP patients P03 and P05; these are descriptive because n = 2.

## Interpretation contract

The main result supports a simple interface-ecology statement only if both co-primary effects point in the expected direction and their uncertainty is reported. Failure of a secondary sensitivity does not automatically negate the main architecture result; it narrows which cell definition or spatial scale is supported. No SLPI-specific claim is part of this analysis.
