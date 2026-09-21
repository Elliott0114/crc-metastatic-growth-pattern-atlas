# Cross-platform malignant-state anchoring specification

Date: 2026-08-23  
Status: post-hoc analysis specification; not a preregistration

## Purpose

Connect the source-defined NC3 state in the 113-gene HGP ISS dataset to a full-transcriptome malignant state in the Ogden et al. paired multiome liver-metastasis cohort. The Ogden cohort has no HGP label and is therefore used to expand state biology, not to validate HGP association.

## Primary state anchor

- ISS population: NC1, NC2, NC3 and NCγ cells; primary analysis uses cells inside the source tumour mask.
- Ogden population: source-annotated epithelial cells from liver metastases.
- Statistical balancing: first sum counts within patient and state, normalize each patient-state pseudobulk to log1p counts per million over the shared panel, and then average patients equally within state.
- Eligibility: at least 20 cells per patient-state pseudobulk and at least five eligible patients per state.
- Shared-gene filter: a gene must be detected in at least 25% of eligible patient-state pseudobulks in both datasets and show a state-centroid logCPM range of at least 0.5 in both datasets.
- Platform alignment: z-score each retained gene across states separately within each dataset.
- Similarity: Pearson correlation across retained genes between each ISS and Ogden state centroid.
- Primary mapping: the Ogden state with the highest full-panel correlation to NC3.

## Interpretability checks

- Rank Ogden states using the six source NC3 markers KRT19, KRT18, CDX2, KRT20, ASCL2 and SLPI.
- Repeat the NC3 mapping after excluding each ISS patient and each Ogden patient in turn, retaining the primary shared-gene set.
- Report the best and second-best correlations and their margin; do not force a unique biological label if the margin is negligible or leave-one-patient-out mappings are unstable.

## Full-transcriptome handoff

After the anchor state is identified, create paired patient pseudobulks for that state versus all other Ogden epithelial states. Include patients with at least 20 anchor-state nuclei and at least 100 other epithelial nuclei. These counts are passed to a paired limma-voom analysis with patient blocking. The resulting program is an Ogden-defined transcriptomic expansion of NC3 and is not described as an independently discovered universal state.

## Claim boundary

This bridge can support similarity, state-program expansion and candidate regulatory interpretation. It cannot prove lineage equivalence, HGP specificity in the Ogden cohort, or causal transcription-factor activity.
