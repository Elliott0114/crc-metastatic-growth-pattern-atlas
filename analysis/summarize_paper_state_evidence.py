#!/usr/bin/env python3
"""Build auditable Phase 2 evidence tables and a manuscript-facing summary.

Analysis date: 2026-09-01
Random seed: 42 (declared for audit consistency; tabulation is deterministic)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


np.random.seed(42)


ROOT = Path(__file__).resolve().parents[1]
OUT = (
    ROOT
    / "analysis_results"
    / "deep_biology_upgrade_2026-08-31"
    / "phase2_mechanistic_specificity"
)

PRINCIPAL = ["Hypoxia", "UPR", "iREC"]
DIRECT_PROGRAMS = {
    "CANELLAS_CORE_HRC": "Core HRC",
    "NUSSE_REGENERATIVE": "Regenerative epithelium",
    "PUBLISHED_PEMT": "Partial EMT",
    "HALLMARK_TNFA_NFKB_OGDEN_COPY": "TNF–NF-κB",
    "HALLMARK_HYPOXIA": "Hypoxia",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE": "UPR",
    "HALLMARK_GLYCOLYSIS": "Glycolysis",
    "HALLMARK_E2F_TARGETS": "E2F targets",
    "HALLMARK_G2M_CHECKPOINT": "G2M checkpoint",
    "HALLMARK_WNT_BETA_CATENIN_SIGNALING": "WNT/β-catenin",
    "HALLMARK_MYC_TARGETS_V1": "MYC targets",
    "HALLMARK_FATTY_ACID_METABOLISM": "Fatty-acid metabolism",
    "PENGWINKLER_PALMITATE_UP_TOP100": "Palmitate response",
}
IDENTITY_DISCOVERY_PROGRAMS = {
    "GOBP_EPITHELIAL_STRUCTURE_MAINTENANCE": "Epithelial structure maintenance",
    "GOBP_HOMOTYPIC_CELL_CELL_ADHESION": "Homotypic cell–cell adhesion",
    "REACTOME_TIGHT_JUNCTION_INTERACTIONS": "Tight-junction interactions",
    "REACTOME_KERATINIZATION": "Keratinization",
}
RIVAL_PROGRAMS = {
    "PAN_EPITHELIAL_7": "Pan-epithelial content",
    "HEPATOCYTE_CONTEXT_6": "Hepatocyte context",
    "HALLMARK_HYPOXIA": "Hallmark hypoxia",
    "OGDEN_HYPOXIA_MP6": "Ogden hypoxia",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE": "UPR",
    "HALLMARK_E2F_TARGETS": "E2F targets",
    "HALLMARK_G2M_CHECKPOINT": "G2M checkpoint",
    "HALLMARK_WNT_BETA_CATENIN_SIGNALING": "WNT/β-catenin",
    "HALLMARK_MYC_TARGETS_V1": "MYC targets",
    "HALLMARK_GLYCOLYSIS": "Glycolysis",
    "HALLMARK_FATTY_ACID_METABOLISM": "Fatty-acid metabolism",
    "PENGWINKLER_PALMITATE_UP_TOP100": "Palmitate response",
    "PENGWINKLER_CAUSAL_AXIS_5": "MYC–proline–COL1A1 axis",
    "GOBP_PROLINE_METABOLIC_PROCESS": "Proline metabolism",
    "REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES": "Collagen biosynthesis",
}
DISCOVERY_TERMS = {
    "REACTOME_FORMATION_OF_THE_CORNIFIED_ENVELOPE": "Cornified-envelope formation",
    "REACTOME_KERATINIZATION": "Keratinization",
    "REACTOME_TIGHT_JUNCTION_INTERACTIONS": "Tight-junction interactions",
    "GOBP_HOMOTYPIC_CELL_CELL_ADHESION": "Homotypic cell–cell adhesion",
    "GOBP_EPITHELIAL_STRUCTURE_MAINTENANCE": "Epithelial structure maintenance",
    "GOBP_PROTEIN_FOLDING_IN_ENDOPLASMIC_RETICULUM": "ER protein folding",
    "REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR": "Unfolded-protein response",
    "GOBP_RESPONSE_TO_OXYGEN_LEVELS": "Response to oxygen levels",
    "GOBP_CENTROMERE_COMPLEX_ASSEMBLY": "Centromere assembly",
    "GOBP_INTERSTRAND_CROSS_LINK_REPAIR": "Interstrand cross-link repair",
}


def fmt(value: float, digits: int = 3) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:.{digits}f}"


def p_fmt(value: float) -> str:
    if value < 0.0001:
        return f"{value:.2e}"
    return f"{value:.4f}"


def main():
    direct_genes = pd.read_csv(OUT / "ogden_direct_state_limma_all.tsv.gz", sep="\t")
    direct_gsea = pd.read_csv(OUT / "ogden_direct_state_gsea_all.tsv.gz", sep="\t")
    # Genes consistently different between REC and all three principal stress states.
    principal_genes = direct_genes[
        direct_genes["comparator_state"].isin(PRINCIPAL)
    ].copy()
    gene_rows: list[dict[str, object]] = []
    for gene, group in principal_genes.groupby("gene"):
        if group["comparator_state"].nunique() != len(PRINCIPAL):
            continue
        up = (group["log2fc_rec_minus_comparator"] > 0) & (group["fdr_bh"] < 0.05)
        down = (group["log2fc_rec_minus_comparator"] < 0) & (group["fdr_bh"] < 0.05)
        direction = "not_consistent"
        if up.all():
            direction = "REC_up_all_three"
        elif down.all():
            direction = "REC_down_all_three"
        gene_rows.append(
            {
                "gene": gene,
                "direction": direction,
                "mean_moderated_t": group["moderated_t"].mean(),
                "minimum_log2fc": group["log2fc_rec_minus_comparator"].min(),
                "maximum_log2fc": group["log2fc_rec_minus_comparator"].max(),
                "maximum_fdr_bh": group["fdr_bh"].max(),
                "hypoxia_log2fc": group.loc[
                    group["comparator_state"].eq("Hypoxia"),
                    "log2fc_rec_minus_comparator",
                ].iloc[0],
                "upr_log2fc": group.loc[
                    group["comparator_state"].eq("UPR"),
                    "log2fc_rec_minus_comparator",
                ].iloc[0],
                "irec_log2fc": group.loc[
                    group["comparator_state"].eq("iREC"),
                    "log2fc_rec_minus_comparator",
                ].iloc[0],
            }
        )
    consistent_genes = pd.DataFrame(gene_rows)
    consistent_genes = consistent_genes[
        consistent_genes["direction"].ne("not_consistent")
    ].sort_values("mean_moderated_t", ascending=False)
    consistent_genes.to_csv(
        OUT / "ogden_rec_principal_consistent_genes.tsv", sep="\t", index=False
    )

    # Frozen programme evidence for the three principal comparisons.
    targeted_key = direct_gsea[
        direct_gsea["analysis_mode"].eq("targeted_predeclared")
        & direct_gsea["comparator_state"].isin(PRINCIPAL)
        & direct_gsea["term"].isin(DIRECT_PROGRAMS)
    ].copy()
    targeted_key["program_label"] = targeted_key["term"].map(DIRECT_PROGRAMS)
    targeted_key["evidence_class"] = "targeted"
    discovery_identity = direct_gsea[
        direct_gsea["comparator_state"].isin(PRINCIPAL)
        & direct_gsea["term"].isin(IDENTITY_DISCOVERY_PROGRAMS)
        & direct_gsea["analysis_mode"].str.startswith("full_")
    ].copy()
    discovery_identity["program_label"] = discovery_identity["term"].map(
        IDENTITY_DISCOVERY_PROGRAMS
    )
    discovery_identity["evidence_class"] = "transcriptome-wide discovery"
    direct_key = pd.concat([targeted_key, discovery_identity], ignore_index=True)
    direct_key["comparator_order"] = direct_key["comparator_state"].map(
        {name: index for index, name in enumerate(PRINCIPAL)}
    )
    combined_program_order = {
        **DIRECT_PROGRAMS,
        **IDENTITY_DISCOVERY_PROGRAMS,
    }
    direct_key["program_order"] = direct_key["term"].map(
        {name: index for index, name in enumerate(combined_program_order)}
    )
    direct_key = direct_key.sort_values(["program_order", "comparator_order"])
    direct_key.to_csv(
        OUT / "table_phase2_direct_program_gsea.tsv", sep="\t", index=False
    )

    discovery_key = direct_gsea[
        direct_gsea["comparator_state"].isin(PRINCIPAL)
        & direct_gsea["term"].isin(DISCOVERY_TERMS)
    ].copy()
    discovery_key["term_label"] = discovery_key["term"].map(DISCOVERY_TERMS)
    discovery_key.to_csv(
        OUT / "table_phase2_recurrent_discovery_terms.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
