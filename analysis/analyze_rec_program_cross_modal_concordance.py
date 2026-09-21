#!/usr/bin/env python3
"""Summarize gene-level REC-program concordance across HGP scRNA and Visium."""

from __future__ import annotations

import platform
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
OGDEN = ROOT / "analysis_results" / "ogden_anchor_program" / "selected_anchor_program_top50.tsv"
SCRNA = ROOT / "analysis_results" / "rec_program_hgp_scrna_projection" / "patient_gene_pseudobulk.tsv"
SPATIAL = ROOT / "analysis_results" / "rec_program_hgp_spatial_projection" / "patient_region_gene_pseudobulk.tsv.gz"
SPEC = ROOT / "metadata" / "rec_program_cross_modal_concordance_spec_2026-08-23.md"
OUTPUT = ROOT / "analysis_results" / "rec_program_cross_modal_concordance"

MECHANISTIC_MODULES = {
    "hypoxia_mp6_overlap": [
        "ADM", "ANGPTL4", "ANKRD37", "EGLN3", "LDHA", "NDRG1", "NDUFA4L2", "SLC16A3",
    ],
    "regenerative_overlap": ["ADM", "APOL1", "DUOXA2", "ISG15", "OAS1", "SLC16A3"],
    "experimental_ap1_target_overlap": ["CYP3A5", "FHL2", "GSN", "KRT80", "PLAUR"],
    "nfkb_regulon_overlap": ["ABCG1", "BIRC3", "DUSP5", "MXD1", "PLAUR", "SDCBP2"],
}


def pairwise_superiority(rhgp: np.ndarray, dhgp: np.ndarray) -> float:
    comparisons = [
        1.0 if r_value > d_value else 0.5 if r_value == d_value else 0.0
        for r_value in rhgp
        for d_value in dhgp
    ]
    return float(np.mean(comparisons))


def patient_gene_effects(table: pd.DataFrame, prefix: str) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for gene, frame in table.groupby("gene", sort=False):
        rhgp = frame.loc[frame["hgp"] == "rHGP", "log2_cpm"].to_numpy(dtype=float)
        dhgp = frame.loc[frame["hgp"] == "dHGP", "log2_cpm"].to_numpy(dtype=float)
        if len(rhgp) == 0 or len(dhgp) == 0:
            continue
        rows.append({
            "gene": gene,
            f"{prefix}_n_rhgp": len(rhgp),
            f"{prefix}_n_dhgp": len(dhgp),
            f"{prefix}_rhgp_mean_log2_cpm": float(rhgp.mean()),
            f"{prefix}_dhgp_mean_log2_cpm": float(dhgp.mean()),
            f"{prefix}_rhgp_minus_dhgp": float(rhgp.mean() - dhgp.mean()),
            f"{prefix}_pairwise_superiority": pairwise_superiority(rhgp, dhgp),
            f"{prefix}_all_rhgp_above_all_dhgp": bool(rhgp.min() > dhgp.max()),
        })
    return pd.DataFrame(rows)


def zscore(values: pd.Series) -> pd.Series:
    standard_deviation = float(values.std(ddof=0))
    if not np.isfinite(standard_deviation) or standard_deviation == 0:
        return pd.Series(np.zeros(len(values)), index=values.index)
    return (values - float(values.mean())) / standard_deviation


def correlation_row(
    frame: pd.DataFrame,
    first: str,
    second: str,
    label: str,
) -> dict[str, object]:
    rho, p_value = spearmanr(frame[first], frame[second])
    return {
        "comparison": label,
        "n_genes": len(frame),
        "spearman_rho": float(rho),
        "nominal_p_value_descriptive": float(p_value),
        "interpretation_boundary": "genes were source-selected and correlated; P value is descriptive",
    }


def main() -> None:
    ogden = pd.read_csv(OGDEN, sep="\t").sort_values("program_rank")
    ogden["gene"] = ogden["gene"].astype(str).str.upper()
    scrna = pd.read_csv(SCRNA, sep="\t")
    scrna["gene"] = scrna["gene"].astype(str).str.upper()
    spatial = pd.read_csv(SPATIAL, sep="\t")
    spatial["gene"] = spatial["gene"].astype(str).str.upper()
    spatial = spatial[
        (spatial["qc_min_detected_genes"] == 200)
        & (spatial["restriction"] == "all_tumour_side")
        & (spatial["region"] == "near_interface_0_500um")
    ].copy()

    scrna_effects = patient_gene_effects(scrna, "scrna")
    spatial_effects = patient_gene_effects(spatial, "spatial")
    genes = (
        ogden[["gene", "program_rank", "logFC", "fraction_patients_positive"]]
        .rename(columns={
            "logFC": "ogden_rec_logfc",
            "fraction_patients_positive": "ogden_fraction_patients_positive",
        })
        .merge(scrna_effects, on="gene", how="left", validate="one_to_one")
        .merge(spatial_effects, on="gene", how="left", validate="one_to_one")
    )
    common = genes.dropna(
        subset=["scrna_rhgp_minus_dhgp", "spatial_rhgp_minus_dhgp"]
    ).copy()
    common["scrna_positive"] = common["scrna_rhgp_minus_dhgp"] > 0
    common["spatial_positive"] = common["spatial_rhgp_minus_dhgp"] > 0
    common["cross_modal_positive"] = common["scrna_positive"] & common["spatial_positive"]
    common["cross_modal_opposite"] = common["scrna_positive"] != common["spatial_positive"]
    common["scrna_effect_z"] = zscore(common["scrna_rhgp_minus_dhgp"])
    common["spatial_effect_z"] = zscore(common["spatial_rhgp_minus_dhgp"])
    common["cross_modal_support_score"] = (
        common["scrna_effect_z"] + common["spatial_effect_z"]
    ) / 2
    common["cross_modal_support_rank"] = common["cross_modal_support_score"].rank(
        ascending=False, method="min"
    ).astype(int)
    genes = genes.merge(
        common[[
            "gene", "scrna_positive", "spatial_positive", "cross_modal_positive",
            "cross_modal_opposite", "scrna_effect_z", "spatial_effect_z",
            "cross_modal_support_score", "cross_modal_support_rank",
        ]],
        on="gene",
        how="left",
        validate="one_to_one",
    )

    correlations = pd.DataFrame([
        correlation_row(
            common,
            "scrna_rhgp_minus_dhgp",
            "spatial_rhgp_minus_dhgp",
            "HGP scRNA effect versus HGP near-interface spatial effect",
        ),
        correlation_row(
            common,
            "ogden_rec_logfc",
            "scrna_rhgp_minus_dhgp",
            "Ogden REC effect versus HGP scRNA effect",
        ),
        correlation_row(
            common,
            "ogden_rec_logfc",
            "spatial_rhgp_minus_dhgp",
            "Ogden REC effect versus HGP near-interface spatial effect",
        ),
    ])

    module_rows: list[dict[str, object]] = []
    for module_name, module_genes in MECHANISTIC_MODULES.items():
        selected = common[common["gene"].isin(module_genes)].copy()
        module_rows.append({
            "module": module_name,
            "n_defined": len(module_genes),
            "n_common": len(selected),
            "genes": ",".join(selected.sort_values("program_rank")["gene"]),
            "median_scrna_effect": float(selected["scrna_rhgp_minus_dhgp"].median()),
            "median_spatial_effect": float(selected["spatial_rhgp_minus_dhgp"].median()),
            "scrna_positive_fraction": float(selected["scrna_positive"].mean()),
            "spatial_positive_fraction": float(selected["spatial_positive"].mean()),
            "cross_modal_positive_fraction": float(selected["cross_modal_positive"].mean()),
            "median_cross_modal_support_score": float(selected["cross_modal_support_score"].median()),
        })
    modules = pd.DataFrame(module_rows).sort_values(
        ["cross_modal_positive_fraction", "median_cross_modal_support_score"],
        ascending=False,
    )

    scrna_patients = set(scrna["patient"].astype(str))
    spatial_patients = set(spatial["patient"].astype(str))
    overlap = sorted(scrna_patients.intersection(spatial_patients))
    overlap_table = pd.DataFrame({
        "category": ["scrna_only", "spatial_only", "shared"],
        "patients": [
            ",".join(sorted(scrna_patients - spatial_patients)),
            ",".join(sorted(spatial_patients - scrna_patients)),
            ",".join(overlap),
        ],
        "n_patients": [
            len(scrna_patients - spatial_patients),
            len(spatial_patients - scrna_patients),
            len(overlap),
        ],
    })

    summary = pd.DataFrame([{
        "n_program_genes": len(ogden),
        "n_common_hgp_genes": len(common),
        "scrna_positive_n": int(common["scrna_positive"].sum()),
        "scrna_positive_fraction": float(common["scrna_positive"].mean()),
        "spatial_positive_n": int(common["spatial_positive"].sum()),
        "spatial_positive_fraction": float(common["spatial_positive"].mean()),
        "cross_modal_positive_n": int(common["cross_modal_positive"].sum()),
        "cross_modal_positive_fraction": float(common["cross_modal_positive"].mean()),
        "cross_modal_opposite_n": int(common["cross_modal_opposite"].sum()),
        "shared_patient_n": len(overlap),
        "scrna_spatial_effect_spearman_rho": correlations.iloc[0]["spearman_rho"],
    }])

    ranked_core = common.sort_values(
        ["cross_modal_positive", "cross_modal_support_score"], ascending=False
    ).copy()
    positive_core = ranked_core[ranked_core["cross_modal_positive"]].copy()
    top_core_genes = positive_core.head(15)["gene"].tolist()
    top_core = positive_core.head(15)[[
        "gene", "program_rank", "ogden_rec_logfc", "scrna_rhgp_minus_dhgp",
        "spatial_rhgp_minus_dhgp", "scrna_pairwise_superiority",
        "spatial_pairwise_superiority", "cross_modal_support_score",
        "cross_modal_support_rank",
    ]]

    OUTPUT.mkdir(parents=True, exist_ok=True)
    outputs = {
        "gene_level_concordance.tsv": genes.sort_values("program_rank"),
        "cross_modal_positive_core.tsv": positive_core,
        "top15_interpretable_core.tsv": top_core,
        "mechanistic_subset_summary.tsv": modules,
        "effect_correlations.tsv": correlations,
        "patient_overlap.tsv": overlap_table,
        "summary.tsv": summary,
        "run_info.tsv": pd.DataFrame([{
            "analysis_date": "2026-08-23",
            "python_version": platform.python_version(),
            "numpy_version": np.__version__,
            "pandas_version": pd.__version__,
            "statistical_unit": "patient_within_modality",
            "spatial_definition": "q200_all_tumour_side_near_0_500um",
        }]),
    }
    for filename, frame in outputs.items():
        frame.to_csv(OUTPUT / filename, sep="\t", index=False, encoding="utf-8")

    main_rho = float(correlations.iloc[0]["spearman_rho"])
    lines = [
        "# Cross-modal concordance of the frozen REC program",
        "",
        f"Of {len(common)} genes measurable in both HGP modalities, "
        f"{int(common['scrna_positive'].sum())} ({common['scrna_positive'].mean():.1%}) were higher in rHGP "
        f"epithelial pseudobulks and {int(common['spatial_positive'].sum())} "
        f"({common['spatial_positive'].mean():.1%}) were higher in the rHGP near-interface spatial band.",
        "",
        f"A total of {int(common['cross_modal_positive'].sum())} genes "
        f"({common['cross_modal_positive'].mean():.1%}) were positive in both modalities. Gene effects "
        f"correlated at Spearman rho={main_rho:+.3f}.",
        "",
        f"The 15 highest cross-modal positive genes were: {', '.join(top_core_genes)}.",
        "",
        "## Mechanistic readout",
        "",
    ]
    for row in modules.itertuples(index=False):
        lines.append(
            f"- {row.module}: {row.n_common} genes; {row.cross_modal_positive_fraction:.1%} positive "
            f"in both modalities; median effects {row.median_scrna_effect:+.3f} in scRNA and "
            f"{row.median_spatial_effect:+.3f} in spatial pseudobulks."
        )
    lines.extend([
        "",
        "## Boundary",
        "",
        f"Four patients were shared between modalities ({', '.join(overlap)}). The concordance therefore "
        "localizes a consistent epithelial/spatial signal but is not independent cohort replication. "
        "The descriptive core does not replace the frozen 50-gene program in any validation analysis.",
        "",
        "## Reproducibility",
        "",
        f"Specification: `{SPEC.relative_to(ROOT)}`.",
    ])
    (OUTPUT / "analysis_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    manifest = [
        "# Analysis output manifest",
        "",
        "Generated by `scripts/analyze_rec_program_cross_modal_concordance.py`.",
        "",
        "| File | Purpose |",
        "|---|---|",
        "| `analysis_summary.md` | Human-readable cross-modal result |",
        "| `gene_level_concordance.tsv` | All frozen genes and modality-specific effects |",
        "| `cross_modal_positive_core.tsv` | Genes positive in both HGP modalities |",
        "| `top15_interpretable_core.tsv` | Compact ranked core for interpretation |",
        "| `mechanistic_subset_summary.tsv` | Fixed pathway/regulon subset directions |",
        "| `effect_correlations.tsv` | Descriptive rank correlations |",
        "| `patient_overlap.tsv` | Cross-modality patient overlap audit |",
        "| `summary.tsv` | Key counts and fractions |",
        "| `run_info.tsv` | Reproducibility settings |",
    ]
    (OUTPUT / "_analysis_outputs.md").write_text("\n".join(manifest) + "\n", encoding="utf-8")
    print(summary.to_string(index=False))
    print(modules.to_string(index=False))
    print(top_core.to_string(index=False))


if __name__ == "__main__":
    main()
