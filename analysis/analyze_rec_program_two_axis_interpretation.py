#!/usr/bin/env python3
"""Integrate HGP and metastatic-outgrowth effects for the frozen REC program."""

from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
HGP_INPUT = (
    ROOT
    / "analysis_results"
    / "rec_program_cross_modal_concordance"
    / "gene_level_concordance.tsv"
)
OUTGROWTH_INPUT = (
    ROOT
    / "analysis_results"
    / "gse294385_rec_program_extension"
    / "gene_level_paired_effects.tsv"
)
PROGRAM_INPUT = (
    ROOT
    / "analysis_results"
    / "ogden_anchor_program"
    / "selected_anchor_program_top50.tsv"
)
SPEC = ROOT / "metadata" / "rec_program_two_axis_interpretation_spec_2026-08-23.md"
OUTPUT = ROOT / "analysis_results" / "rec_program_two_axis_interpretation"

MODULE_COLUMNS = {
    "experimental_ap1_target_overlap": "member_experimental_ap1_target_overlap",
    "hypoxia_mp6_overlap": "member_hypoxia_mp6_overlap",
    "regenerative_overlap": "member_regenerative_overlap",
    "nfkb_regulon_overlap": "member_nfkb_regulon_overlap",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def classify(row: pd.Series) -> str:
    if bool(row["hgp_supported"]) and bool(row["outgrowth_supported"]):
        return "shared_hgp_outgrowth"
    if bool(row["hgp_supported"]):
        return "hgp_biased"
    if bool(row["outgrowth_supported"]):
        return "outgrowth_biased"
    return "unresolved_or_opposed"


def main() -> None:
    hgp = pd.read_csv(HGP_INPUT, sep="\t")
    outgrowth = pd.read_csv(OUTGROWTH_INPUT, sep="\t")
    program = pd.read_csv(PROGRAM_INPUT, sep="\t")

    for table in (hgp, outgrowth, program):
        table["gene"] = table["gene"].astype(str).str.upper()

    selected_outgrowth = outgrowth[outgrowth["member_rec_equal_gene"]].copy()
    merged = (
        program[["gene", "program_rank", "logFC", "fraction_patients_positive"]]
        .rename(columns={"logFC": "ogden_rec_logfc"})
        .merge(
            hgp[[
                "gene",
                "scrna_rhgp_minus_dhgp",
                "spatial_rhgp_minus_dhgp",
                "scrna_positive",
                "spatial_positive",
                "cross_modal_positive",
            ]],
            on="gene",
            how="inner",
            validate="one_to_one",
        )
        .merge(
            selected_outgrowth[[
                "gene",
                "mean_macro_minus_micro",
                "median_macro_minus_micro",
                "positive_patient_n",
                "positive_patient_fraction",
                "paired_standardized_effect",
                "exact_signflip_p_two_sided",
                "BH_FDR_all_selected_genes",
                *MODULE_COLUMNS.values(),
            ]],
            on="gene",
            how="inner",
            validate="one_to_one",
        )
    )

    merged["hgp_supported"] = merged["cross_modal_positive"].astype(bool)
    merged["outgrowth_supported"] = (
        (merged["mean_macro_minus_micro"] > 0)
        & (merged["positive_patient_n"] >= 8)
    )
    merged["two_axis_category"] = merged.apply(classify, axis=1)

    category_order = [
        "shared_hgp_outgrowth",
        "hgp_biased",
        "outgrowth_biased",
        "unresolved_or_opposed",
    ]
    merged["two_axis_category"] = pd.Categorical(
        merged["two_axis_category"], categories=category_order, ordered=True
    )
    merged = merged.sort_values(
        ["two_axis_category", "positive_patient_n", "mean_macro_minus_micro", "program_rank"],
        ascending=[True, False, False, True],
    )

    category_rows: list[dict[str, object]] = []
    for category in category_order:
        selected = merged[merged["two_axis_category"] == category]
        category_rows.append({
            "two_axis_category": category,
            "n_genes": len(selected),
            "genes": ",".join(selected["gene"].tolist()),
            "median_scrna_rhgp_minus_dhgp": float(
                selected["scrna_rhgp_minus_dhgp"].median()
            ),
            "median_spatial_rhgp_minus_dhgp": float(
                selected["spatial_rhgp_minus_dhgp"].median()
            ),
            "median_macro_minus_micro": float(
                selected["mean_macro_minus_micro"].median()
            ),
        })
    category_summary = pd.DataFrame(category_rows)

    module_rows: list[dict[str, object]] = []
    for module, column in MODULE_COLUMNS.items():
        selected = merged[merged[column].astype(bool)].copy()
        counts = selected["two_axis_category"].value_counts().reindex(
            category_order, fill_value=0
        )
        module_rows.append({
            "module": module,
            "n_genes": len(selected),
            "genes": ",".join(selected.sort_values("program_rank")["gene"]),
            "shared_hgp_outgrowth_n": int(counts["shared_hgp_outgrowth"]),
            "hgp_biased_n": int(counts["hgp_biased"]),
            "outgrowth_biased_n": int(counts["outgrowth_biased"]),
            "unresolved_or_opposed_n": int(counts["unresolved_or_opposed"]),
            "median_scrna_rhgp_minus_dhgp": float(
                selected["scrna_rhgp_minus_dhgp"].median()
            ),
            "median_spatial_rhgp_minus_dhgp": float(
                selected["spatial_rhgp_minus_dhgp"].median()
            ),
            "median_macro_minus_micro": float(
                selected["mean_macro_minus_micro"].median()
            ),
        })
    module_summary = pd.DataFrame(module_rows)

    OUTPUT.mkdir(parents=True, exist_ok=True)
    merged.to_csv(OUTPUT / "gene_two_axis_classification.tsv", sep="\t", index=False)
    category_summary.to_csv(OUTPUT / "category_summary.tsv", sep="\t", index=False)
    module_summary.to_csv(OUTPUT / "mechanistic_module_summary.tsv", sep="\t", index=False)

    provenance = {
        "analysis_date": "2026-08-23",
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "statistical_unit": "patient within each source analysis",
        "classification_role": "descriptive integration",
        "outgrowth_support_rule": "mean_macro_minus_micro > 0 and positive_patient_n >= 8 of 11",
        "inputs": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (HGP_INPUT, OUTGROWTH_INPUT, PROGRAM_INPUT, SPEC)
        },
    }
    with (OUTPUT / "provenance.json").open("w", encoding="utf-8") as handle:
        json.dump(provenance, handle, ensure_ascii=False, indent=2)

    shared = merged[merged["two_axis_category"] == "shared_hgp_outgrowth"]
    hgp_biased = merged[merged["two_axis_category"] == "hgp_biased"]
    outgrowth_biased = merged[merged["two_axis_category"] == "outgrowth_biased"]
    lines = [
        "# Two-axis interpretation of the frozen REC program",
        "",
        f"The analysis retained {len(merged)} genes measurable in the two HGP modalities and the paired micro-to-macro dataset.",
        "",
        f"- Shared HGP/outgrowth: {len(shared)} genes ({', '.join(shared['gene'])}).",
        f"- HGP-biased: {len(hgp_biased)} genes ({', '.join(hgp_biased['gene'])}).",
        f"- Outgrowth-biased: {len(outgrowth_biased)} genes ({', '.join(outgrowth_biased['gene'])}).",
        f"- Unresolved or opposed: {len(merged) - len(shared) - len(hgp_biased) - len(outgrowth_biased)} genes.",
        "",
        "## Mechanistic module pattern",
        "",
    ]
    for row in module_summary.itertuples(index=False):
        lines.append(
            f"- {row.module}: shared={row.shared_hgp_outgrowth_n}, "
            f"HGP-biased={row.hgp_biased_n}, outgrowth-biased={row.outgrowth_biased_n}, "
            f"unresolved/opposed={row.unresolved_or_opposed_n}; median macro-minus-micro "
            f"effect {row.median_macro_minus_micro:+.3f}."
        )
    lines.extend([
        "",
        "## Interpretation boundary",
        "",
        "The categories are a direction-and-consistency summary, not fitted latent states or causal trajectories. HGP and outgrowth were measured in different datasets, and GSE294385 has no HGP labels.",
    ])
    (OUTPUT / "analysis_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    manifest = pd.DataFrame({
        "file": [
            "gene_two_axis_classification.tsv",
            "category_summary.tsv",
            "mechanistic_module_summary.tsv",
            "provenance.json",
            "analysis_summary.md",
        ],
        "role": [
            "gene-level two-axis evidence",
            "four-category summary",
            "mechanistic subset summary",
            "input and environment provenance",
            "human-readable result summary",
        ],
    })
    manifest.to_csv(OUTPUT / "_analysis_outputs.md", sep="\t", index=False)


if __name__ == "__main__":
    main()
