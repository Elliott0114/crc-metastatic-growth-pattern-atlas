#!/usr/bin/env python3
"""Prepare traceable ISS cell coordinates for the publication validation figure.

This script performs data extraction only. Plotting and all visual exports remain
in the R figure-production workflow.
"""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
H5AD = (
    ROOT
    / "data_sources"
    / "EscrivaConde_BJC_2026_HGP_ISS"
    / "prepared_cells_panel_2.h5ad"
)
STATE = ROOT / "analysis_results" / "slpi_bjc2026_state_composition"
OUT = ROOT / "analysis_results" / "slpi_publication_figure_augmentation"

MIXED_PATIENTS = {"P03", "P05"}
NEOPLASTIC_CLUSTERS = {"NC1", "NC2", "NC3", "NCγ"}
MIN_NEOPLASTIC_CELLS = 1_000


def select_representative_samples() -> pd.DataFrame:
    patients = pd.read_csv(STATE / "patient_state_endpoints.tsv", sep="\t")
    samples = pd.read_csv(STATE / "sample_state_endpoints.tsv", sep="\t")

    independent_patients = patients.loc[
        ~patients["canonical_patient"].isin(MIXED_PATIENTS)
    ].copy()
    group_means = independent_patients.groupby("hgp", observed=True)[
        "nc3_fraction_neoplastic"
    ].mean()

    candidates = samples.loc[
        ~samples["canonical_patient"].isin(MIXED_PATIENTS)
        & samples["n_neoplastic_cells"].ge(MIN_NEOPLASTIC_CELLS)
    ].copy()
    candidates["group_patient_mean_nc3_fraction"] = candidates["hgp"].map(
        group_means
    )
    candidates["absolute_distance_from_group_mean"] = (
        candidates["nc3_fraction_neoplastic"]
        - candidates["group_patient_mean_nc3_fraction"]
    ).abs()

    selected = (
        candidates.sort_values(
            [
                "hgp",
                "absolute_distance_from_group_mean",
                "n_neoplastic_cells",
                "Sample",
            ],
            ascending=[True, True, False, True],
        )
        .groupby("hgp", observed=True, as_index=False)
        .head(1)
        .sort_values("hgp")
        .copy()
    )
    selected["display_hgp"] = selected["hgp"].map(
        {"EHGP": "dHGP", "RHGP": "rHGP"}
    )
    selected["selection_rule"] = (
        "Independent-patient sample with >=1000 neoplastic cells and NC3 fraction "
        "closest to the group patient-level mean"
    )
    columns = [
        "Sample",
        "canonical_patient",
        "hgp",
        "display_hgp",
        "n_neoplastic_cells",
        "n_NC3_cells",
        "nc3_fraction_neoplastic",
        "group_patient_mean_nc3_fraction",
        "absolute_distance_from_group_mean",
        "selection_rule",
    ]
    return selected[columns]


def extract_spatial_cells(selected: pd.DataFrame) -> pd.DataFrame:
    adata = ad.read_h5ad(H5AD, backed="r")
    selected_samples = set(selected["Sample"])
    mask = adata.obs["Sample"].astype(str).isin(selected_samples).to_numpy()
    obs = adata.obs.loc[
        mask,
        [
            "Sample",
            "Growth pattern",
            "Tumor",
            "Tumor front",
            "Region",
            "Clusters",
        ],
    ].copy()
    coordinates = np.asarray(adata.obsm["spatial"])[mask]
    obs.insert(0, "cell_id", obs.index.astype(str))
    obs["x"] = coordinates[:, 0]
    obs["y"] = coordinates[:, 1]
    obs["display_hgp"] = obs["Growth pattern"].astype(str).map(
        {"EHGP": "dHGP", "RHGP": "rHGP"}
    )
    tumour = obs["Tumor"].astype(bool)
    clusters = obs["Clusters"].astype(str)
    obs["display_class"] = np.select(
        [
            tumour & clusters.eq("NC3"),
            tumour & clusters.isin(NEOPLASTIC_CLUSTERS),
            tumour,
        ],
        ["NC3", "Other neoplastic", "Tumour-region context"],
        default="Non-tumour context",
    )
    return obs.reset_index(drop=True)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    selected = select_representative_samples()
    if set(selected["display_hgp"]) != {"dHGP", "rHGP"}:
        raise RuntimeError("Expected one representative sample for each HGP group")
    cells = extract_spatial_cells(selected)
    if set(cells["Sample"].astype(str)) != set(selected["Sample"]):
        raise RuntimeError("Spatial extraction did not recover both selected samples")

    selected.to_csv(OUT / "representative_selection.tsv", sep="\t", index=False)
    cells.to_csv(
        OUT / "representative_spatial_cells.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    metadata = {
        "input": str(H5AD.relative_to(ROOT)),
        "selection_input": str(
            (STATE / "sample_state_endpoints.tsv").relative_to(ROOT)
        ),
        "excluded_mixed_patients": sorted(MIXED_PATIENTS),
        "minimum_neoplastic_cells": MIN_NEOPLASTIC_CELLS,
        "selected_samples": selected["Sample"].tolist(),
        "n_exported_cells": int(len(cells)),
        "visual_backend": "R (this script performs data extraction only)",
    }
    (OUT / "run_info.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(selected.to_string(index=False))
    print(f"Exported {len(cells):,} cells for the R publication figure")


if __name__ == "__main__":
    main()
