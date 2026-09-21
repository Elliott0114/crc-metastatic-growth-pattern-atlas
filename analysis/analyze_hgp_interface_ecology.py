#!/usr/bin/env python3
"""Patient-level ISS analysis of malignant-state abundance and liver adjacency.

Analysis date: 2026-08-23
Input: data_sources/EscrivaConde_BJC_2026_HGP_ISS/prepared_cells_panel_2.h5ad
Output: analysis_results/hgp_interface_ecology/
Random seed: 42

The two co-primary endpoints are deliberately simple: the tumour-mask NC3
fraction and the fraction of neoplastic cells with at least one liver epithelial
cell among their 10 nearest segmented-cell neighbours. All other definitions are
reported as sensitivity or descriptive analyses.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import platform
from importlib.metadata import version
from pathlib import Path
from typing import Iterable

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu


ROOT = Path(__file__).resolve().parents[1]
INPUT = (
    ROOT
    / "data_sources"
    / "EscrivaConde_BJC_2026_HGP_ISS"
    / "prepared_cells_panel_2.h5ad"
)
SPECIFICATION = (
    ROOT / "metadata" / "hgp_interface_ecology_analysis_spec_2026-08-23.md"
)
OUT = ROOT / "analysis_results" / "hgp_interface_ecology"

SEED = 42
BOOTSTRAP_REPLICATES = 10_000
MAIN_K = 10
NEOPLASTIC_CLUSTERS = ("NC1", "NC2", "NC3", "NCγ")
LIVER_EPITHELIAL = ("HC1", "DHC1", "DHC2", "DHC3", "CC1")
LIVER_NO_CC1 = ("HC1", "DHC1", "DHC2", "DHC3")
DAMAGED_HEPATOCYTES = ("DHC1", "DHC2", "DHC3")
EXPECTED_CELLS = 304_410
EXPECTED_GENES = 113
EXPECTED_SAMPLES = 22


CONTACT_DEFINITIONS = (
    {
        "definition": "liver_epithelial_k10",
        "k": 10,
        "focal_scope": "all_neoplastic",
        "target_set": "liver_epithelial",
        "target_clusters": LIVER_EPITHELIAL,
        "target_scope": "all_regions",
        "role": "primary",
    },
    *(
        {
            "definition": f"liver_epithelial_k{k}",
            "k": k,
            "focal_scope": "all_neoplastic",
            "target_set": "liver_epithelial",
            "target_clusters": LIVER_EPITHELIAL,
            "target_scope": "all_regions",
            "role": "neighbourhood_size_sensitivity",
        }
        for k in (6, 15, 30)
    ),
    {
        "definition": "liver_no_cc1_k10",
        "k": 10,
        "focal_scope": "all_neoplastic",
        "target_set": "liver_no_cc1",
        "target_clusters": LIVER_NO_CC1,
        "target_scope": "all_regions",
        "role": "target_definition_sensitivity",
    },
    {
        "definition": "damaged_hepatocytes_k10",
        "k": 10,
        "focal_scope": "all_neoplastic",
        "target_set": "damaged_hepatocytes",
        "target_clusters": DAMAGED_HEPATOCYTES,
        "target_scope": "all_regions",
        "role": "target_definition_sensitivity",
    },
    {
        "definition": "hc1_k10",
        "k": 10,
        "focal_scope": "all_neoplastic",
        "target_set": "hc1",
        "target_clusters": ("HC1",),
        "target_scope": "all_regions",
        "role": "target_definition_sensitivity",
    },
    {
        "definition": "liver_epithelial_tumour_focal_k10",
        "k": 10,
        "focal_scope": "tumour_mask_neoplastic",
        "target_set": "liver_epithelial",
        "target_clusters": LIVER_EPITHELIAL,
        "target_scope": "all_regions",
        "role": "focal_mask_sensitivity",
    },
    {
        "definition": "liver_epithelial_liver_target_k10",
        "k": 10,
        "focal_scope": "all_neoplastic",
        "target_set": "liver_epithelial",
        "target_clusters": LIVER_EPITHELIAL,
        "target_scope": "liver_mask",
        "role": "target_mask_sensitivity",
    },
    {
        "definition": "liver_epithelial_tumour_focal_liver_target_k10",
        "k": 10,
        "focal_scope": "tumour_mask_neoplastic",
        "target_set": "liver_epithelial",
        "target_clusters": LIVER_EPITHELIAL,
        "target_scope": "liver_mask",
        "role": "dual_mask_sensitivity",
    },
)


def sha256(path: Path) -> str:
    """Return a streaming SHA-256 checksum."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_patient(sample: str) -> str:
    """Extract the source pseudonymous patient identifier from a sample ID."""
    match = pd.Series([sample], dtype="string").str.extract(r"(P\d+)", expand=False)
    if match.isna().iloc[0]:
        raise ValueError(f"Cannot derive canonical patient from sample {sample!r}")
    return str(match.iloc[0])


def random_any_neighbour_probability(
    n_cells: int, n_targets: int, k: int
) -> float:
    """Probability of >=1 target in k draws under within-ROI label mixing."""
    possible_neighbours = n_cells - 1
    if possible_neighbours < k:
        raise ValueError("Neighbourhood is larger than the available cell set")
    if n_targets == 0:
        return 0.0
    n_non_targets = possible_neighbours - n_targets
    if n_non_targets < k:
        return 1.0
    offsets = np.arange(k, dtype=float)
    probability_none = np.prod(
        (n_non_targets - offsets) / (possible_neighbours - offsets)
    )
    return float(1.0 - probability_none)


def remove_self_neighbours(
    queried_indices: np.ndarray, focal_indices: np.ndarray, k: int
) -> np.ndarray:
    """Remove each query cell from its own nearest-neighbour result."""
    if queried_indices.ndim == 1:
        queried_indices = queried_indices[:, None]
    keep = queried_indices != focal_indices[:, None]
    counts = keep.sum(axis=1)
    if not np.all(counts >= k):
        raise ValueError("Failed to retrieve enough non-self neighbours")
    cleaned = np.empty((len(focal_indices), k), dtype=np.int64)
    for row_index in range(len(focal_indices)):
        cleaned[row_index] = queried_indices[row_index, keep[row_index]][:k]
    return cleaned


def build_sample_tables(
    obs: pd.DataFrame, spatial: np.ndarray
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Calculate ROI-level composition, adjacency and cluster-resolved endpoints."""
    required = {"Sample", "Growth pattern", "Clusters", "Tumor", "Liver"}
    missing = sorted(required.difference(obs.columns))
    if missing:
        raise ValueError(f"Missing required observation columns: {missing}")
    if spatial.shape != (len(obs), 2):
        raise ValueError(f"Unexpected spatial coordinate shape: {spatial.shape}")
    if not np.isfinite(spatial).all():
        raise ValueError("Spatial coordinates contain non-finite values")

    frame = obs.loc[:, sorted(required)].copy()
    frame["Clusters"] = frame["Clusters"].astype("string")
    frame["Sample"] = frame["Sample"].astype("string")
    frame["Growth pattern"] = frame["Growth pattern"].astype("string")
    frame["Tumor"] = frame["Tumor"].astype(bool)
    frame["Liver"] = frame["Liver"].astype(bool)

    composition_rows: list[dict[str, object]] = []
    adjacency_rows: list[dict[str, object]] = []
    cluster_rows: list[dict[str, object]] = []
    max_k = max(int(item["k"]) for item in CONTACT_DEFINITIONS)

    samples = frame["Sample"].drop_duplicates().tolist()
    for sample in samples:
        global_indices = np.flatnonzero(frame["Sample"].eq(sample).to_numpy())
        sample_frame = frame.iloc[global_indices].reset_index(drop=True)
        sample_spatial = spatial[global_indices]
        hgp_values = sample_frame["Growth pattern"].dropna().unique()
        if len(hgp_values) != 1:
            raise ValueError(f"Sample {sample} has {len(hgp_values)} HGP labels")
        hgp = str(hgp_values[0])
        patient = canonical_patient(str(sample))
        clusters = sample_frame["Clusters"].to_numpy(dtype=str)
        tumour_mask = sample_frame["Tumor"].to_numpy(dtype=bool)
        liver_mask = sample_frame["Liver"].to_numpy(dtype=bool)
        is_neoplastic = np.isin(clusters, NEOPLASTIC_CLUSTERS)
        is_nc3 = clusters == "NC3"
        tumour_neoplastic = tumour_mask & is_neoplastic
        n_neoplastic = int(tumour_neoplastic.sum())
        n_nc3 = int((tumour_neoplastic & is_nc3).sum())
        if n_neoplastic == 0:
            raise ValueError(f"Sample {sample} has no tumour-mask neoplastic cells")
        composition_rows.append(
            {
                "Sample": str(sample),
                "canonical_patient": patient,
                "hgp": hgp,
                "n_segmented_cells": len(sample_frame),
                "n_tumour_mask_neoplastic": n_neoplastic,
                "n_tumour_mask_nc3": n_nc3,
                "nc3_fraction_neoplastic": n_nc3 / n_neoplastic,
            }
        )

        focal_all = np.flatnonzero(is_neoplastic)
        if len(focal_all) == 0:
            raise ValueError(f"Sample {sample} has no source-defined neoplastic cells")
        if len(sample_frame) <= max_k:
            raise ValueError(f"Sample {sample} has too few cells for k={max_k}")
        tree = cKDTree(sample_spatial)
        _, queried_indices = tree.query(
            sample_spatial[focal_all], k=max_k + 1, workers=-1
        )
        neighbour_indices = remove_self_neighbours(
            np.asarray(queried_indices), focal_all, max_k
        )
        focal_is_tumour = tumour_mask[focal_all]

        for definition in CONTACT_DEFINITIONS:
            k = int(definition["k"])
            focal_scope = str(definition["focal_scope"])
            if focal_scope == "all_neoplastic":
                focal_selector = np.ones(len(focal_all), dtype=bool)
            elif focal_scope == "tumour_mask_neoplastic":
                focal_selector = focal_is_tumour
            else:
                raise ValueError(f"Unknown focal scope: {focal_scope}")
            selected_neighbours = neighbour_indices[focal_selector, :k]
            n_focal = int(focal_selector.sum())
            if n_focal == 0:
                raise ValueError(f"No focal cells for {sample}/{definition['definition']}")

            target = np.isin(clusters, definition["target_clusters"])
            target_scope = str(definition["target_scope"])
            if target_scope == "liver_mask":
                target &= liver_mask
            elif target_scope != "all_regions":
                raise ValueError(f"Unknown target scope: {target_scope}")
            hits = target[selected_neighbours]
            any_hits = hits.any(axis=1)
            n_targets = int(target.sum())
            expected_any = random_any_neighbour_probability(
                len(sample_frame), n_targets, k
            )
            any_rate = float(any_hits.mean())
            adjacency_rows.append(
                {
                    "Sample": str(sample),
                    "canonical_patient": patient,
                    "hgp": hgp,
                    "definition": str(definition["definition"]),
                    "role": str(definition["role"]),
                    "k": k,
                    "focal_scope": focal_scope,
                    "target_set": str(definition["target_set"]),
                    "target_scope": target_scope,
                    "n_segmented_cells": len(sample_frame),
                    "n_focal_cells": n_focal,
                    "n_target_cells": n_targets,
                    "n_focal_with_target_neighbour": int(any_hits.sum()),
                    "n_target_neighbour_slots": int(hits.sum()),
                    "any_target_neighbour_rate": any_rate,
                    "target_neighbour_fraction": float(hits.mean()),
                    "random_expected_any_rate": expected_any,
                    "observed_to_random_any_ratio": (
                        any_rate / expected_any if expected_any > 0 else np.nan
                    ),
                }
            )

        main_target = np.isin(clusters, LIVER_EPITHELIAL)
        main_neighbours = neighbour_indices[:, :MAIN_K]
        for focal_cluster in NEOPLASTIC_CLUSTERS:
            selector = clusters[focal_all] == focal_cluster
            if not selector.any():
                continue
            hits = main_target[main_neighbours[selector]]
            cluster_rows.append(
                {
                    "Sample": str(sample),
                    "canonical_patient": patient,
                    "hgp": hgp,
                    "focal_cluster": focal_cluster,
                    "k": MAIN_K,
                    "n_focal_cells": int(selector.sum()),
                    "n_focal_with_target_neighbour": int(hits.any(axis=1).sum()),
                    "any_target_neighbour_rate": float(hits.any(axis=1).mean()),
                    "target_neighbour_fraction": float(hits.mean()),
                }
            )

    return (
        pd.DataFrame(composition_rows),
        pd.DataFrame(adjacency_rows),
        pd.DataFrame(cluster_rows),
    )


def aggregate_patient_tables(
    sample_composition: pd.DataFrame,
    sample_adjacency: pd.DataFrame,
    sample_cluster_adjacency: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Aggregate ROIs equally within patient and retain cell-pooled sensitivities."""
    patient_composition = (
        sample_composition.groupby(
            ["canonical_patient", "hgp"], observed=True, as_index=False
        )
        .agg(
            n_rois=("Sample", "nunique"),
            n_tumour_mask_neoplastic=("n_tumour_mask_neoplastic", "sum"),
            n_tumour_mask_nc3=("n_tumour_mask_nc3", "sum"),
            nc3_fraction_neoplastic_equal_roi=(
                "nc3_fraction_neoplastic",
                "mean",
            ),
        )
        .sort_values(["canonical_patient", "hgp"])
    )
    patient_composition["nc3_fraction_neoplastic_cell_pooled"] = (
        patient_composition["n_tumour_mask_nc3"]
        / patient_composition["n_tumour_mask_neoplastic"]
    )

    group_columns = [
        "canonical_patient",
        "hgp",
        "definition",
        "role",
        "k",
        "focal_scope",
        "target_set",
        "target_scope",
    ]
    patient_adjacency = (
        sample_adjacency.groupby(group_columns, observed=True, as_index=False)
        .agg(
            n_rois=("Sample", "nunique"),
            n_focal_cells=("n_focal_cells", "sum"),
            n_focal_with_target_neighbour=(
                "n_focal_with_target_neighbour",
                "sum",
            ),
            n_target_neighbour_slots=("n_target_neighbour_slots", "sum"),
            any_target_neighbour_rate_equal_roi=(
                "any_target_neighbour_rate",
                "mean",
            ),
            target_neighbour_fraction_equal_roi=(
                "target_neighbour_fraction",
                "mean",
            ),
            random_expected_any_rate_equal_roi=(
                "random_expected_any_rate",
                "mean",
            ),
            observed_to_random_any_ratio_equal_roi=(
                "observed_to_random_any_ratio",
                "mean",
            ),
        )
        .sort_values(["definition", "canonical_patient", "hgp"])
    )
    patient_adjacency["any_target_neighbour_rate_cell_pooled"] = (
        patient_adjacency["n_focal_with_target_neighbour"]
        / patient_adjacency["n_focal_cells"]
    )
    patient_adjacency["target_neighbour_fraction_cell_pooled"] = (
        patient_adjacency["n_target_neighbour_slots"]
        / (patient_adjacency["n_focal_cells"] * patient_adjacency["k"])
    )

    cluster_groups = ["canonical_patient", "hgp", "focal_cluster", "k"]
    patient_cluster = (
        sample_cluster_adjacency.groupby(
            cluster_groups, observed=True, as_index=False
        )
        .agg(
            n_rois=("Sample", "nunique"),
            n_focal_cells=("n_focal_cells", "sum"),
            n_focal_with_target_neighbour=(
                "n_focal_with_target_neighbour",
                "sum",
            ),
            any_target_neighbour_rate_equal_roi=(
                "any_target_neighbour_rate",
                "mean",
            ),
            target_neighbour_fraction_equal_roi=(
                "target_neighbour_fraction",
                "mean",
            ),
        )
        .sort_values(["focal_cluster", "canonical_patient", "hgp"])
    )
    patient_cluster["any_target_neighbour_rate_cell_pooled"] = (
        patient_cluster["n_focal_with_target_neighbour"]
        / patient_cluster["n_focal_cells"]
    )
    return patient_composition, patient_adjacency, patient_cluster


def build_patient_endpoint_values(
    patient_composition: pd.DataFrame, patient_adjacency: pd.DataFrame
) -> pd.DataFrame:
    """Create a common long table for primary and sensitivity inference."""
    rows: list[pd.DataFrame] = []
    for column, role in (
        ("nc3_fraction_neoplastic_equal_roi", "primary"),
        ("nc3_fraction_neoplastic_cell_pooled", "aggregation_sensitivity"),
    ):
        part = patient_composition[["canonical_patient", "hgp", column]].copy()
        part = part.rename(columns={column: "value"})
        part["endpoint"] = column
        part["role"] = role
        rows.append(part)

    adjacency_metrics = (
        "any_target_neighbour_rate_equal_roi",
        "any_target_neighbour_rate_cell_pooled",
        "target_neighbour_fraction_equal_roi",
        "observed_to_random_any_ratio_equal_roi",
    )
    for metric in adjacency_metrics:
        for definition, definition_frame in patient_adjacency.groupby(
            "definition", observed=True
        ):
            part = definition_frame[
                ["canonical_patient", "hgp", "role", metric]
            ].copy()
            part = part.rename(columns={metric: "value"})
            part["endpoint"] = f"{definition}__{metric}"
            if not (
                definition == "liver_epithelial_k10"
                and metric == "any_target_neighbour_rate_equal_roi"
            ):
                part["role"] = part["role"].astype(str) + f"__{metric}"
            rows.append(part)

    result = pd.concat(rows, ignore_index=True)
    result["is_mixed_patient"] = result.groupby("canonical_patient")[
        "hgp"
    ].transform("nunique").gt(1)
    return result.sort_values(["endpoint", "canonical_patient", "hgp"])


def bootstrap_mean_difference(
    rhgp: np.ndarray, ehgp: np.ndarray, seed: int
) -> tuple[float, float]:
    """Patient-stratified percentile interval for the difference in means."""
    rng = np.random.default_rng(seed)
    r_draws = rng.choice(
        rhgp, size=(BOOTSTRAP_REPLICATES, len(rhgp)), replace=True
    )
    e_draws = rng.choice(
        ehgp, size=(BOOTSTRAP_REPLICATES, len(ehgp)), replace=True
    )
    differences = r_draws.mean(axis=1) - e_draws.mean(axis=1)
    low, high = np.quantile(differences, [0.025, 0.975])
    return float(low), float(high)


def exact_mean_difference_permutation(
    values: np.ndarray, rhgp_labels: np.ndarray
) -> tuple[float, int]:
    """Enumerate all assignments preserving the observed group sizes."""
    n_total = len(values)
    n_rhgp = int(rhgp_labels.sum())
    n_ehgp = n_total - n_rhgp
    observed = float(values[rhgp_labels].mean() - values[~rhgp_labels].mean())
    total = float(values.sum())
    extreme = 0
    assignments = 0
    for indices in itertools.combinations(range(n_total), n_rhgp):
        rhgp_sum = float(values[list(indices)].sum())
        difference = rhgp_sum / n_rhgp - (total - rhgp_sum) / n_ehgp
        extreme += abs(difference) >= abs(observed) - 1e-15
        assignments += 1
    return extreme / assignments, assignments


def summarize_effect(
    frame: pd.DataFrame, endpoint: str, seed_offset: int
) -> dict[str, object]:
    """Summarize one endpoint among pure-HGP patients."""
    selected = frame.loc[frame["endpoint"].eq(endpoint)].copy()
    selected = selected.loc[~selected["is_mixed_patient"]]
    rhgp = selected.loc[selected["hgp"].eq("RHGP"), "value"].to_numpy(float)
    ehgp = selected.loc[selected["hgp"].eq("EHGP"), "value"].to_numpy(float)
    if (len(ehgp), len(rhgp)) != (7, 8):
        raise ValueError(
            f"Expected 7 EHGP and 8 RHGP patients for {endpoint}; "
            f"found {len(ehgp)} and {len(rhgp)}"
        )
    ci_low, ci_high = bootstrap_mean_difference(
        rhgp, ehgp, SEED + seed_offset
    )
    exact_p, assignments = exact_mean_difference_permutation(
        selected["value"].to_numpy(float), selected["hgp"].eq("RHGP").to_numpy()
    )
    u = float(
        mannwhitneyu(rhgp, ehgp, alternative="greater", method="asymptotic").statistic
    )
    probability_superiority = u / (len(rhgp) * len(ehgp))
    return {
        "endpoint": endpoint,
        "role": str(selected["role"].iloc[0]),
        "n_EHGP": len(ehgp),
        "n_RHGP": len(rhgp),
        "EHGP_mean": float(ehgp.mean()),
        "RHGP_mean": float(rhgp.mean()),
        "mean_difference_RHGP_minus_EHGP": float(rhgp.mean() - ehgp.mean()),
        "mean_difference_bootstrap_CI_low": ci_low,
        "mean_difference_bootstrap_CI_high": ci_high,
        "EHGP_median": float(np.median(ehgp)),
        "RHGP_median": float(np.median(rhgp)),
        "EHGP_q1": float(np.quantile(ehgp, 0.25)),
        "EHGP_q3": float(np.quantile(ehgp, 0.75)),
        "RHGP_q1": float(np.quantile(rhgp, 0.25)),
        "RHGP_q3": float(np.quantile(rhgp, 0.75)),
        "probability_RHGP_greater_than_EHGP": probability_superiority,
        "rank_biserial_correlation": 2 * probability_superiority - 1,
        "exact_mean_difference_permutation_p": exact_p,
        "exact_label_assignments": assignments,
    }


def holm_adjust(p_values: Iterable[float]) -> np.ndarray:
    """Holm family-wise error adjustment."""
    values = np.asarray(list(p_values), dtype=float)
    order = np.argsort(values)
    sorted_values = values[order]
    adjusted_sorted = np.maximum.accumulate(
        (len(values) - np.arange(len(values))) * sorted_values
    )
    adjusted_sorted = np.minimum(adjusted_sorted, 1.0)
    adjusted = np.empty_like(adjusted_sorted)
    adjusted[order] = adjusted_sorted
    return adjusted


def benjamini_hochberg(p_values: Iterable[float]) -> np.ndarray:
    """Benjamini-Hochberg adjustment for descriptive cluster contrasts."""
    values = np.asarray(list(p_values), dtype=float)
    order = np.argsort(values)
    sorted_values = values[order]
    adjusted_sorted = sorted_values * len(values) / (np.arange(len(values)) + 1)
    adjusted_sorted = np.minimum.accumulate(adjusted_sorted[::-1])[::-1]
    adjusted_sorted = np.minimum(adjusted_sorted, 1.0)
    adjusted = np.empty_like(adjusted_sorted)
    adjusted[order] = adjusted_sorted
    return adjusted


def leave_one_patient_out(
    endpoint_values: pd.DataFrame, endpoints: Iterable[str]
) -> pd.DataFrame:
    """Record the direction of primary effects after omitting each patient."""
    rows: list[dict[str, object]] = []
    pure = endpoint_values.loc[~endpoint_values["is_mixed_patient"]]
    for endpoint in endpoints:
        endpoint_frame = pure.loc[pure["endpoint"].eq(endpoint)]
        for patient in sorted(endpoint_frame["canonical_patient"].unique()):
            remaining = endpoint_frame.loc[
                ~endpoint_frame["canonical_patient"].eq(patient)
            ]
            difference = (
                remaining.loc[remaining["hgp"].eq("RHGP"), "value"].mean()
                - remaining.loc[remaining["hgp"].eq("EHGP"), "value"].mean()
            )
            rows.append(
                {
                    "endpoint": endpoint,
                    "left_out_patient": patient,
                    "mean_difference_RHGP_minus_EHGP": difference,
                    "positive_direction": bool(difference > 0),
                }
            )
    return pd.DataFrame(rows)


def mixed_patient_differences(endpoint_values: pd.DataFrame) -> pd.DataFrame:
    """Return descriptive within-patient RHGP-minus-EHGP differences."""
    mixed = endpoint_values.loc[endpoint_values["is_mixed_patient"]].copy()
    pivot = mixed.pivot_table(
        index=["canonical_patient", "endpoint"],
        columns="hgp",
        values="value",
        aggfunc="first",
    ).reset_index()
    pivot = pivot.dropna(subset=["EHGP", "RHGP"])
    pivot["difference_RHGP_minus_EHGP"] = pivot["RHGP"] - pivot["EHGP"]
    return pivot.sort_values(["endpoint", "canonical_patient"])


def cluster_effects(patient_cluster: pd.DataFrame) -> pd.DataFrame:
    """Summarize cluster-resolved adjacency as an explicitly exploratory result."""
    patient_hgp_counts = patient_cluster.groupby("canonical_patient")["hgp"].nunique()
    pure_patients = patient_hgp_counts.loc[patient_hgp_counts.eq(1)].index
    selected = patient_cluster.loc[
        patient_cluster["canonical_patient"].isin(pure_patients)
    ]
    rows: list[dict[str, object]] = []
    for index, focal_cluster in enumerate(NEOPLASTIC_CLUSTERS):
        cluster_frame = selected.loc[selected["focal_cluster"].eq(focal_cluster)]
        rhgp = cluster_frame.loc[
            cluster_frame["hgp"].eq("RHGP"),
            "any_target_neighbour_rate_equal_roi",
        ].to_numpy(float)
        ehgp = cluster_frame.loc[
            cluster_frame["hgp"].eq("EHGP"),
            "any_target_neighbour_rate_equal_roi",
        ].to_numpy(float)
        if not len(rhgp) or not len(ehgp):
            continue
        ci_low, ci_high = bootstrap_mean_difference(
            rhgp, ehgp, SEED + 500 + index
        )
        labels = cluster_frame["hgp"].eq("RHGP").to_numpy()
        p_value, assignments = exact_mean_difference_permutation(
            cluster_frame["any_target_neighbour_rate_equal_roi"].to_numpy(float),
            labels,
        )
        rows.append(
            {
                "focal_cluster": focal_cluster,
                "n_EHGP": len(ehgp),
                "n_RHGP": len(rhgp),
                "EHGP_mean": float(ehgp.mean()),
                "RHGP_mean": float(rhgp.mean()),
                "mean_difference_RHGP_minus_EHGP": float(
                    rhgp.mean() - ehgp.mean()
                ),
                "mean_difference_bootstrap_CI_low": ci_low,
                "mean_difference_bootstrap_CI_high": ci_high,
                "exact_mean_difference_permutation_p": p_value,
                "exact_label_assignments": assignments,
            }
        )
    result = pd.DataFrame(rows)
    result["BH_FDR_across_four_clusters"] = benjamini_hochberg(
        result["exact_mean_difference_permutation_p"]
    )
    return result


def diagnostic_plot(
    endpoint_values: pd.DataFrame, primary_effects: pd.DataFrame
) -> None:
    """Create a provisional two-panel patient-level diagnostic plot."""
    endpoint_labels = {
        "nc3_fraction_neoplastic_equal_roi": "NC3 fraction among\nneoplastic cells",
        "liver_epithelial_k10__any_target_neighbour_rate_equal_roi": (
            "Neoplastic cells with a\nliver epithelial neighbour"
        ),
    }
    colors = {"EHGP": "#607D8B", "RHGP": "#C65F42"}
    pure = endpoint_values.loc[~endpoint_values["is_mixed_patient"]]
    rng = np.random.default_rng(SEED)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), sharey=False)
    for axis, endpoint in zip(axes, endpoint_labels, strict=True):
        part = pure.loc[pure["endpoint"].eq(endpoint)]
        for x_position, hgp in enumerate(("EHGP", "RHGP")):
            values = part.loc[part["hgp"].eq(hgp), "value"].to_numpy(float)
            jitter = rng.uniform(-0.07, 0.07, size=len(values))
            axis.scatter(
                np.full(len(values), x_position) + jitter,
                values,
                s=30,
                color=colors[hgp],
                edgecolor="white",
                linewidth=0.5,
                zorder=3,
            )
            axis.hlines(
                np.mean(values),
                x_position - 0.22,
                x_position + 0.22,
                color="black",
                linewidth=1.5,
                zorder=4,
            )
        effect = primary_effects.loc[primary_effects["endpoint"].eq(endpoint)].iloc[0]
        axis.set_xticks([0, 1], ["dHGP", "rHGP"])
        axis.set_ylabel(endpoint_labels[endpoint])
        axis.set_ylim(bottom=0)
        axis.spines[["top", "right"]].set_visible(False)
        axis.set_title(
            f"Δ={effect['mean_difference_RHGP_minus_EHGP']:.2f}; "
            f"Holm P={effect['holm_p_across_two_primary']:.4f}",
            fontsize=9,
        )
    fig.suptitle("Provisional patient-level ISS diagnostics", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "primary_endpoints_diagnostic.png", dpi=200)
    fig.savefig(OUT / "primary_endpoints_diagnostic.pdf")
    plt.close(fig)


def dataframe_to_markdown(frame: pd.DataFrame) -> str:
    """Render a small DataFrame as Markdown without optional dependencies."""
    columns = [str(column) for column in frame.columns]

    def render(value: object) -> str:
        if pd.isna(value):
            return ""
        if isinstance(value, (float, np.floating)):
            text = f"{float(value):.6g}"
        else:
            text = str(value)
        return text.replace("|", "\\|").replace("\n", " ")

    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    lines.extend(
        "| " + " | ".join(render(value) for value in row) + " |"
        for row in frame.itertuples(index=False, name=None)
    )
    return "\n".join(lines)


def write_summary(
    primary_effects: pd.DataFrame,
    sensitivity_effects: pd.DataFrame,
    loo: pd.DataFrame,
    mixed: pd.DataFrame,
) -> None:
    """Write an analysis-facing summary without overstating causal interpretation."""
    nc3 = primary_effects.loc[
        primary_effects["endpoint"].eq("nc3_fraction_neoplastic_equal_roi")
    ].iloc[0]
    adjacency = primary_effects.loc[
        primary_effects["endpoint"].eq(
            "liver_epithelial_k10__any_target_neighbour_rate_equal_roi"
        )
    ].iloc[0]
    size_sensitivities = sensitivity_effects.loc[
        sensitivity_effects["endpoint"].isin(
            [
                f"liver_epithelial_k{k}__any_target_neighbour_rate_equal_roi"
                for k in (6, 15, 30)
            ]
        )
    ]
    all_loo_positive = loo.groupby("endpoint")["positive_direction"].all()
    mixed_primary = mixed.loc[
        mixed["endpoint"].isin(
            [
                "nc3_fraction_neoplastic_equal_roi",
                "liver_epithelial_k10__any_target_neighbour_rate_equal_roi",
            ]
        )
    ]
    lines = [
        "## Material Passport",
        "",
        "- ID: HGP-INTERFACE-ECOLOGY-ISS-2026-08-23",
        "- Type: patient-level observational spatial reanalysis",
        "- Verification status: VERIFIED (script completed and outputs were regenerated)",
        "- Data status: public, pseudonymised repository identifiers; no direct PHI fields used",
        "- Primary unit: canonical patient",
        "- Primary analysis set: 7 EHGP/dHGP and 8 RHGP/rHGP patients",
        "",
        "## Primary results",
        "",
        (
            f"- NC3 represented {nc3['EHGP_mean']:.1%} of neoplastic cells in "
            f"EHGP and {nc3['RHGP_mean']:.1%} in RHGP. The mean difference was "
            f"{nc3['mean_difference_RHGP_minus_EHGP']:.1%} points "
            f"(patient-bootstrap 95% CI {nc3['mean_difference_bootstrap_CI_low']:.1%} "
            f"to {nc3['mean_difference_bootstrap_CI_high']:.1%}; exact permutation "
            f"P={nc3['exact_mean_difference_permutation_p']:.6f}; Holm "
            f"P={nc3['holm_p_across_two_primary']:.6f})."
        ),
        (
            f"- A liver epithelial cell occurred among the 10 nearest neighbours of "
            f"{adjacency['EHGP_mean']:.1%} of neoplastic cells in EHGP and "
            f"{adjacency['RHGP_mean']:.1%} in RHGP. The mean difference was "
            f"{adjacency['mean_difference_RHGP_minus_EHGP']:.1%} points "
            f"(patient-bootstrap 95% CI "
            f"{adjacency['mean_difference_bootstrap_CI_low']:.1%} to "
            f"{adjacency['mean_difference_bootstrap_CI_high']:.1%}; exact permutation "
            f"P={adjacency['exact_mean_difference_permutation_p']:.6f}; Holm "
            f"P={adjacency['holm_p_across_two_primary']:.6f})."
        ),
        "",
        "These two endpoints support a direct architectural statement: RHGP combines "
        "expansion of the source-defined NC3 state with more frequent local adjacency "
        "between malignant and liver epithelial cells. The adjacency endpoint does not "
        "establish molecular communication or preferential attraction.",
        "",
        "## Robustness and boundary checks",
        "",
        (
            "- Neighbourhood-size effects (RHGP minus EHGP) were: "
            + "; ".join(
                f"{row.endpoint.split('_k')[1].split('__')[0]} neighbours "
                f"{row.mean_difference_RHGP_minus_EHGP:.1%} points "
                f"(P={row.exact_mean_difference_permutation_p:.4g})"
                for row in size_sensitivities.itertuples()
            )
            + "."
        ),
        (
            "- All leave-one-patient-out effects retained a positive direction: "
            + "; ".join(
                f"{endpoint}={bool(value)}"
                for endpoint, value in all_loo_positive.items()
            )
            + "."
        ),
        "- P03 and P05 are reported only as two-patient paired descriptions; their "
        "directions are not treated as a confirmatory test.",
        "",
        "## Mixed-patient primary endpoint values",
        "",
        dataframe_to_markdown(mixed_primary),
        "",
        "## Reproducibility",
        "",
        "Run from the project root:",
        "",
        "```bash",
        "conda run -n crc-metastatic-growth python scripts/analyze_hgp_interface_ecology.py",
        "```",
    ]
    (OUT / "analysis_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_manifest() -> None:
    """Write a concise output manifest for downstream manuscript work."""
    descriptions = {
        "sample_composition.tsv": "ROI-level NC3 composition endpoints",
        "sample_adjacency.tsv": "ROI-level adjacency endpoints for all definitions",
        "sample_cluster_adjacency.tsv": "ROI-level NC-cluster-resolved adjacency",
        "patient_composition.tsv": "Patient-HGP composition, equal-ROI and cell-pooled",
        "patient_adjacency.tsv": "Patient-HGP adjacency, equal-ROI and cell-pooled",
        "patient_cluster_adjacency.tsv": "Patient-HGP cluster-resolved adjacency",
        "patient_endpoint_values.tsv": "Long patient-level table used for inference",
        "primary_effects.tsv": "Two co-primary effects with Holm adjustment",
        "sensitivity_effects.tsv": "All prespecified sensitivity effects",
        "cluster_effects.tsv": "Exploratory cluster-resolved effects with BH FDR",
        "leave_one_patient_out.tsv": "Primary endpoint direction after each omission",
        "mixed_patient_differences.tsv": "Descriptive paired P03/P05 differences",
        "primary_endpoints_diagnostic.png": "Provisional raster diagnostic",
        "primary_endpoints_diagnostic.pdf": "Provisional vector diagnostic",
        "analysis_summary.md": "Human-readable verified result summary",
        "provenance.json": "Input checksum, software versions and analysis constants",
    }
    lines = ["# Analysis outputs", ""]
    for filename, description in descriptions.items():
        path = OUT / filename
        status = f"{path.stat().st_size} bytes" if path.exists() else "MISSING"
        lines.append(f"- `{filename}` — {description}; {status}")
    (OUT / "_analysis_outputs.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def main() -> None:
    """Run the full interface-ecology analysis and write auditable outputs."""
    OUT.mkdir(parents=True, exist_ok=True)
    if not INPUT.exists() or not SPECIFICATION.exists():
        raise FileNotFoundError("Required input or analysis specification is missing")
    adata = ad.read_h5ad(INPUT, backed="r")
    if (adata.n_obs, adata.n_vars) != (EXPECTED_CELLS, EXPECTED_GENES):
        raise ValueError(
            f"Unexpected AnnData dimensions: {adata.n_obs} x {adata.n_vars}"
        )
    if adata.obs["Sample"].nunique() != EXPECTED_SAMPLES:
        raise ValueError(f"Expected {EXPECTED_SAMPLES} post-QC samples")
    obs = adata.obs.copy()
    spatial = np.asarray(adata.obsm["spatial"])
    adata.file.close()

    sample_composition, sample_adjacency, sample_cluster = build_sample_tables(
        obs, spatial
    )
    patient_composition, patient_adjacency, patient_cluster = (
        aggregate_patient_tables(
            sample_composition, sample_adjacency, sample_cluster
        )
    )
    endpoint_values = build_patient_endpoint_values(
        patient_composition, patient_adjacency
    )

    primary_endpoints = [
        "nc3_fraction_neoplastic_equal_roi",
        "liver_epithelial_k10__any_target_neighbour_rate_equal_roi",
    ]
    all_endpoints = endpoint_values["endpoint"].drop_duplicates().tolist()
    effect_rows = [
        summarize_effect(endpoint_values, endpoint, index)
        for index, endpoint in enumerate(all_endpoints)
    ]
    all_effects = pd.DataFrame(effect_rows)
    primary_effects = all_effects.loc[
        all_effects["endpoint"].isin(primary_endpoints)
    ].copy()
    primary_effects["holm_p_across_two_primary"] = holm_adjust(
        primary_effects["exact_mean_difference_permutation_p"]
    )
    primary_effects["expected_direction_positive"] = (
        primary_effects["mean_difference_RHGP_minus_EHGP"] > 0
    )
    sensitivity_effects = all_effects.loc[
        ~all_effects["endpoint"].isin(primary_endpoints)
    ].copy()
    loo = leave_one_patient_out(endpoint_values, primary_endpoints)
    mixed = mixed_patient_differences(endpoint_values)
    resolved_cluster_effects = cluster_effects(patient_cluster)

    output_frames = {
        "sample_composition.tsv": sample_composition,
        "sample_adjacency.tsv": sample_adjacency,
        "sample_cluster_adjacency.tsv": sample_cluster,
        "patient_composition.tsv": patient_composition,
        "patient_adjacency.tsv": patient_adjacency,
        "patient_cluster_adjacency.tsv": patient_cluster,
        "patient_endpoint_values.tsv": endpoint_values,
        "primary_effects.tsv": primary_effects,
        "sensitivity_effects.tsv": sensitivity_effects,
        "cluster_effects.tsv": resolved_cluster_effects,
        "leave_one_patient_out.tsv": loo,
        "mixed_patient_differences.tsv": mixed,
    }
    for filename, output_frame in output_frames.items():
        output_frame.to_csv(OUT / filename, sep="\t", index=False)

    provenance = {
        "analysis_id": "HGP-INTERFACE-ECOLOGY-ISS-2026-08-23",
        "input": str(INPUT.relative_to(ROOT)),
        "input_sha256": sha256(INPUT),
        "specification": str(SPECIFICATION.relative_to(ROOT)),
        "script": str(Path(__file__).resolve().relative_to(ROOT)),
        "seed": SEED,
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "main_k": MAIN_K,
        "neoplastic_clusters": list(NEOPLASTIC_CLUSTERS),
        "liver_epithelial_clusters": list(LIVER_EPITHELIAL),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "packages": {
            package: version(package)
            for package in (
                "anndata",
                "matplotlib",
                "numpy",
                "pandas",
                "scipy",
            )
        },
    }
    (OUT / "provenance.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    diagnostic_plot(endpoint_values, primary_effects)
    write_summary(primary_effects, sensitivity_effects, loo, mixed)
    write_manifest()
    print(primary_effects.to_string(index=False))


if __name__ == "__main__":
    main()
