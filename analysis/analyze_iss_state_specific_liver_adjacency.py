#!/usr/bin/env python3
"""
Analysis: ISS NC3-specific DHC/HC1 adjacency with spatially constrained nulls.
Date: 2026-09-01
Random seed: 42
Python and package versions are written to the provenance JSON.

Neoplastic state labels are exchanged within ROI, anatomical region and both
front-status strata while target cells and neighbourhood geometry remain fixed.
Cells are measurement units; inference and intervals use equal-ROI patient-HGP
aggregates.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import platform
import re
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[1]
INPUT = (
    ROOT
    / "data_sources"
    / "EscrivaConde_BJC_2026_HGP_ISS"
    / "prepared_cells_panel_2.h5ad"
)
SPECIFICATION = ROOT / "metadata" / "deep_biology_phase2_analysis_spec_2026-09-01.md"
OUT = (
    ROOT
    / "analysis_results"
    / "deep_biology_upgrade_2026-08-31"
    / "phase2_mechanistic_specificity"
)

SEED = 42
K = 10
PERMUTATIONS = 1000
BOOTSTRAP_REPLICATES = 10_000
NEOPLASTIC = ("NC1", "NC2", "NC3", "NCγ")
DHC = ("DHC1", "DHC2", "DHC3")
HC1 = ("HC1",)
METRICS = (
    "dhc_any_neighbour",
    "dhc_neighbour_slot_fraction",
    "hc1_any_neighbour",
    "hc1_neighbour_slot_fraction",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_patient(sample: str) -> str:
    match = re.search(r"P\d+", sample)
    if match is None:
        raise ValueError(f"Cannot derive canonical patient from {sample!r}")
    return match.group(0)


def remove_self_neighbours(
    queried: np.ndarray, focal_indices: np.ndarray, k: int
) -> np.ndarray:
    if queried.ndim == 1:
        queried = queried[:, None]
    keep = queried != focal_indices[:, None]
    if np.any(keep.sum(axis=1) < k):
        raise ValueError("Failed to retrieve enough non-self neighbours")
    cleaned = np.empty((len(focal_indices), k), dtype=np.int64)
    for row in range(len(focal_indices)):
        cleaned[row] = queried[row, keep[row]][:k]
    return cleaned


def bh_adjust(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(float)
    result = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    if not valid.any():
        return result
    selected = p[valid]
    order = np.argsort(selected)
    ranked = selected[order] * len(selected) / (np.arange(len(selected)) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(ranked)
    adjusted[order] = np.minimum(ranked, 1.0)
    result[valid] = adjusted
    return result


def exact_signflip(values: np.ndarray) -> tuple[float, int]:
    values = np.asarray(values, dtype=float)
    signs = np.asarray(list(itertools.product((-1.0, 1.0), repeat=len(values))))
    null = np.mean(signs * values[None, :], axis=1)
    observed = float(np.mean(values))
    p_value = float(np.mean(np.abs(null) >= abs(observed) - 1e-15))
    return p_value, len(null)


def exact_group_permutation(values: np.ndarray, rhgp: np.ndarray) -> tuple[float, int]:
    values = np.asarray(values, dtype=float)
    rhgp = np.asarray(rhgp, dtype=bool)
    observed = float(values[rhgp].mean() - values[~rhgp].mean())
    total = float(values.sum())
    n_rhgp = int(rhgp.sum())
    n_ehgp = len(values) - n_rhgp
    extreme = 0
    assignments = 0
    for indices in itertools.combinations(range(len(values)), n_rhgp):
        rhgp_sum = float(values[list(indices)].sum())
        difference = rhgp_sum / n_rhgp - (total - rhgp_sum) / n_ehgp
        extreme += abs(difference) >= abs(observed) - 1e-15
        assignments += 1
    return extreme / assignments, assignments


def bootstrap_mean_ci(values: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    draws = rng.choice(
        values,
        size=(BOOTSTRAP_REPLICATES, len(values)),
        replace=True,
    ).mean(axis=1)
    low, high = np.quantile(draws, [0.025, 0.975])
    return float(low), float(high)


def bootstrap_group_difference_ci(
    rhgp_values: np.ndarray,
    ehgp_values: np.ndarray,
    rng: np.random.Generator,
) -> tuple[float, float]:
    rhgp_draws = rng.choice(
        rhgp_values,
        size=(BOOTSTRAP_REPLICATES, len(rhgp_values)),
        replace=True,
    ).mean(axis=1)
    ehgp_draws = rng.choice(
        ehgp_values,
        size=(BOOTSTRAP_REPLICATES, len(ehgp_values)),
        replace=True,
    ).mean(axis=1)
    low, high = np.quantile(rhgp_draws - ehgp_draws, [0.025, 0.975])
    return float(low), float(high)


def constrained_null_sums(
    outcomes: np.ndarray,
    strata: np.ndarray,
    nc3: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """Return permuted NC3 outcome sums, preserving NC3 count in every stratum."""
    selected_sums = np.zeros((PERMUTATIONS, outcomes.shape[1]), dtype=np.float64)
    for stratum in pd.unique(strata):
        indices = np.flatnonzero(strata == stratum)
        n_selected = int(nc3[indices].sum())
        if n_selected == 0:
            continue
        stratum_outcomes = outcomes[indices]
        stratum_total = stratum_outcomes.sum(axis=0)
        if n_selected == len(indices):
            selected_sums += stratum_total
            continue

        choose_complement = n_selected > len(indices) / 2
        n_draw = len(indices) - n_selected if choose_complement else n_selected
        for start in range(0, PERMUTATIONS, 100):
            stop = min(start + 100, PERMUTATIONS)
            random_values = rng.random((stop - start, len(indices)))
            chosen = np.argpartition(random_values, n_draw - 1, axis=1)[:, :n_draw]
            drawn_sums = stratum_outcomes[chosen].sum(axis=1)
            if choose_complement:
                selected_sums[start:stop] += stratum_total - drawn_sums
            else:
                selected_sums[start:stop] += drawn_sums
    return selected_sums


def analyze_sample(
    sample: str,
    frame: pd.DataFrame,
    coordinates: np.ndarray,
    rng: np.random.Generator,
) -> tuple[list[dict[str, object]], dict[str, np.ndarray]]:
    clusters = frame["Clusters"].astype(str).to_numpy()
    is_neoplastic = np.isin(clusters, NEOPLASTIC)
    focal_indices = np.flatnonzero(is_neoplastic)
    if len(focal_indices) <= K:
        raise ValueError(f"Sample {sample} has too few neoplastic cells")
    nc3 = clusters[focal_indices] == "NC3"
    if nc3.sum() == 0 or (~nc3).sum() == 0:
        raise ValueError(f"Sample {sample} lacks NC3 or comparator neoplastic cells")

    tree = cKDTree(coordinates)
    _, queried = tree.query(coordinates[focal_indices], k=K + 1, workers=-1)
    neighbours = remove_self_neighbours(np.asarray(queried), focal_indices, K)
    dhc_hits = np.isin(clusters[neighbours], DHC)
    hc1_hits = np.isin(clusters[neighbours], HC1)
    outcomes = np.column_stack(
        [
            dhc_hits.any(axis=1).astype(float),
            dhc_hits.mean(axis=1),
            hc1_hits.any(axis=1).astype(float),
            hc1_hits.mean(axis=1),
        ]
    )

    region = frame.loc[is_neoplastic, "Region"].astype(str)
    liver_front = frame.loc[is_neoplastic, "Liver front"].astype(str)
    tumour_front = frame.loc[is_neoplastic, "Tumor front"].astype(str)
    strata = (
        region.str.cat(liver_front, sep="|").str.cat(tumour_front, sep="|").to_numpy()
    )
    selected_null_sums = constrained_null_sums(outcomes, strata, nc3, rng)
    total_sums = outcomes.sum(axis=0)
    null_contrasts = (
        selected_null_sums / nc3.sum()
        - (total_sums[None, :] - selected_null_sums) / (~nc3).sum()
    )

    hgp_values = frame["Growth pattern"].astype(str).unique()
    if len(hgp_values) != 1:
        raise ValueError(f"Sample {sample} has multiple HGP labels")
    hgp = hgp_values[0]
    patient = canonical_patient(sample)
    rows: list[dict[str, object]] = []
    null_by_metric: dict[str, np.ndarray] = {}
    for metric_index, metric in enumerate(METRICS):
        nc3_mean = float(outcomes[nc3, metric_index].mean())
        other_mean = float(outcomes[~nc3, metric_index].mean())
        observed = nc3_mean - other_mean
        null = null_contrasts[:, metric_index]
        null_mean = float(null.mean())
        centered_observed = observed - null_mean
        empirical_p = float(
            (1 + np.sum(np.abs(null - null_mean) >= abs(centered_observed) - 1e-15))
            / (PERMUTATIONS + 1)
        )
        rows.append(
            {
                "sample": sample,
                "canonical_patient": patient,
                "hgp": hgp,
                "metric": metric,
                "k": K,
                "n_neoplastic": len(focal_indices),
                "n_nc3": int(nc3.sum()),
                "n_other_neoplastic": int((~nc3).sum()),
                "n_exchange_strata": int(pd.Series(strata).nunique()),
                "nc3_mean": nc3_mean,
                "other_neoplastic_mean": other_mean,
                "nc3_minus_other": observed,
                "spatial_null_mean": null_mean,
                "spatial_null_sd": float(null.std(ddof=1)),
                "spatial_null_residual": centered_observed,
                "spatial_null_empirical_p_two_sided": empirical_p,
            }
        )
        null_by_metric[metric] = null
    return rows, null_by_metric


def main() -> None:
    if not INPUT.exists() or not SPECIFICATION.exists():
        raise FileNotFoundError("Missing ISS AnnData or Phase 2 specification")
    adata = ad.read_h5ad(INPUT)
    required_obs = {
        "Sample",
        "Growth pattern",
        "Clusters",
        "Region",
        "Liver front",
        "Tumor front",
    }
    missing = sorted(required_obs.difference(adata.obs.columns))
    if missing:
        raise ValueError(f"Missing ISS observation columns: {missing}")
    if "spatial" not in adata.obsm:
        raise ValueError("ISS AnnData has no spatial coordinates")
    coordinates_all = np.asarray(adata.obsm["spatial"], dtype=float)
    if coordinates_all.shape != (adata.n_obs, 2) or not np.isfinite(coordinates_all).all():
        raise ValueError("Invalid ISS spatial coordinate matrix")

    obs = adata.obs.loc[:, sorted(required_obs)].copy().reset_index(drop=True)
    rng = np.random.default_rng(SEED)
    sample_rows: list[dict[str, object]] = []
    sample_nulls: dict[tuple[str, str], np.ndarray] = {}
    for sample in obs["Sample"].astype(str).drop_duplicates():
        selector = obs["Sample"].astype(str).eq(sample).to_numpy()
        rows, nulls = analyze_sample(
            sample,
            obs.loc[selector].reset_index(drop=True),
            coordinates_all[selector],
            rng,
        )
        sample_rows.extend(rows)
        for metric, null in nulls.items():
            sample_nulls[(sample, metric)] = null
        print(f"Completed spatial null for {sample}", flush=True)

    sample_table = pd.DataFrame(sample_rows)
    patient_rows: list[dict[str, object]] = []
    patient_nulls: dict[tuple[str, str, str], np.ndarray] = {}
    for keys, group in sample_table.groupby(
        ["canonical_patient", "hgp", "metric"], observed=True
    ):
        patient, hgp, metric = keys
        samples = group["sample"].astype(str).tolist()
        null = np.vstack([sample_nulls[(sample, metric)] for sample in samples]).mean(axis=0)
        observed = float(group["nc3_minus_other"].mean())
        null_mean = float(null.mean())
        residual = observed - null_mean
        empirical_p = float(
            (1 + np.sum(np.abs(null - null_mean) >= abs(residual) - 1e-15))
            / (PERMUTATIONS + 1)
        )
        patient_rows.append(
            {
                "canonical_patient": patient,
                "hgp": hgp,
                "metric": metric,
                "n_rois": len(samples),
                "sample_ids": ";".join(samples),
                "nc3_minus_other_equal_roi": observed,
                "spatial_null_mean_equal_roi": null_mean,
                "spatial_null_residual_equal_roi": residual,
                "spatial_null_empirical_p_two_sided": empirical_p,
            }
        )
        patient_nulls[(patient, hgp, metric)] = null
    patient_table = pd.DataFrame(patient_rows)
    hgp_counts = patient_table.groupby("canonical_patient")["hgp"].nunique()
    mixed_patients = set(hgp_counts[hgp_counts > 1].index)
    patient_table["is_mixed_patient"] = patient_table["canonical_patient"].isin(
        mixed_patients
    )

    null_draw_rows: list[pd.DataFrame] = []
    for (patient, hgp, metric), values in patient_nulls.items():
        null_draw_rows.append(
            pd.DataFrame(
                {
                    "canonical_patient": patient,
                    "hgp": hgp,
                    "metric": metric,
                    "permutation": np.arange(1, PERMUTATIONS + 1),
                    "null_nc3_minus_other_equal_roi": values,
                }
            )
        )
    null_draws = pd.concat(null_draw_rows, ignore_index=True)

    summary_rows: list[dict[str, object]] = []
    pure = patient_table.loc[~patient_table["is_mixed_patient"]].copy()
    summary_rng = np.random.default_rng(SEED + 10_000)
    for metric in METRICS:
        metric_frame = pure.loc[pure["metric"].eq(metric)]
        for value_type, column in (
            ("observed_nc3_minus_other", "nc3_minus_other_equal_roi"),
            ("spatial_null_residual", "spatial_null_residual_equal_roi"),
        ):
            for hgp in ("EHGP", "RHGP"):
                values = metric_frame.loc[metric_frame["hgp"].eq(hgp), column].to_numpy(float)
                ci_low, ci_high = bootstrap_mean_ci(values, summary_rng)
                p_value, assignments = exact_signflip(values)
                summary_rows.append(
                    {
                        "metric": metric,
                        "value_type": value_type,
                        "test": "within_hgp_mean_vs_zero",
                        "hgp": hgp,
                        "n_ehgp": len(values) if hgp == "EHGP" else np.nan,
                        "n_rhgp": len(values) if hgp == "RHGP" else np.nan,
                        "effect": float(values.mean()),
                        "bootstrap_ci_lower": ci_low,
                        "bootstrap_ci_upper": ci_high,
                        "exact_p_two_sided": p_value,
                        "n_exact_assignments": assignments,
                    }
                )

            rhgp_values = metric_frame.loc[
                metric_frame["hgp"].eq("RHGP"), column
            ].to_numpy(float)
            ehgp_values = metric_frame.loc[
                metric_frame["hgp"].eq("EHGP"), column
            ].to_numpy(float)
            ci_low, ci_high = bootstrap_group_difference_ci(
                rhgp_values, ehgp_values, summary_rng
            )
            labels = metric_frame["hgp"].eq("RHGP").to_numpy()
            p_value, assignments = exact_group_permutation(
                metric_frame[column].to_numpy(float), labels
            )
            summary_rows.append(
                {
                    "metric": metric,
                    "value_type": value_type,
                    "test": "rhgp_minus_ehgp",
                    "hgp": "RHGP_minus_EHGP",
                    "n_ehgp": len(ehgp_values),
                    "n_rhgp": len(rhgp_values),
                    "effect": float(rhgp_values.mean() - ehgp_values.mean()),
                    "bootstrap_ci_lower": ci_low,
                    "bootstrap_ci_upper": ci_high,
                    "exact_p_two_sided": p_value,
                    "n_exact_assignments": assignments,
                }
            )
    summary = pd.DataFrame(summary_rows)
    summary["fdr_bh_within_test_family"] = summary.groupby(
        ["value_type", "test", "hgp"], observed=True
    )["exact_p_two_sided"].transform(lambda values: bh_adjust(values))

    mixed = patient_table.loc[patient_table["is_mixed_patient"]].copy()
    mixed_wide = mixed.pivot(
        index=["canonical_patient", "metric"],
        columns="hgp",
        values=[
            "nc3_minus_other_equal_roi",
            "spatial_null_residual_equal_roi",
        ],
    )
    mixed_wide.columns = ["__".join(column) for column in mixed_wide.columns]
    mixed_wide = mixed_wide.reset_index()
    mixed_wide["observed_rhgp_minus_ehgp"] = (
        mixed_wide["nc3_minus_other_equal_roi__RHGP"]
        - mixed_wide["nc3_minus_other_equal_roi__EHGP"]
    )
    mixed_wide["null_residual_rhgp_minus_ehgp"] = (
        mixed_wide["spatial_null_residual_equal_roi__RHGP"]
        - mixed_wide["spatial_null_residual_equal_roi__EHGP"]
    )

    OUT.mkdir(parents=True, exist_ok=True)
    sample_table.to_csv(
        OUT / "iss_state_specific_neighbour_sample.tsv", sep="\t", index=False
    )
    patient_table.to_csv(
        OUT / "iss_state_specific_neighbour_patient_hgp.tsv", sep="\t", index=False
    )
    null_draws.to_csv(
        OUT / "iss_state_specific_neighbour_null_draws.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    summary.to_csv(
        OUT / "iss_state_specific_neighbour_group_summary.tsv", sep="\t", index=False
    )
    mixed_wide.to_csv(
        OUT / "iss_state_specific_neighbour_mixed_pairs.tsv", sep="\t", index=False
    )

    provenance = {
        "analysis_date": "2026-09-01",
        "python_version": platform.python_version(),
        "anndata_version": version("anndata"),
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "scipy_version": version("scipy"),
        "seed": SEED,
        "k_nearest_neighbours": K,
        "state_label_permutations": PERMUTATIONS,
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "exchange_strata": ["Sample", "Region", "Liver front", "Tumor front"],
        "inferential_unit": "canonical pseudonymous patient; ROIs equally weighted",
        "input_sha256": sha256(INPUT),
        "specification_sha256": sha256(SPECIFICATION),
    }
    with (OUT / "iss_state_specific_neighbour_provenance.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(provenance, handle, ensure_ascii=False, indent=2)

    print(summary.to_string(index=False))
    print(mixed_wide.to_string(index=False))


if __name__ == "__main__":
    main()
