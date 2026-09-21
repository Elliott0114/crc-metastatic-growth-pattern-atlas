#!/usr/bin/env python3
"""
Analysis: ISS neoplastic-state homotypy and damaged-liver dual interfaces.
Date: 2026-09-01
Random seed: 42
Python and key package versions are written to the provenance JSON.

Neoplastic-state labels are exchanged within ROI, anatomical region and both
front-status strata. Coordinates, neighbourhood geometry, state abundance in
each stratum and all non-neoplastic labels remain fixed. Cells are measurement
units; inference uses equal-ROI patient-HGP aggregates.
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
SPECIFICATION = (
    ROOT / "metadata" / "iss_nc3_cohesive_dual_interface_spec_2026-09-01.md"
)
OUT = (
    ROOT
    / "analysis_results"
    / "deep_biology_upgrade_2026-08-31"
    / "phase2_mechanistic_specificity"
    / "iss_cohesive_dual_interface"
)

SEED = 42
NEIGHBOURHOOD_SIZES = (5, 10, 20)
PRIMARY_K = 10
PERMUTATIONS = 1_000
BOOTSTRAP_REPLICATES = 10_000
NEOPLASTIC = ("NC1", "NC2", "NC3", "NCγ")
STATE_TO_CODE = {state: index for index, state in enumerate(NEOPLASTIC)}
DHC = ("DHC1", "DHC2", "DHC3")
HC1 = ("HC1",)
METRICS = (
    "homotypic_any_neighbour",
    "homotypic_neighbour_slot_fraction",
    "dhc_dual_interface",
    "hc1_dual_interface",
)

np.random.seed(SEED)


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


def bootstrap_mean_ci(
    values: np.ndarray, rng: np.random.Generator
) -> tuple[float, float]:
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


def endpoint_means(
    focal_codes: np.ndarray,
    neighbour_codes: np.ndarray,
    dhc_hits: np.ndarray,
    hc1_hits: np.ndarray,
) -> np.ndarray:
    same_state = neighbour_codes == focal_codes[:, None]
    homotypic_any = same_state.any(axis=1)
    outcomes = np.column_stack(
        [
            homotypic_any.astype(float),
            same_state.mean(axis=1),
            (homotypic_any & dhc_hits.any(axis=1)).astype(float),
            (homotypic_any & hc1_hits.any(axis=1)).astype(float),
        ]
    )
    means = np.full((len(NEOPLASTIC), len(METRICS)), np.nan, dtype=float)
    for code in range(len(NEOPLASTIC)):
        selected = focal_codes == code
        if selected.any():
            means[code] = outcomes[selected].mean(axis=0)
    return means


def analyze_sample(
    sample: str,
    frame: pd.DataFrame,
    coordinates: np.ndarray,
    rng: np.random.Generator,
) -> tuple[list[dict[str, object]], dict[tuple[int, str, str], np.ndarray]]:
    clusters = frame["Clusters"].astype(str).to_numpy()
    is_neoplastic = np.isin(clusters, NEOPLASTIC)
    focal_global = np.flatnonzero(is_neoplastic)
    if len(focal_global) <= max(NEIGHBOURHOOD_SIZES):
        raise ValueError(f"Sample {sample} has too few neoplastic cells")

    focal_states = clusters[focal_global]
    if set(focal_states) != set(NEOPLASTIC):
        raise ValueError(f"Sample {sample} lacks one or more neoplastic states")
    original_codes = np.asarray([STATE_TO_CODE[value] for value in focal_states])

    tree = cKDTree(coordinates)
    _, queried = tree.query(
        coordinates[focal_global],
        k=max(NEIGHBOURHOOD_SIZES) + 1,
        workers=-1,
    )
    neighbours_global = remove_self_neighbours(
        np.asarray(queried), focal_global, max(NEIGHBOURHOOD_SIZES)
    )
    global_to_neoplastic = np.full(len(frame), -1, dtype=np.int64)
    global_to_neoplastic[focal_global] = np.arange(len(focal_global))
    neighbours_local = global_to_neoplastic[neighbours_global]
    safe_neighbours_local = np.maximum(neighbours_local, 0)
    dhc_hits_all = np.isin(clusters[neighbours_global], DHC)
    hc1_hits_all = np.isin(clusters[neighbours_global], HC1)

    neoplastic_frame = frame.loc[is_neoplastic]
    strata = (
        neoplastic_frame["Region"]
        .astype(str)
        .str.cat(neoplastic_frame["Liver front"].astype(str), sep="|")
        .str.cat(neoplastic_frame["Tumor front"].astype(str), sep="|")
        .to_numpy()
    )
    stratum_indices = [
        np.flatnonzero(strata == stratum) for stratum in pd.unique(strata)
    ]

    observed = np.empty(
        (len(NEIGHBOURHOOD_SIZES), len(NEOPLASTIC), len(METRICS)), dtype=float
    )
    null = np.empty(
        (
            PERMUTATIONS,
            len(NEIGHBOURHOOD_SIZES),
            len(NEOPLASTIC),
            len(METRICS),
        ),
        dtype=np.float32,
    )

    def calculate(codes: np.ndarray, k_index: int, k: int) -> np.ndarray:
        neighbour_codes = codes[safe_neighbours_local[:, :k]]
        neighbour_codes = np.where(
            neighbours_local[:, :k] >= 0, neighbour_codes, -1
        )
        return endpoint_means(
            codes,
            neighbour_codes,
            dhc_hits_all[:, :k],
            hc1_hits_all[:, :k],
        )

    for k_index, k in enumerate(NEIGHBOURHOOD_SIZES):
        observed[k_index] = calculate(original_codes, k_index, k)

    for permutation in range(PERMUTATIONS):
        permuted_codes = original_codes.copy()
        for indices in stratum_indices:
            permuted_codes[indices] = rng.permutation(original_codes[indices])
        for k_index, k in enumerate(NEIGHBOURHOOD_SIZES):
            null[permutation, k_index] = calculate(permuted_codes, k_index, k)

    hgp_values = frame["Growth pattern"].astype(str).unique()
    if len(hgp_values) != 1:
        raise ValueError(f"Sample {sample} has multiple HGP labels")
    hgp = hgp_values[0]
    patient = canonical_patient(sample)

    rows: list[dict[str, object]] = []
    null_by_endpoint: dict[tuple[int, str, str], np.ndarray] = {}
    for k_index, k in enumerate(NEIGHBOURHOOD_SIZES):
        for state_index, state in enumerate(NEOPLASTIC):
            n_state = int((original_codes == state_index).sum())
            for metric_index, metric in enumerate(METRICS):
                draws = null[:, k_index, state_index, metric_index].astype(float)
                observed_value = float(observed[k_index, state_index, metric_index])
                null_mean = float(draws.mean())
                residual = observed_value - null_mean
                empirical_p = float(
                    (
                        1
                        + np.sum(
                            np.abs(draws - null_mean) >= abs(residual) - 1e-15
                        )
                    )
                    / (PERMUTATIONS + 1)
                )
                rows.append(
                    {
                        "sample": sample,
                        "canonical_patient": patient,
                        "hgp": hgp,
                        "k": k,
                        "state": state,
                        "metric": metric,
                        "n_neoplastic": len(focal_global),
                        "n_state": n_state,
                        "n_exchange_strata": len(stratum_indices),
                        "observed_proportion": observed_value,
                        "spatial_null_mean": null_mean,
                        "spatial_null_sd": float(draws.std(ddof=1)),
                        "spatial_null_residual": residual,
                        "spatial_null_empirical_p_two_sided": empirical_p,
                    }
                )
                null_by_endpoint[(k, state, metric)] = draws
    return rows, null_by_endpoint


def build_group_summaries(
    patient_table: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pure = patient_table.loc[~patient_table["is_mixed_patient"]].copy()
    rng = np.random.default_rng(SEED + 20_000)
    summary_rows: list[dict[str, object]] = []
    descriptive_rows: list[dict[str, object]] = []
    columns = (
        ("observed_proportion", "observed_proportion_equal_roi"),
        ("spatial_null_residual", "spatial_null_residual_equal_roi"),
    )

    for (k, state, metric), endpoint in pure.groupby(
        ["k", "state", "metric"], observed=True
    ):
        for value_type, column in columns:
            for hgp in ("EHGP", "RHGP"):
                values = endpoint.loc[endpoint["hgp"].eq(hgp), column].to_numpy(float)
                ci_low, ci_high = bootstrap_mean_ci(values, rng)
                descriptive_rows.append(
                    {
                        "k": k,
                        "state": state,
                        "metric": metric,
                        "value_type": value_type,
                        "hgp": hgp,
                        "n_patients": len(values),
                        "effect": float(values.mean()),
                        "bootstrap_ci_lower": ci_low,
                        "bootstrap_ci_upper": ci_high,
                        "n_positive": int((values > 0).sum()),
                        "n_negative": int((values < 0).sum()),
                        "n_zero": int((values == 0).sum()),
                    }
                )
                if value_type == "spatial_null_residual":
                    p_value, assignments = exact_signflip(values)
                    summary_rows.append(
                        {
                            "k": k,
                            "state": state,
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

            rhgp_values = endpoint.loc[endpoint["hgp"].eq("RHGP"), column].to_numpy(
                float
            )
            ehgp_values = endpoint.loc[endpoint["hgp"].eq("EHGP"), column].to_numpy(
                float
            )
            ci_low, ci_high = bootstrap_group_difference_ci(
                rhgp_values, ehgp_values, rng
            )
            labels = endpoint["hgp"].eq("RHGP").to_numpy()
            p_value, assignments = exact_group_permutation(
                endpoint[column].to_numpy(float), labels
            )
            summary_rows.append(
                {
                    "k": k,
                    "state": state,
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
    summary["fdr_bh_across_metrics_within_state"] = summary.groupby(
        ["k", "state", "value_type", "test", "hgp"], observed=True
    )["exact_p_two_sided"].transform(bh_adjust)
    summary["fdr_bh_across_states_within_metric"] = summary.groupby(
        ["k", "metric", "value_type", "test", "hgp"], observed=True
    )["exact_p_two_sided"].transform(bh_adjust)
    return summary, pd.DataFrame(descriptive_rows)


def main() -> None:
    if not INPUT.exists() or not SPECIFICATION.exists():
        raise FileNotFoundError("Missing ISS AnnData or frozen analysis specification")

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
    if coordinates_all.shape != (adata.n_obs, 2):
        raise ValueError("Unexpected ISS spatial coordinate matrix shape")
    if not np.isfinite(coordinates_all).all():
        raise ValueError("ISS spatial coordinates contain non-finite values")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("ISS observation or variable identifiers are not unique")

    obs = adata.obs.loc[:, sorted(required_obs)].copy().reset_index(drop=True)
    if obs.isna().any().any():
        raise ValueError("Required ISS metadata contain missing values")
    rng = np.random.default_rng(SEED)
    sample_rows: list[dict[str, object]] = []
    sample_nulls: dict[tuple[str, int, str, str], np.ndarray] = {}
    for sample in obs["Sample"].astype(str).drop_duplicates():
        selector = obs["Sample"].astype(str).eq(sample).to_numpy()
        rows, nulls = analyze_sample(
            sample,
            obs.loc[selector].reset_index(drop=True),
            coordinates_all[selector],
            rng,
        )
        sample_rows.extend(rows)
        for (k, state, metric), values in nulls.items():
            sample_nulls[(sample, k, state, metric)] = values
        print(f"Completed cohesive-interface spatial null for {sample}", flush=True)

    sample_table = pd.DataFrame(sample_rows)
    patient_rows: list[dict[str, object]] = []
    patient_nulls: dict[tuple[str, str, int, str, str], np.ndarray] = {}
    for keys, group in sample_table.groupby(
        ["canonical_patient", "hgp", "k", "state", "metric"], observed=True
    ):
        patient, hgp, k, state, metric = keys
        samples = group["sample"].astype(str).tolist()
        null = np.vstack(
            [sample_nulls[(sample, k, state, metric)] for sample in samples]
        ).mean(axis=0)
        observed = float(group["observed_proportion"].mean())
        null_mean = float(null.mean())
        residual = observed - null_mean
        empirical_p = float(
            (
                1
                + np.sum(
                    np.abs(null - null_mean) >= abs(residual) - 1e-15
                )
            )
            / (PERMUTATIONS + 1)
        )
        patient_rows.append(
            {
                "canonical_patient": patient,
                "hgp": hgp,
                "k": k,
                "state": state,
                "metric": metric,
                "n_rois": len(samples),
                "sample_ids": ";".join(samples),
                "observed_proportion_equal_roi": observed,
                "spatial_null_mean_equal_roi": null_mean,
                "spatial_null_residual_equal_roi": residual,
                "spatial_null_empirical_p_two_sided": empirical_p,
            }
        )
        patient_nulls[(patient, hgp, k, state, metric)] = null
    patient_table = pd.DataFrame(patient_rows)
    hgp_counts = patient_table.groupby("canonical_patient")["hgp"].nunique()
    mixed_patients = set(hgp_counts[hgp_counts > 1].index)
    patient_table["is_mixed_patient"] = patient_table["canonical_patient"].isin(
        mixed_patients
    )

    summary, descriptive = build_group_summaries(patient_table)

    mixed = patient_table.loc[patient_table["is_mixed_patient"]].copy()
    mixed_wide = mixed.pivot(
        index=["canonical_patient", "k", "state", "metric"],
        columns="hgp",
        values=["observed_proportion_equal_roi", "spatial_null_residual_equal_roi"],
    )
    mixed_wide.columns = ["__".join(column) for column in mixed_wide.columns]
    mixed_wide = mixed_wide.reset_index()
    mixed_wide["observed_rhgp_minus_ehgp"] = (
        mixed_wide["observed_proportion_equal_roi__RHGP"]
        - mixed_wide["observed_proportion_equal_roi__EHGP"]
    )
    mixed_wide["residual_rhgp_minus_ehgp"] = (
        mixed_wide["spatial_null_residual_equal_roi__RHGP"]
        - mixed_wide["spatial_null_residual_equal_roi__EHGP"]
    )

    primary = summary.loc[
        summary["k"].eq(PRIMARY_K)
        & summary["state"].eq("NC3")
        & summary["metric"].isin(
            ["homotypic_any_neighbour", "dhc_dual_interface", "hc1_dual_interface"]
        )
    ].copy()

    OUT.mkdir(parents=True, exist_ok=True)
    sample_table.to_csv(
        OUT / "iss_cohesive_dual_interface_sample.tsv", sep="\t", index=False
    )
    patient_table.to_csv(
        OUT / "iss_cohesive_dual_interface_patient_hgp.tsv", sep="\t", index=False
    )
    summary.to_csv(
        OUT / "iss_cohesive_dual_interface_group_summary.tsv", sep="\t", index=False
    )
    summary.to_csv(OUT / "iss_cohesive_dual_interface_group_summary.csv", index=False)
    descriptive.to_csv(
        OUT / "iss_cohesive_dual_interface_group_descriptive.tsv",
        sep="\t",
        index=False,
    )
    mixed_wide.to_csv(
        OUT / "iss_cohesive_dual_interface_mixed_pairs.tsv", sep="\t", index=False
    )
    primary.to_csv(
        OUT / "iss_cohesive_dual_interface_primary_summary.tsv",
        sep="\t",
        index=False,
    )

    provenance = {
        "analysis_date": "2026-09-01",
        "python_version": platform.python_version(),
        "anndata_version": version("anndata"),
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "scipy_version": version("scipy"),
        "seed": SEED,
        "neighbourhood_sizes": list(NEIGHBOURHOOD_SIZES),
        "primary_k": PRIMARY_K,
        "state_label_permutations": PERMUTATIONS,
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "exchange_strata": ["Sample", "Region", "Liver front", "Tumor front"],
        "inferential_unit": "canonical pseudonymous patient; ROIs equally weighted",
        "input_sha256": sha256(INPUT),
        "specification_sha256": sha256(SPECIFICATION),
    }
    with (OUT / "iss_cohesive_dual_interface_provenance.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(provenance, handle, ensure_ascii=False, indent=2)

    manifest = """# Analysis Outputs

Generated: 2026-09-01  
Study type: patient-level spatial group comparison with constrained label exchange

## Tables

- `iss_cohesive_dual_interface_sample.tsv` -- ROI-level observed and null estimates.
- `iss_cohesive_dual_interface_patient_hgp.tsv` -- Equal-ROI patient-HGP estimates.
- `iss_cohesive_dual_interface_group_summary.tsv` / `.csv` -- Exact tests and bootstrap intervals.
- `iss_cohesive_dual_interface_group_descriptive.tsv` -- HGP-specific means and patient consistency.
- `iss_cohesive_dual_interface_mixed_pairs.tsv` -- Paired descriptive mixed-HGP patients.
- `iss_cohesive_dual_interface_primary_summary.tsv` -- NC3 reader-facing k=10 endpoints.

## Provenance

- `iss_cohesive_dual_interface_provenance.json` -- Input/specification hashes and software versions.
"""
    with (OUT / "_analysis_outputs.md").open("w", encoding="utf-8") as handle:
        handle.write(manifest)

    print("\nPrimary NC3 k=10 results")
    print(primary.to_string(index=False))


if __name__ == "__main__":
    main()
