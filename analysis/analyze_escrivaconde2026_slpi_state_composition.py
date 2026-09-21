#!/usr/bin/env python3
"""Explain the BJC 2026 SLPI-HGP direction by neoplastic-state composition."""

from __future__ import annotations

import hashlib
import itertools
import json
import platform
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import mannwhitneyu, norm
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
INPUT = (
    ROOT
    / "data_sources"
    / "EscrivaConde_BJC_2026_HGP_ISS"
    / "prepared_cells_panel_2.h5ad"
)
CONTRACT = (
    ROOT
    / "metadata"
    / "slpi_bjc2026_state_composition_contract_2026-08-08.yml"
)
OUT = ROOT / "analysis_results" / "slpi_bjc2026_state_composition"
MIXED_PATIENTS = {"P03", "P05"}
CLUSTERS = ["NC1", "NC2", "NC3", "NCγ"]
GENES = ["SLPI", "KRT18", "GSN"]
BOOTSTRAP_SEED = 20260808
BOOTSTRAP_REPLICATES = 10_000
SOURCE_CODE_COMMIT = "06776e1b4f6d52619834e8588bfb9672eea7918d"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def dense_vector(x: object) -> np.ndarray:
    if sparse.issparse(x):
        return np.asarray(x.toarray()).ravel()
    return np.asarray(x).ravel()


def build_sample_endpoints(adata: ad.AnnData) -> pd.DataFrame:
    if adata.raw is None:
        raise ValueError("The AnnData object lacks adata.raw counts")
    missing = sorted(set(GENES).difference(adata.raw.var_names))
    if missing:
        raise ValueError(f"Missing locked genes: {missing}")

    obs = adata.obs[
        ["Sample", "Growth pattern", "Clusters", "Tumor"]
    ].copy()
    obs["canonical_patient"] = obs["Sample"].astype(str).str.extract(r"(P\d+)")
    obs["Tumor"] = obs["Tumor"].astype(bool)
    for gene in GENES:
        index = adata.raw.var_names.get_loc(gene)
        obs[gene] = dense_vector(adata.raw.X[:, index])

    rows: list[dict[str, object]] = []
    for sample, sample_frame in obs.groupby("Sample", observed=True):
        tumour = sample_frame.loc[sample_frame["Tumor"]].copy()
        neoplastic = tumour.loc[tumour["Clusters"].astype(str).isin(CLUSTERS)].copy()
        if tumour.empty or neoplastic.empty:
            raise ValueError(f"Sample {sample} lacks tumour or neoplastic cells")
        row: dict[str, object] = {
            "Sample": str(sample),
            "canonical_patient": str(sample_frame["canonical_patient"].iloc[0]),
            "hgp": str(sample_frame["Growth pattern"].iloc[0]),
            "n_tumour_cells": len(tumour),
            "n_neoplastic_cells": len(neoplastic),
            "n_NC3_cells": int(neoplastic["Clusters"].astype(str).eq("NC3").sum()),
        }
        row["nc3_fraction_all_tumour"] = row["n_NC3_cells"] / row["n_tumour_cells"]
        row["nc3_fraction_neoplastic"] = row["n_NC3_cells"] / row["n_neoplastic_cells"]

        for cluster in CLUSTERS:
            cluster_frame = neoplastic.loc[
                neoplastic["Clusters"].astype(str).eq(cluster)
            ]
            if cluster_frame.empty:
                raise ValueError(f"Sample {sample} lacks locked cluster {cluster}")
            row[f"n_{cluster}_cells"] = len(cluster_frame)
            row[f"prop_{cluster}"] = len(cluster_frame) / len(neoplastic)
            for gene in GENES:
                row[f"{gene}_{cluster}_mean"] = float(cluster_frame[gene].mean())

        for gene in GENES:
            row[f"{gene}_neoplastic_mean"] = float(neoplastic[gene].mean())
        rows.append(row)

    result = pd.DataFrame(rows)
    if result["canonical_patient"].isna().any():
        raise ValueError("Failed to derive canonical patient IDs")
    return result.sort_values(["canonical_patient", "hgp", "Sample"])


def aggregate_patient_hgp(sample_endpoints: pd.DataFrame) -> pd.DataFrame:
    metric_columns = [
        column
        for column in sample_endpoints.columns
        if column
        not in {
            "Sample",
            "canonical_patient",
            "hgp",
            "n_tumour_cells",
            "n_neoplastic_cells",
            "n_NC3_cells",
        }
    ]
    aggregations: dict[str, tuple[str, str]] = {
        "n_rois": ("Sample", "nunique"),
        "n_tumour_cells": ("n_tumour_cells", "sum"),
        "n_neoplastic_cells": ("n_neoplastic_cells", "sum"),
        "n_NC3_cells": ("n_NC3_cells", "sum"),
    }
    aggregations.update({column: (column, "mean") for column in metric_columns})
    result = (
        sample_endpoints.groupby(
            ["canonical_patient", "hgp"], observed=True, as_index=False
        )
        .agg(**aggregations)
        .sort_values(["canonical_patient", "hgp"])
    )
    return result


def add_standardised_endpoints(
    patient_hgp: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    independent = patient_hgp.loc[
        ~patient_hgp["canonical_patient"].isin(MIXED_PATIENTS)
    ].copy()
    weights = independent[[f"prop_{cluster}" for cluster in CLUSTERS]].mean()
    weights = weights / weights.sum()
    weight_table = pd.DataFrame(
        {
            "cluster": CLUSTERS,
            "common_reference_weight": [weights[f"prop_{cluster}"] for cluster in CLUSTERS],
        }
    )
    result = patient_hgp.copy()
    for gene in GENES:
        result[f"{gene}_standardised_mean"] = sum(
            weight_table.loc[
                weight_table["cluster"].eq(cluster), "common_reference_weight"
            ].iloc[0]
            * result[f"{gene}_{cluster}_mean"]
            for cluster in CLUSTERS
        )
    return result, weight_table


def hodges_lehmann(rhgp: np.ndarray, ehgp: np.ndarray) -> float:
    return float(np.median(rhgp[:, None] - ehgp[None, :]))


def exact_rank_permutation(
    values: np.ndarray, observed_rhgp: np.ndarray
) -> tuple[float, int]:
    n_total = len(values)
    n_r = int(observed_rhgp.sum())
    n_e = n_total - n_r
    ranks = pd.Series(values).rank(method="average").to_numpy()
    observed_u = float(ranks[observed_rhgp].sum() - n_r * (n_r + 1) / 2)
    center = n_r * n_e / 2
    observed_deviation = abs(observed_u - center)
    extreme = 0
    assignments = 0
    for indices in itertools.combinations(range(n_total), n_r):
        u = float(ranks[list(indices)].sum() - n_r * (n_r + 1) / 2)
        extreme += abs(u - center) >= observed_deviation - 1e-12
        assignments += 1
    return extreme / assignments, assignments


def bootstrap_hl_interval(
    rhgp: np.ndarray, ehgp: np.ndarray, seed_offset: int
) -> tuple[float, float]:
    rng = np.random.default_rng(BOOTSTRAP_SEED + seed_offset)
    estimates = np.empty(BOOTSTRAP_REPLICATES, dtype=float)
    for index in range(BOOTSTRAP_REPLICATES):
        r = rng.choice(rhgp, size=len(rhgp), replace=True)
        e = rng.choice(ehgp, size=len(ehgp), replace=True)
        estimates[index] = hodges_lehmann(r, e)
    lower, upper = np.quantile(estimates, [0.025, 0.975])
    return float(lower), float(upper)


def summarize_endpoint(
    frame: pd.DataFrame, endpoint: str, seed_offset: int
) -> dict[str, object]:
    values = frame[endpoint].astype(float)
    if values.isna().any():
        raise ValueError(f"Missing values in {endpoint}")
    rhgp = frame.loc[frame["hgp"].eq("RHGP"), endpoint].to_numpy(float)
    ehgp = frame.loc[frame["hgp"].eq("EHGP"), endpoint].to_numpy(float)
    if (len(ehgp), len(rhgp)) != (7, 8):
        raise ValueError(f"Unexpected groups for {endpoint}: {len(ehgp)}/{len(rhgp)}")
    u_result = mannwhitneyu(rhgp, ehgp, alternative="two-sided", method="asymptotic")
    u = float(mannwhitneyu(rhgp, ehgp, alternative="greater").statistic)
    exact_p, assignments = exact_rank_permutation(
        values.to_numpy(), frame["hgp"].eq("RHGP").to_numpy()
    )
    ci_low, ci_high = bootstrap_hl_interval(rhgp, ehgp, seed_offset)
    return {
        "endpoint": endpoint,
        "n_EHGP": len(ehgp),
        "n_RHGP": len(rhgp),
        "EHGP_mean": float(ehgp.mean()),
        "RHGP_mean": float(rhgp.mean()),
        "mean_difference_RHGP_minus_EHGP": float(rhgp.mean() - ehgp.mean()),
        "EHGP_median": float(np.median(ehgp)),
        "RHGP_median": float(np.median(rhgp)),
        "hodges_lehmann_RHGP_minus_EHGP": hodges_lehmann(rhgp, ehgp),
        "hodges_lehmann_bootstrap_CI_low": ci_low,
        "hodges_lehmann_bootstrap_CI_high": ci_high,
        "rank_biserial_correlation": float(2 * u / (len(rhgp) * len(ehgp)) - 1),
        "mannwhitney_asymptotic_two_sided_p": float(u_result.pvalue),
        "exact_rank_permutation_two_sided_p": exact_p,
        "exact_label_assignments": assignments,
    }


def leave_one_patient_out(frame: pd.DataFrame, endpoints: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for endpoint in endpoints:
        for patient in sorted(frame["canonical_patient"].unique()):
            remaining = frame.loc[~frame["canonical_patient"].eq(patient)]
            rhgp = remaining.loc[remaining["hgp"].eq("RHGP"), endpoint].to_numpy(float)
            ehgp = remaining.loc[remaining["hgp"].eq("EHGP"), endpoint].to_numpy(float)
            effect = hodges_lehmann(rhgp, ehgp)
            rows.append(
                {
                    "endpoint": endpoint,
                    "left_out_patient": patient,
                    "n_EHGP": len(ehgp),
                    "n_RHGP": len(rhgp),
                    "hodges_lehmann_RHGP_minus_EHGP": effect,
                    "positive_direction": effect > 0,
                }
            )
    return pd.DataFrame(rows)


def decompose_gene(
    frame: pd.DataFrame,
    gene: str,
    weights: pd.DataFrame,
    seed_offset: int,
) -> tuple[dict[str, object], pd.DataFrame]:
    raw_endpoint = f"{gene}_neoplastic_mean"
    standard_endpoint = f"{gene}_standardised_mean"
    rhgp = frame.loc[frame["hgp"].eq("RHGP")].copy()
    ehgp = frame.loc[frame["hgp"].eq("EHGP")].copy()
    observed = float(rhgp[raw_endpoint].mean() - ehgp[raw_endpoint].mean())
    within = float(rhgp[standard_endpoint].mean() - ehgp[standard_endpoint].mean())
    composition = observed - within
    rng = np.random.default_rng(BOOTSTRAP_SEED + seed_offset)
    bootstrap_rows: list[dict[str, object]] = []
    for index in range(BOOTSTRAP_REPLICATES):
        r = rhgp.iloc[rng.integers(0, len(rhgp), size=len(rhgp))]
        e = ehgp.iloc[rng.integers(0, len(ehgp), size=len(ehgp))]
        boot_observed = float(r[raw_endpoint].mean() - e[raw_endpoint].mean())
        boot_within = float(r[standard_endpoint].mean() - e[standard_endpoint].mean())
        bootstrap_rows.append(
            {
                "gene": gene,
                "bootstrap_iteration": index + 1,
                "observed_effect": boot_observed,
                "within_state_effect": boot_within,
                "composition_associated_effect": boot_observed - boot_within,
            }
        )
    bootstrap = pd.DataFrame(bootstrap_rows)

    def interval(column: str) -> tuple[float, float]:
        low, high = bootstrap[column].quantile([0.025, 0.975])
        return float(low), float(high)

    observed_ci = interval("observed_effect")
    within_ci = interval("within_state_effect")
    composition_ci = interval("composition_associated_effect")
    attenuation = (
        1 - abs(within) / abs(observed) if not np.isclose(observed, 0) else np.nan
    )
    composition_fraction = (
        composition / observed if not np.isclose(observed, 0) else np.nan
    )
    row = {
        "gene": gene,
        "reference_weights": ";".join(
            f"{row.cluster}={row.common_reference_weight:.8f}"
            for row in weights.itertuples(index=False)
        ),
        "observed_effect": observed,
        "observed_bootstrap_CI_low": observed_ci[0],
        "observed_bootstrap_CI_high": observed_ci[1],
        "within_state_effect": within,
        "within_state_bootstrap_CI_low": within_ci[0],
        "within_state_bootstrap_CI_high": within_ci[1],
        "composition_associated_effect": composition,
        "composition_associated_bootstrap_CI_low": composition_ci[0],
        "composition_associated_bootstrap_CI_high": composition_ci[1],
        "absolute_effect_attenuation_fraction": attenuation,
        "composition_associated_fraction_of_observed": composition_fraction,
    }
    return row, bootstrap


def mixed_patient_paired(patient_hgp: pd.DataFrame, endpoints: list[str]) -> pd.DataFrame:
    mixed = patient_hgp.loc[
        patient_hgp["canonical_patient"].isin(MIXED_PATIENTS),
        ["canonical_patient", "hgp", *endpoints],
    ].copy()
    rows: list[dict[str, object]] = []
    for patient in sorted(MIXED_PATIENTS):
        subset = mixed.loc[mixed["canonical_patient"].eq(patient)].set_index("hgp")
        for endpoint in endpoints:
            rows.append(
                {
                    "canonical_patient": patient,
                    "endpoint": endpoint,
                    "EHGP": float(subset.loc["EHGP", endpoint]),
                    "RHGP": float(subset.loc["RHGP", endpoint]),
                    "RHGP_minus_EHGP": float(
                        subset.loc["RHGP", endpoint] - subset.loc["EHGP", endpoint]
                    ),
                }
            )
    return pd.DataFrame(rows)


def fit_mixed_models(sample_endpoints: pd.DataFrame) -> pd.DataFrame:
    import statsmodels.formula.api as smf

    data = sample_endpoints.copy()
    data["rhgp"] = data["hgp"].eq("RHGP").astype(int)
    data["logit_nc3_fraction_neoplastic"] = np.log(
        (data["n_NC3_cells"] + 0.5)
        / (data["n_neoplastic_cells"] - data["n_NC3_cells"] + 0.5)
    )
    data["log1p_SLPI_NC3_mean"] = np.log1p(data["SLPI_NC3_mean"])
    specifications = [
        ("logit_NC3_fraction_random_patient_intercept", "logit_nc3_fraction_neoplastic"),
        ("log1p_within_NC3_SLPI_random_patient_intercept", "log1p_SLPI_NC3_mean"),
    ]
    rows: list[dict[str, object]] = []
    for model_name, endpoint in specifications:
        model = smf.mixedlm(
            f"{endpoint} ~ rhgp", data=data, groups=data["canonical_patient"]
        )
        fit = None
        errors: list[str] = []
        optimizer = ""
        for candidate in ["lbfgs", "powell", "nm"]:
            try:
                fit = model.fit(
                    reml=False, method=candidate, maxiter=1000, disp=False
                )
                optimizer = candidate
                break
            except (np.linalg.LinAlgError, ValueError) as error:
                errors.append(f"{candidate}:{type(error).__name__}:{error}")
        if fit is None:
            rows.append(
                {
                    "model": model_name,
                    "n_samples": len(data),
                    "n_patients": data["canonical_patient"].nunique(),
                    "optimizer": "failed",
                    "RHGP_coefficient": np.nan,
                    "standard_error": np.nan,
                    "wald_z": np.nan,
                    "wald_two_sided_p": np.nan,
                    "converged": False,
                    "random_intercept_variance": np.nan,
                    "fit_note": " | ".join(errors),
                }
            )
            continue
        coefficient = float(fit.params["rhgp"])
        standard_error = float(fit.bse["rhgp"])
        z_value = coefficient / standard_error
        rows.append(
            {
                "model": model_name,
                "n_samples": len(data),
                "n_patients": data["canonical_patient"].nunique(),
                "optimizer": optimizer,
                "RHGP_coefficient": coefficient,
                "standard_error": standard_error,
                "wald_z": z_value,
                "wald_two_sided_p": float(2 * norm.sf(abs(z_value))),
                "converged": bool(fit.converged),
                "random_intercept_variance": float(fit.cov_re.iloc[0, 0]),
                "fit_note": " | ".join(errors),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(INPUT)
    if (adata.n_obs, adata.n_vars) != (304_410, 113):
        raise ValueError(f"Unexpected dimensions: {adata.n_obs} x {adata.n_vars}")
    if adata.obs["Sample"].nunique() != 22:
        raise ValueError("Expected 22 post-QC samples")

    sample_endpoints = build_sample_endpoints(adata)
    patient_hgp = aggregate_patient_hgp(sample_endpoints)
    patient_hgp, weights = add_standardised_endpoints(patient_hgp)
    independent = patient_hgp.loc[
        ~patient_hgp["canonical_patient"].isin(MIXED_PATIENTS)
    ].copy()
    if independent["canonical_patient"].nunique() != 15:
        raise ValueError("Expected 15 independent patients")

    primary_endpoints = [
        "nc3_fraction_neoplastic",
        "SLPI_NC3_mean",
        "SLPI_standardised_mean",
        "SLPI_neoplastic_mean",
        "nc3_fraction_all_tumour",
    ]
    marker_endpoints = [
        endpoint
        for gene in GENES
        for endpoint in (
            f"{gene}_NC3_mean",
            f"{gene}_standardised_mean",
            f"{gene}_neoplastic_mean",
        )
    ]
    cluster_slpi_endpoints = [f"SLPI_{cluster}_mean" for cluster in CLUSTERS]
    all_endpoints = list(
        dict.fromkeys([*primary_endpoints, *marker_endpoints, *cluster_slpi_endpoints])
    )

    summaries = pd.DataFrame(
        [
            summarize_endpoint(independent, endpoint, seed_offset=index * 10_000)
            for index, endpoint in enumerate(all_endpoints)
        ]
    )
    loo = leave_one_patient_out(independent, all_endpoints)
    loo_summary = (
        loo.groupby("endpoint", as_index=False)
        .agg(
            loo_runs=("left_out_patient", "nunique"),
            loo_positive_runs=("positive_direction", "sum"),
            loo_positive_fraction=("positive_direction", "mean"),
            loo_min_hodges_lehmann=("hodges_lehmann_RHGP_minus_EHGP", "min"),
            loo_max_hodges_lehmann=("hodges_lehmann_RHGP_minus_EHGP", "max"),
        )
    )
    summaries = summaries.merge(loo_summary, on="endpoint", how="left")

    cluster_slpi = summaries.loc[
        summaries["endpoint"].isin(cluster_slpi_endpoints)
    ].copy()
    cluster_slpi["BH_FDR_across_four_clusters"] = multipletests(
        cluster_slpi["exact_rank_permutation_two_sided_p"], method="fdr_bh"
    )[1]

    decomposition_rows: list[dict[str, object]] = []
    bootstrap_frames: list[pd.DataFrame] = []
    for index, gene in enumerate(GENES):
        row, bootstrap = decompose_gene(
            independent, gene, weights, seed_offset=500_000 + index * 10_000
        )
        decomposition_rows.append(row)
        bootstrap_frames.append(bootstrap)
    decomposition = pd.DataFrame(decomposition_rows)
    decomposition_bootstrap = pd.concat(bootstrap_frames, ignore_index=True)

    paired_endpoints = list(
        dict.fromkeys(
            [
                "nc3_fraction_neoplastic",
                "nc3_fraction_all_tumour",
                *[
                    endpoint
                    for gene in GENES
                    for endpoint in (
                        f"{gene}_neoplastic_mean",
                        f"{gene}_NC3_mean",
                        f"{gene}_standardised_mean",
                    )
                ],
            ]
        )
    )
    paired = mixed_patient_paired(patient_hgp, paired_endpoints)
    mixed_models = fit_mixed_models(sample_endpoints)

    def summary_row(endpoint: str) -> pd.Series:
        return summaries.loc[summaries["endpoint"].eq(endpoint)].iloc[0]

    nc3 = summary_row("nc3_fraction_neoplastic")
    within_nc3 = summary_row("SLPI_NC3_mean")
    standardised = summary_row("SLPI_standardised_mean")
    raw = summary_row("SLPI_neoplastic_mean")
    slpi_decomposition = decomposition.loc[decomposition["gene"].eq("SLPI")].iloc[0]
    attenuation_gate = bool(
        standardised["mean_difference_RHGP_minus_EHGP"] <= 0
        or abs(standardised["mean_difference_RHGP_minus_EHGP"])
        <= 0.5 * abs(raw["mean_difference_RHGP_minus_EHGP"])
    )
    composition_supported = bool(
        nc3["mean_difference_RHGP_minus_EHGP"] > 0
        and raw["mean_difference_RHGP_minus_EHGP"] > 0
        and attenuation_gate
        and abs(slpi_decomposition["composition_associated_effect"])
        > abs(slpi_decomposition["within_state_effect"])
    )
    intrinsic_supported = bool(
        within_nc3["mean_difference_RHGP_minus_EHGP"] > 0
        and within_nc3["exact_rank_permutation_two_sided_p"] < 0.05
        and abs(standardised["mean_difference_RHGP_minus_EHGP"])
        >= 0.5 * abs(raw["mean_difference_RHGP_minus_EHGP"])
    )
    verdict = (
        "STATE_COMPOSITION_SUPPORTED"
        if composition_supported
        else "STATE_INTRINSIC_SUPPORTED"
        if intrinsic_supported
        else "UNRESOLVED_OR_MIXED"
    )
    decision = pd.DataFrame(
        [
            {
                "contract": str(CONTRACT.relative_to(ROOT)),
                "verdict": verdict,
                "state_composition_supported": composition_supported,
                "state_intrinsic_supported": intrinsic_supported,
                "NC3_fraction_direction_positive": bool(
                    nc3["mean_difference_RHGP_minus_EHGP"] > 0
                ),
                "NC3_fraction_exact_p_lt_0_05": bool(
                    nc3["exact_rank_permutation_two_sided_p"] < 0.05
                ),
                "raw_neoplastic_SLPI_direction_positive": bool(
                    raw["mean_difference_RHGP_minus_EHGP"] > 0
                ),
                "within_NC3_SLPI_direction_positive": bool(
                    within_nc3["mean_difference_RHGP_minus_EHGP"] > 0
                ),
                "within_NC3_SLPI_exact_p_lt_0_05": bool(
                    within_nc3["exact_rank_permutation_two_sided_p"] < 0.05
                ),
                "composition_standardised_attenuation_gate": attenuation_gate,
                "composition_component_larger_than_within_state_component": bool(
                    abs(slpi_decomposition["composition_associated_effect"])
                    > abs(slpi_decomposition["within_state_effect"])
                ),
                "core_HRC_projection_eligible": False,
                "main_text_role": "external_explanatory_boundary_not_novel_discovery",
            }
        ]
    )

    sample_endpoints.to_csv(OUT / "sample_state_endpoints.tsv", sep="\t", index=False)
    patient_hgp.to_csv(OUT / "patient_state_endpoints.tsv", sep="\t", index=False)
    weights.to_csv(OUT / "common_reference_cluster_weights.tsv", sep="\t", index=False)
    summaries.to_csv(OUT / "endpoint_summary.tsv", sep="\t", index=False)
    cluster_slpi.to_csv(OUT / "within_cluster_SLPI_effects.tsv", sep="\t", index=False)
    loo.to_csv(OUT / "leave_one_patient_out.tsv", sep="\t", index=False)
    decomposition.to_csv(OUT / "marker_decomposition.tsv", sep="\t", index=False)
    decomposition_bootstrap.to_csv(
        OUT / "marker_decomposition_bootstrap_replicates.tsv", sep="\t", index=False
    )
    paired.to_csv(OUT / "mixed_patient_paired_sensitivity.tsv", sep="\t", index=False)
    mixed_models.to_csv(OUT / "sample_level_mixed_models.tsv", sep="\t", index=False)
    decision.to_csv(OUT / "validation_decision.tsv", sep="\t", index=False)

    run_info = {
        "input": str(INPUT.relative_to(ROOT)),
        "input_sha256": sha256(INPUT),
        "contract": str(CONTRACT.relative_to(ROOT)),
        "doi": "10.1038/s41416-026-03567-y",
        "source_code_commit": SOURCE_CODE_COMMIT,
        "python": platform.python_version(),
        "anndata": version("anndata"),
        "numpy": version("numpy"),
        "pandas": version("pandas"),
        "scipy": version("scipy"),
        "statsmodels": version("statsmodels"),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_samples": int(sample_endpoints["Sample"].nunique()),
        "n_canonical_patients": int(patient_hgp["canonical_patient"].nunique()),
        "n_independent_patients": int(independent["canonical_patient"].nunique()),
        "bootstrap_seed": BOOTSTRAP_SEED,
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
    }
    (OUT / "run_info.json").write_text(
        json.dumps(run_info, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(decision.to_string(index=False))
    print(
        summaries.loc[summaries["endpoint"].isin(primary_endpoints)].to_string(
            index=False
        )
    )
    print(decomposition.to_string(index=False))


if __name__ == "__main__":
    main()
