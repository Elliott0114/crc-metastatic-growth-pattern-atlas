#!/usr/bin/env python3
"""
Analysis: within-cell covariance of ATAC-informed regulons and junction genes.
Date: 2026-09-01
Random seed: 42
Python and key package versions are written to the provenance JSON.
"""

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
from scipy import sparse, stats


SEED = 42
BOOTSTRAP_ITERATIONS = 10_000
MIN_CORRELATION_CELLS = 5
SUMMARY_CELL_THRESHOLDS = (10, 20, 30)
STATES = ("REC", "Hypoxia", "UPR", "iREC")
METHODS = ("matched_residual", "matched_unadjusted", "gene_z_residual")

ROOT = Path(__file__).resolve().parents[1]
INPUT = (
    ROOT
    / "data_sources"
    / "Ogden_2025_CRLM_multiome"
    / "CRCLM_multiome_GEX_decontaminated.h5ad"
)
SPEC = ROOT / "metadata" / "ogden_atac_regulatory_junction_bridge_spec_2026-09-01.md"
OUT = (
    ROOT
    / "analysis_results"
    / "deep_biology_upgrade_2026-08-31"
    / "phase2_mechanistic_specificity"
    / "ogden_regulatory_junction_bridge"
)
LOCKED = OUT / "ogden_regulatory_junction_locked_genes.tsv"
REGISTRY = OUT / "ogden_regulatory_junction_pair_registry.tsv"
LOCK_PROVENANCE = OUT / "ogden_regulatory_junction_lock_provenance.json"

np.random.seed(SEED)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_seed(label: str) -> int:
    digest = hashlib.sha256(label.encode("utf-8")).digest()
    return SEED + int.from_bytes(digest[:4], "little")


def bh_adjust(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    finite = values.notna()
    p = values.loc[finite].to_numpy(dtype=float)
    if not len(p):
        return result
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    result.loc[finite] = restored
    return result


def residualize(values: np.ndarray, covariates: np.ndarray) -> np.ndarray:
    columns = [np.ones(len(values), dtype=float)]
    for column in covariates.T:
        sd = float(np.std(column, ddof=0))
        if np.isfinite(sd) and sd > 0:
            columns.append((column - float(np.mean(column))) / sd)
    design = np.column_stack(columns)
    beta = np.linalg.lstsq(design, values, rcond=None)[0]
    return values - design @ beta


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < MIN_CORRELATION_CELLS:
        return np.nan
    x_keep = x[keep]
    y_keep = y[keep]
    if np.ptp(x_keep) == 0 or np.ptp(y_keep) == 0:
        return np.nan
    return float(stats.spearmanr(x_keep, y_keep).statistic)


def fisher_z(rho: np.ndarray | float) -> np.ndarray | float:
    return np.arctanh(np.clip(rho, -0.999999, 0.999999))


def exact_sign_flip_p(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan
    observed = abs(float(np.mean(values)))
    exceed = 0
    total = 2 ** len(values)
    for signs in itertools.product((-1.0, 1.0), repeat=len(values)):
        permuted = abs(float(np.mean(values * np.asarray(signs))))
        if permuted >= observed - 1e-12:
            exceed += 1
    return exceed / total


def bootstrap_meta_rho(values: np.ndarray, label: str) -> tuple[float, float]:
    z_values = np.asarray(fisher_z(values), dtype=float)
    rng = np.random.default_rng(stable_seed(label))
    draws = rng.choice(
        z_values,
        size=(BOOTSTRAP_ITERATIONS, len(z_values)),
        replace=True,
    )
    estimates = np.tanh(draws.mean(axis=1))
    return tuple(np.quantile(estimates, (0.025, 0.975)).tolist())


def bootstrap_mean(values: np.ndarray, label: str) -> tuple[float, float]:
    rng = np.random.default_rng(stable_seed(label))
    draws = rng.choice(
        values,
        size=(BOOTSTRAP_ITERATIONS, len(values)),
        replace=True,
    ).mean(axis=1)
    return tuple(np.quantile(draws, (0.025, 0.975)).tolist())


def sparse_row_mean(matrix: sparse.csr_matrix, positions: list[int]) -> np.ndarray:
    return np.asarray(matrix[:, positions].mean(axis=1)).ravel()


def gene_balanced_z_score(
    matrix: sparse.csr_matrix, positions: list[int]
) -> tuple[np.ndarray, int]:
    selected = matrix[:, positions].tocsr()
    gene_mean = np.asarray(selected.mean(axis=0)).ravel()
    gene_second_moment = np.asarray(selected.power(2).mean(axis=0)).ravel()
    gene_variance = np.maximum(gene_second_moment - gene_mean**2, 0)
    valid = gene_variance > 0
    if valid.sum() < 5:
        raise ValueError("Fewer than five non-constant genes in gene-z score")
    selected = selected[:, valid]
    gene_mean = gene_mean[valid]
    gene_sd = np.sqrt(gene_variance[valid])
    weighted = selected.multiply(1.0 / gene_sd)
    score = np.asarray(weighted.mean(axis=1)).ravel() - np.mean(gene_mean / gene_sd)
    return score, int(valid.sum())


def patient_equal_summary(
    correlations: pd.DataFrame, threshold: int
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    eligible = correlations.loc[correlations["n_cells"] >= threshold].copy()
    group_columns = [
        "pair_id",
        "regulon_id",
        "candidate_class",
        "method",
        "state",
    ]
    for keys, group in eligible.groupby(group_columns, sort=True, observed=True):
        values = group["spearman_rho"].dropna().to_numpy(dtype=float)
        if not len(values):
            continue
        z_values = np.asarray(fisher_z(values), dtype=float)
        ci_low, ci_high = bootstrap_meta_rho(
            values, f"summary:{threshold}:{':'.join(map(str, keys))}"
        )
        rows.append(
            {
                **dict(zip(group_columns, keys)),
                "minimum_cells": threshold,
                "n_patients": len(values),
                "n_positive": int((values > 0).sum()),
                "n_negative": int((values < 0).sum()),
                "n_zero": int((values == 0).sum()),
                "median_rho": float(np.median(values)),
                "minimum_rho": float(np.min(values)),
                "maximum_rho": float(np.max(values)),
                "patient_equal_meta_rho": float(np.tanh(np.mean(z_values))),
                "bootstrap_ci_low": ci_low,
                "bootstrap_ci_high": ci_high,
                "exact_sign_flip_p": exact_sign_flip_p(z_values),
            }
        )
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    result["sign_flip_fdr_across_regulons"] = result.groupby(
        ["state", "method", "minimum_cells"], observed=True
    )["exact_sign_flip_p"].transform(bh_adjust)
    return result


def rec_comparator_summary(
    correlations: pd.DataFrame, threshold: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    eligible = correlations.loc[correlations["n_cells"] >= threshold].copy()
    summary_rows: list[dict[str, object]] = []
    patient_rows: list[dict[str, object]] = []
    for keys, group in eligible.groupby(
        ["pair_id", "regulon_id", "candidate_class", "method"],
        sort=True,
        observed=True,
    ):
        pair_id, regulon_id, candidate_class, method = keys
        rec = group.loc[group["state"].eq("REC"), ["patient", "spearman_rho"]]
        rec = rec.rename(columns={"spearman_rho": "rec_rho"})
        for comparator in ("Hypoxia", "UPR", "iREC"):
            other = group.loc[
                group["state"].eq(comparator), ["patient", "spearman_rho"]
            ].rename(columns={"spearman_rho": "comparator_rho"})
            paired = rec.merge(other, on="patient", how="inner").dropna()
            if paired.empty:
                continue
            paired["rho_difference"] = paired["rec_rho"] - paired["comparator_rho"]
            paired["fisher_z_difference"] = fisher_z(
                paired["rec_rho"].to_numpy()
            ) - fisher_z(paired["comparator_rho"].to_numpy())
            for row in paired.itertuples(index=False):
                patient_rows.append(
                    {
                        "pair_id": pair_id,
                        "regulon_id": regulon_id,
                        "candidate_class": candidate_class,
                        "method": method,
                        "minimum_cells": threshold,
                        "comparator_state": comparator,
                        "patient": row.patient,
                        "rec_rho": row.rec_rho,
                        "comparator_rho": row.comparator_rho,
                        "rho_difference": row.rho_difference,
                        "fisher_z_difference": row.fisher_z_difference,
                    }
                )
            rho_diff = paired["rho_difference"].to_numpy(dtype=float)
            z_diff = paired["fisher_z_difference"].to_numpy(dtype=float)
            ci_low, ci_high = bootstrap_mean(
                rho_diff,
                f"comparison:{threshold}:{pair_id}:{method}:{comparator}",
            )
            summary_rows.append(
                {
                    "pair_id": pair_id,
                    "regulon_id": regulon_id,
                    "candidate_class": candidate_class,
                    "method": method,
                    "minimum_cells": threshold,
                    "comparator_state": comparator,
                    "n_paired_patients": len(paired),
                    "n_rec_greater": int((rho_diff > 0).sum()),
                    "n_rec_lower": int((rho_diff < 0).sum()),
                    "mean_rec_rho": float(paired["rec_rho"].mean()),
                    "mean_comparator_rho": float(paired["comparator_rho"].mean()),
                    "mean_rho_difference": float(np.mean(rho_diff)),
                    "median_rho_difference": float(np.median(rho_diff)),
                    "bootstrap_ci_low": ci_low,
                    "bootstrap_ci_high": ci_high,
                    "mean_fisher_z_difference": float(np.mean(z_diff)),
                    "exact_paired_sign_flip_p": exact_sign_flip_p(z_diff),
                }
            )
    summary = pd.DataFrame(summary_rows)
    if not summary.empty:
        summary["paired_sign_flip_fdr_across_regulons"] = summary.groupby(
            ["method", "minimum_cells", "comparator_state"], observed=True
        )["exact_paired_sign_flip_p"].transform(bh_adjust)
    return summary, pd.DataFrame(patient_rows)


def main() -> None:
    required = (INPUT, SPEC, LOCKED, REGISTRY, LOCK_PROVENANCE)
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing inputs: {missing}")

    lock_provenance = json.loads(LOCK_PROVENANCE.read_text(encoding="utf-8"))
    if lock_provenance["status"] != "frozen_before_regulon_junction_covariance_outcomes":
        raise ValueError("Regulon-junction lock has unexpected status")
    if lock_provenance["specification_sha256"] != sha256(SPEC):
        raise ValueError("Specification changed after locking")
    if lock_provenance["locked_genes_sha256"] != sha256(LOCKED):
        raise ValueError("Locked gene table hash mismatch")
    if lock_provenance["pair_registry_sha256"] != sha256(REGISTRY):
        raise ValueError("Pair registry hash mismatch")

    locked = pd.read_csv(LOCKED, sep="\t", dtype={"gene": str})
    registry = pd.read_csv(REGISTRY, sep="\t")
    pair_metadata = registry[
        ["pair_id", "regulon_id", "candidate_class"]
    ].drop_duplicates().sort_values("pair_id")

    adata = ad.read_h5ad(INPUT)
    required_obs = {"Patient", "Therapy", "Cell_type", "Cell_subtype"}
    if missing_obs := sorted(required_obs.difference(adata.obs.columns)):
        raise ValueError(f"Missing observation fields: {missing_obs}")
    if adata.var_names.has_duplicates or not adata.obs_names.is_unique:
        raise ValueError("Gene and cell identifiers must be unique")

    epithelial_mask = (
        adata.obs["Cell_type"].astype("string").eq("Epithelial").to_numpy()
    )
    epithelial_obs = adata.obs.loc[
        epithelial_mask, ["Patient", "Therapy", "Cell_subtype"]
    ].copy()
    epithelial_obs.insert(0, "cell_id", adata.obs_names[epithelial_mask].astype(str))
    epithelial_obs = epithelial_obs.reset_index(drop=True)
    for column in ("Patient", "Therapy", "Cell_subtype"):
        epithelial_obs[column] = epithelial_obs[column].astype("string")

    matrix = adata.X[epithelial_mask]
    if not sparse.issparse(matrix):
        matrix = sparse.csr_matrix(np.asarray(matrix, dtype=np.float64))
    else:
        matrix = sparse.csr_matrix(matrix, dtype=np.float64)
    if matrix.data.size and (
        not np.isfinite(matrix.data).all() or (matrix.data < 0).any()
    ):
        raise ValueError("Expression matrix must contain finite non-negative values")
    library_size = np.asarray(matrix.sum(axis=1)).ravel()
    detected_genes = np.asarray(matrix.getnnz(axis=1)).ravel()
    if np.any(library_size <= 0):
        raise ValueError("One or more epithelial cells has zero expression total")

    scale = 10_000.0 / library_size
    lognorm = matrix.multiply(scale[:, None]).tocsr()
    np.log1p(lognorm.data, out=lognorm.data)
    genes = pd.Index(adata.var_names.astype(str).str.upper())
    gene_to_position = {gene: index for index, gene in enumerate(genes)}

    selected_mask = epithelial_obs["Cell_subtype"].isin(STATES).to_numpy()
    selected_obs = epithelial_obs.loc[selected_mask].reset_index(drop=True)
    selected_matrix = matrix[selected_mask].tocsr()
    selected_lognorm = lognorm[selected_mask].tocsr()
    technical_covariates = np.column_stack(
        [
            np.log1p(library_size[selected_mask]),
            np.log1p(detected_genes[selected_mask]),
        ]
    )
    group_indices = selected_obs.groupby(
        ["Patient", "Cell_subtype"], sort=True, observed=True
    ).indices

    correlation_rows: list[dict[str, object]] = []
    qc_rows: list[dict[str, object]] = []
    for meta in pair_metadata.itertuples(index=False):
        pair_lock = locked.loc[locked["pair_id"].eq(meta.pair_id)].copy()
        side_scores: dict[str, dict[str, np.ndarray | int | str]] = {}
        side_target_sets: dict[str, set[str]] = {}
        side_control_sets: dict[str, set[str]] = {}
        for side in ("a", "b"):
            side_lock = pair_lock.loc[pair_lock["side"].eq(side)]
            program_values = side_lock["program_id"].drop_duplicates().tolist()
            if len(program_values) != 1:
                raise ValueError(f"Ambiguous program for {meta.pair_id} side {side}")
            target_genes = side_lock.loc[
                side_lock["gene_role"].eq("target"), "gene"
            ].tolist()
            control_genes = side_lock.loc[
                side_lock["gene_role"].eq("control"), "gene"
            ].tolist()
            missing_genes = sorted(
                set(target_genes + control_genes).difference(gene_to_position)
            )
            if missing_genes:
                raise ValueError(f"Locked genes absent from matrix: {missing_genes[:5]}")
            target_positions = [gene_to_position[gene] for gene in target_genes]
            control_positions = [gene_to_position[gene] for gene in control_genes]
            target_mean = sparse_row_mean(selected_lognorm, target_positions)
            control_mean = sparse_row_mean(selected_lognorm, control_positions)
            matched_score = target_mean - control_mean
            gene_z_score, z_gene_count = gene_balanced_z_score(
                selected_lognorm, target_positions
            )
            target_detection = np.asarray(
                (selected_matrix[:, target_positions] > 0).mean(axis=1)
            ).ravel()
            side_scores[side] = {
                "program_id": program_values[0],
                "matched": matched_score,
                "gene_z": gene_z_score,
                "detection": target_detection,
                "target_count": len(target_positions),
                "control_count": len(control_positions),
                "z_gene_count": z_gene_count,
            }
            side_target_sets[side] = set(target_genes)
            side_control_sets[side] = set(control_genes)

        if side_target_sets["a"].intersection(side_target_sets["b"]):
            raise ValueError(f"Target overlap remains in {meta.pair_id}")
        if side_control_sets["a"].intersection(side_control_sets["b"]):
            raise ValueError(f"Control overlap remains in {meta.pair_id}")

        residuals: dict[str, np.ndarray] = {}
        for side in ("a", "b"):
            for score_type in ("matched", "gene_z"):
                values = np.asarray(side_scores[side][score_type], dtype=float)
                result = np.full(len(values), np.nan, dtype=float)
                for indices in group_indices.values():
                    index = np.asarray(indices, dtype=int)
                    if len(index) >= MIN_CORRELATION_CELLS:
                        result[index] = residualize(
                            values[index], technical_covariates[index]
                        )
                residuals[f"{side}_{score_type}"] = result

        method_arrays = {
            "matched_residual": (
                residuals["a_matched"],
                residuals["b_matched"],
            ),
            "matched_unadjusted": (
                np.asarray(side_scores["a"]["matched"], dtype=float),
                np.asarray(side_scores["b"]["matched"], dtype=float),
            ),
            "gene_z_residual": (
                residuals["a_gene_z"],
                residuals["b_gene_z"],
            ),
        }

        for (patient, state), indices in group_indices.items():
            index = np.asarray(indices, dtype=int)
            therapy_values = selected_obs.loc[index, "Therapy"].dropna().unique()
            if len(therapy_values) != 1:
                raise ValueError(f"Inconsistent therapy for {patient} {state}")
            for method, (values_a, values_b) in method_arrays.items():
                correlation_rows.append(
                    {
                        "pair_id": meta.pair_id,
                        "regulon_id": meta.regulon_id,
                        "candidate_class": meta.candidate_class,
                        "method": method,
                        "patient": str(patient),
                        "therapy": str(therapy_values[0]),
                        "state": str(state),
                        "n_cells": len(index),
                        "spearman_rho": safe_spearman(
                            values_a[index], values_b[index]
                        ),
                    }
                )

        for state in STATES:
            state_index = np.flatnonzero(
                selected_obs["Cell_subtype"].eq(state).to_numpy()
            )
            if not len(state_index):
                continue
            for side in ("a", "b"):
                matched = np.asarray(side_scores[side]["matched"], dtype=float)
                matched_residual = residuals[f"{side}_matched"]
                finite = np.isfinite(matched_residual[state_index])
                residual_index = state_index[finite]
                qc_rows.append(
                    {
                        "pair_id": meta.pair_id,
                        "regulon_id": meta.regulon_id,
                        "candidate_class": meta.candidate_class,
                        "state": state,
                        "side": side,
                        "program_id": side_scores[side]["program_id"],
                        "n_cells": len(state_index),
                        "target_gene_count": side_scores[side]["target_count"],
                        "matched_control_gene_count": side_scores[side]["control_count"],
                        "gene_z_nonconstant_gene_count": side_scores[side]["z_gene_count"],
                        "median_target_detection_fraction": float(
                            np.median(
                                np.asarray(side_scores[side]["detection"])[state_index]
                            )
                        ),
                        "matched_score_sd": float(np.std(matched[state_index])),
                        "matched_vs_log_library_spearman": safe_spearman(
                            matched[state_index], technical_covariates[state_index, 0]
                        ),
                        "residual_vs_log_library_spearman": safe_spearman(
                            matched_residual[residual_index],
                            technical_covariates[residual_index, 0],
                        ),
                        "residual_vs_log_detected_spearman": safe_spearman(
                            matched_residual[residual_index],
                            technical_covariates[residual_index, 1],
                        ),
                    }
                )
        print(f"Completed regulon-junction pair {meta.regulon_id}", flush=True)

    correlations = pd.DataFrame(correlation_rows).sort_values(
        ["regulon_id", "method", "state", "patient"]
    )
    summaries = pd.concat(
        [
            patient_equal_summary(correlations, threshold)
            for threshold in SUMMARY_CELL_THRESHOLDS
        ],
        ignore_index=True,
    ).sort_values(["minimum_cells", "method", "state", "regulon_id"])
    comparison_parts = [
        rec_comparator_summary(correlations, threshold)
        for threshold in SUMMARY_CELL_THRESHOLDS
    ]
    comparisons = pd.concat(
        [part[0] for part in comparison_parts], ignore_index=True
    ).sort_values(
        ["minimum_cells", "method", "comparator_state", "regulon_id"]
    )
    comparison_patients = pd.concat(
        [part[1] for part in comparison_parts], ignore_index=True
    ).sort_values(
        [
            "minimum_cells",
            "method",
            "comparator_state",
            "regulon_id",
            "patient",
        ]
    )
    score_qc = pd.DataFrame(qc_rows).sort_values(
        ["regulon_id", "state", "side"]
    )

    primary = summaries.loc[
        summaries["minimum_cells"].eq(20)
        & summaries["method"].eq("matched_residual")
        & summaries["state"].eq("REC")
        & summaries["candidate_class"].isin(
            [
                "primary_AP1",
                "inflammatory_comparator",
                "mechanotransduction_comparator",
                "lineage_WNT_comparator",
            ]
        )
    ].sort_values("patient_equal_meta_rho", ascending=False)

    outputs = {
        "ogden_regulatory_junction_patient_state_correlations.tsv": correlations,
        "ogden_regulatory_junction_patient_equal_summaries.tsv": summaries,
        "ogden_regulatory_junction_rec_vs_comparator.tsv": comparisons,
        "ogden_regulatory_junction_rec_vs_comparator_patient.tsv": comparison_patients,
        "ogden_regulatory_junction_score_qc.tsv": score_qc,
        "ogden_regulatory_junction_primary_summary.tsv": primary,
    }
    for filename, frame in outputs.items():
        frame.to_csv(OUT / filename, sep="\t", index=False)

    provenance = {
        "analysis_date": "2026-09-01",
        "python_version": platform.python_version(),
        "anndata_version": version("anndata"),
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "scipy_version": version("scipy"),
        "seed": SEED,
        "bootstrap_iterations": BOOTSTRAP_ITERATIONS,
        "minimum_cells_sensitivities": list(SUMMARY_CELL_THRESHOLDS),
        "states": list(STATES),
        "normalization": "log1p counts per 10,000",
        "primary_score": "target mean minus disjoint expression-bin-matched control mean",
        "technical_residualization": "within patient-state OLS on log1p library size and detected genes",
        "input_sha256": sha256(INPUT),
        "specification_sha256": sha256(SPEC),
        "locked_genes_sha256": sha256(LOCKED),
        "pair_registry_sha256": sha256(REGISTRY),
    }
    with (OUT / "ogden_regulatory_junction_coexistence_provenance.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(provenance, handle, ensure_ascii=False, indent=2)

    manifest_path = OUT / "_analysis_outputs.md"
    existing = manifest_path.read_text(encoding="utf-8")
    addition = """

## Within-cell covariance tables

- `ogden_regulatory_junction_patient_state_correlations.tsv` -- Patient-state correlations.
- `ogden_regulatory_junction_patient_equal_summaries.tsv` -- Patient-equal estimates, intervals and exact tests.
- `ogden_regulatory_junction_rec_vs_comparator.tsv` -- Paired REC-versus-state correlation differences.
- `ogden_regulatory_junction_rec_vs_comparator_patient.tsv` -- Contributing paired patients.
- `ogden_regulatory_junction_score_qc.tsv` -- Coverage and technical-covariance checks.
- `ogden_regulatory_junction_primary_summary.tsv` -- Reader-facing candidate and control regulons.
- `ogden_regulatory_junction_coexistence_provenance.json` -- Analysis provenance.
"""
    manifest_path.write_text(existing.rstrip() + addition, encoding="utf-8")

    print("\nPrimary REC regulon-junction covariance")
    print(primary.to_string(index=False))


if __name__ == "__main__":
    main()
