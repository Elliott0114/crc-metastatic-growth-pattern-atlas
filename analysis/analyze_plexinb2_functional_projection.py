#!/usr/bin/env python3
"""Project frozen REC-state programmes into the GSE267981 Plexin B2 experiment.

The source experiment has one 10x library per condition.  All treatment
contrasts and equal-cell resampling intervals are therefore descriptive and
must not be interpreted as replicated intervention inference.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import mannwhitneyu


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data_sources" / "Borrelli_2024_PlexinB2_KLF4" / "GSE267981"
PHASE2 = (
    ROOT
    / "analysis_results"
    / "deep_biology_upgrade_2026-08-31"
    / "phase2_mechanistic_specificity"
)
PHASE3 = (
    ROOT
    / "analysis_results"
    / "deep_biology_upgrade_2026-08-31"
    / "phase3_external_regulatory_projection"
)
OUTPUT = (
    ROOT
    / "analysis_results"
    / "deep_biology_upgrade_2026-08-31"
    / "phase4_plexinb2_functional_projection"
)

FILES = {
    "rPlexinB2": DATA / "GSM8282884_rPlxnb2_filtered_feature_bc_matrix.h5",
    "vehicle": DATA / "GSM8282885_ctrl_filtered_feature_bc_matrix.h5",
}
CONDITIONS = ("vehicle", "rPlexinB2")
NAMED_GENES = (
    "KLF4",
    "ELF3",
    "GRHL2",
    "EPCAM",
    "ZEB1",
    "LGR5",
    "MKI67",
    "CDH17",
    "TJP1",
    "CLDN2",
    "CLDN3",
    "CLDN4",
    "CLDN7",
    "CRB3",
    "F11R",
)
BASE_SEED = 42
N_RAREFACTION = 500


def decode(values: np.ndarray) -> np.ndarray:
    return np.asarray(
        [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values]
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_seed(*parts: str) -> int:
    payload = "|".join((str(BASE_SEED), *parts)).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:4], "little")


def read_10x_h5(path: Path) -> tuple[sparse.csr_matrix, np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        shape = tuple(int(value) for value in matrix["shape"][:])
        genes_by_cells = sparse.csc_matrix(
            (
                matrix["data"][:],
                matrix["indices"][:],
                matrix["indptr"][:],
            ),
            shape=shape,
        )
        features = matrix["features"]
        names = decode(features["name"][:])
        barcodes = decode(matrix["barcodes"][:])
        if "feature_type" in features:
            feature_type = decode(features["feature_type"][:])
            keep = feature_type == "Gene Expression"
            genes_by_cells = genes_by_cells[keep, :]
            names = names[keep]
    return genes_by_cells.T.tocsr(), names, barcodes


def load_programmes() -> dict[str, list[str]]:
    phase2 = pd.read_csv(PHASE2 / "frozen_gene_sets_long.tsv", sep="\t")
    phase3 = pd.read_csv(
        PHASE3 / "frozen_external_regulatory_gene_sets_long.tsv", sep="\t"
    )
    coexistence = pd.read_csv(
        PHASE2 / "ogden_rec_coexistence_locked_genes.tsv", sep="\t"
    )

    def genes(frame: pd.DataFrame, set_id: str) -> list[str]:
        values = frame.loc[frame["set_id"] == set_id, "gene"].astype(str).str.upper()
        return list(dict.fromkeys(values))

    programmes = {
        "REC_PROGRAM_TOP50": genes(phase2, "FROZEN_REC_TOP50"),
        "REACTOME_TIGHT_JUNCTION_INTERACTIONS": genes(
            phase3, "REACTOME_TIGHT_JUNCTION_INTERACTIONS"
        ),
        "CORE_HRC": genes(phase2, "CANELLAS_CORE_HRC"),
        "PARTIAL_EMT": genes(phase2, "PUBLISHED_PEMT"),
        "HALLMARK_E2F_TARGETS": genes(phase2, "HALLMARK_E2F_TARGETS"),
        "HALLMARK_G2M_CHECKPOINT": genes(phase2, "HALLMARK_G2M_CHECKPOINT"),
    }
    actin = coexistence.loc[
        (coexistence["program_id"] == "ACTIN_TURNOVER_UNION")
        & (coexistence["gene_role"] == "target"),
        "gene",
    ].astype(str).str.upper()
    programmes["ACTIN_TURNOVER"] = list(dict.fromkeys(actin))
    programmes["JUNB_TIGHT_JUNCTION_TARGETS_5_POSTHOC"] = genes(
        phase3, "JUNB_TIGHT_JUNCTION_TARGETS_5"
    )

    rec = programmes["REC_PROGRAM_TOP50"]
    junction = programmes["REACTOME_TIGHT_JUNCTION_INTERACTIONS"]
    overlap = set(rec).intersection(junction)
    programmes["REC_PROGRAM_TOP50_DISJOINT_JUNCTION"] = [
        gene for gene in rec if gene not in overlap
    ]
    programmes["TIGHT_JUNCTION_DISJOINT_REC"] = [
        gene for gene in junction if gene not in overlap
    ]
    return programmes


def combine_libraries() -> tuple[sparse.csr_matrix, np.ndarray, pd.DataFrame]:
    matrices: list[sparse.csr_matrix] = []
    metadata: list[pd.DataFrame] = []
    reference_names: np.ndarray | None = None
    for condition in CONDITIONS:
        matrix, names, barcodes = read_10x_h5(FILES[condition])
        if reference_names is None:
            reference_names = names
        elif not np.array_equal(reference_names, names):
            raise ValueError("The two 10x libraries do not have identical feature order")
        matrices.append(matrix)
        metadata.append(
            pd.DataFrame(
                {
                    "condition": condition,
                    "barcode": [f"{condition}:{barcode}" for barcode in barcodes],
                }
            )
        )
    if reference_names is None:
        raise RuntimeError("No 10x features were loaded")
    counts = sparse.vstack(matrices, format="csr")
    cells = pd.concat(metadata, ignore_index=True)

    detected_cells = np.asarray((counts > 0).sum(axis=0)).ravel()
    keep_genes = detected_cells >= 3
    return counts[:, keep_genes].tocsr(), reference_names[keep_genes], cells


def qc_table(counts: sparse.csr_matrix, genes: np.ndarray, cells: pd.DataFrame) -> pd.DataFrame:
    library_size = np.asarray(counts.sum(axis=1)).ravel()
    n_features = np.diff(counts.indptr)
    mito = np.char.startswith(np.char.lower(genes.astype(str)), "mt-")
    mito_counts = (
        np.asarray(counts[:, mito].sum(axis=1)).ravel() if mito.any() else np.zeros(counts.shape[0])
    )
    pct_mito = np.divide(
        mito_counts * 100.0,
        library_size,
        out=np.zeros_like(mito_counts, dtype=float),
        where=library_size > 0,
    )
    result = cells.copy()
    result["library_size"] = library_size
    result["n_features"] = n_features
    result["pct_mito"] = pct_mito
    result["primary_qc"] = (
        (n_features > 100) & (n_features < 2500) & (pct_mito < 25.0)
    )
    result["strict_qc"] = (
        (n_features > 200) & (n_features < 2500) & (pct_mito < 15.0)
    )
    return result


def log_normalize(counts: sparse.csr_matrix) -> sparse.csr_matrix:
    library_size = np.asarray(counts.sum(axis=1)).ravel()
    scale = np.divide(
        10000.0,
        library_size,
        out=np.zeros_like(library_size, dtype=float),
        where=library_size > 0,
    )
    normalized = sparse.diags(scale).dot(counts).tocsr()
    normalized.data = np.log1p(normalized.data)
    return normalized


def symbol_index(genes: np.ndarray) -> dict[str, list[int]]:
    mapping: dict[str, list[int]] = {}
    for index, gene in enumerate(genes):
        mapping.setdefault(str(gene).upper(), []).append(index)
    return mapping


def dense_gene_values(
    expression: sparse.csr_matrix,
    mapping: dict[str, list[int]],
    symbol: str,
) -> np.ndarray:
    indices = mapping.get(symbol.upper(), [])
    if not indices:
        raise KeyError(symbol)
    values = np.asarray(expression[:, indices].sum(axis=1)).ravel()
    return values


def programme_score(
    expression: sparse.csr_matrix,
    mapping: dict[str, list[int]],
    declared_genes: list[str],
) -> tuple[np.ndarray, list[str], list[str], pd.DataFrame]:
    measurable = [gene for gene in declared_genes if gene in mapping]
    missing = [gene for gene in declared_genes if gene not in mapping]
    if not measurable:
        raise ValueError("Programme has no measurable genes")
    columns = np.column_stack(
        [dense_gene_values(expression, mapping, gene) for gene in measurable]
    )
    means = columns.mean(axis=0)
    standard_deviations = columns.std(axis=0, ddof=0)
    variable = standard_deviations > 0
    retained = [gene for gene, keep in zip(measurable, variable, strict=True) if keep]
    zero_variance = [gene for gene, keep in zip(measurable, variable, strict=True) if not keep]
    if not retained:
        raise ValueError("Programme has no variable measurable genes")
    z_values = (columns[:, variable] - means[variable]) / standard_deviations[variable]
    score = z_values.mean(axis=1)
    gene_effects = pd.DataFrame(
        {
            "gene": retained,
            "expression_mean": means[variable],
            "expression_sd": standard_deviations[variable],
        }
    )
    return score, retained, missing + zero_variance, gene_effects


def effect_summary(values: np.ndarray, condition: np.ndarray) -> dict[str, float | int]:
    treated = values[condition == "rPlexinB2"]
    vehicle = values[condition == "vehicle"]
    auc = mannwhitneyu(treated, vehicle, alternative="two-sided").statistic / (
        treated.size * vehicle.size
    )
    return {
        "n_treated_cells": int(treated.size),
        "n_vehicle_cells": int(vehicle.size),
        "mean_treated": float(treated.mean()),
        "mean_vehicle": float(vehicle.mean()),
        "mean_difference": float(treated.mean() - vehicle.mean()),
        "median_treated": float(np.median(treated)),
        "median_vehicle": float(np.median(vehicle)),
        "median_difference": float(np.median(treated) - np.median(vehicle)),
        "auc_treated_greater_vehicle": float(auc),
    }


def rarefy(
    values: np.ndarray,
    condition: np.ndarray,
    qc_name: str,
    endpoint_id: str,
) -> pd.DataFrame:
    treated = values[condition == "rPlexinB2"]
    vehicle = values[condition == "vehicle"]
    n_cells = min(treated.size, vehicle.size)
    rng = np.random.default_rng(stable_seed(qc_name, endpoint_id, "rarefaction"))

    def sampled_summaries(group: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if group.size == n_cells:
            return (
                np.full(N_RAREFACTION, group.mean(), dtype=float),
                np.full(N_RAREFACTION, np.median(group), dtype=float),
            )
        indices = np.vstack(
            [
                rng.choice(group.size, size=n_cells, replace=False)
                for _ in range(N_RAREFACTION)
            ]
        )
        sampled = group[indices]
        return sampled.mean(axis=1), np.median(sampled, axis=1)

    treated_mean, treated_median = sampled_summaries(treated)
    vehicle_mean, vehicle_median = sampled_summaries(vehicle)
    return pd.DataFrame(
        {
            "qc_definition": qc_name,
            "endpoint_id": endpoint_id,
            "iteration": np.arange(1, N_RAREFACTION + 1),
            "cells_per_condition": n_cells,
            "mean_difference": treated_mean - vehicle_mean,
            "median_difference": treated_median - vehicle_median,
        }
    )


def analyse_qc(
    counts: sparse.csr_matrix,
    genes: np.ndarray,
    qc: pd.DataFrame,
    qc_column: str,
    programmes: dict[str, list[str]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    selected = qc[qc_column].to_numpy(dtype=bool)
    expression = log_normalize(counts[selected, :])
    condition = qc.loc[selected, "condition"].to_numpy()
    mapping = symbol_index(genes)
    endpoint_rows: list[dict[str, object]] = []
    coverage_rows: list[dict[str, object]] = []
    rarefaction_rows: list[pd.DataFrame] = []
    gene_effect_rows: list[pd.DataFrame] = []

    for endpoint_id, declared in programmes.items():
        score, retained, missing, gene_reference = programme_score(
            expression, mapping, declared
        )
        summary = effect_summary(score, condition)
        endpoint_rows.append(
            {
                "qc_definition": qc_column,
                "endpoint_type": "programme",
                "endpoint_id": endpoint_id,
                "declared_gene_count": len(declared),
                "measurable_variable_gene_count": len(retained),
                **summary,
            }
        )
        coverage_rows.append(
            {
                "qc_definition": qc_column,
                "endpoint_id": endpoint_id,
                "declared_gene_count": len(declared),
                "measurable_variable_gene_count": len(retained),
                "coverage_fraction": len(retained) / len(declared),
                "retained_genes": ";".join(retained),
                "missing_or_zero_variance_genes": ";".join(missing),
            }
        )
        rarefaction_rows.append(rarefy(score, condition, qc_column, endpoint_id))

        gene_reference["qc_definition"] = qc_column
        gene_reference["endpoint_id"] = endpoint_id
        gene_reference["mean_treated"] = [
            dense_gene_values(expression, mapping, gene)[condition == "rPlexinB2"].mean()
            for gene in retained
        ]
        gene_reference["mean_vehicle"] = [
            dense_gene_values(expression, mapping, gene)[condition == "vehicle"].mean()
            for gene in retained
        ]
        gene_reference["mean_difference"] = (
            gene_reference["mean_treated"] - gene_reference["mean_vehicle"]
        )
        gene_effect_rows.append(gene_reference)

    for gene in NAMED_GENES:
        if gene not in mapping:
            coverage_rows.append(
                {
                    "qc_definition": qc_column,
                    "endpoint_id": gene,
                    "declared_gene_count": 1,
                    "measurable_variable_gene_count": 0,
                    "coverage_fraction": 0.0,
                    "retained_genes": "",
                    "missing_or_zero_variance_genes": gene,
                }
            )
            continue
        values = dense_gene_values(expression, mapping, gene)
        endpoint_rows.append(
            {
                "qc_definition": qc_column,
                "endpoint_type": "named_gene",
                "endpoint_id": gene,
                "declared_gene_count": 1,
                "measurable_variable_gene_count": int(values.std(ddof=0) > 0),
                **effect_summary(values, condition),
            }
        )
        rarefaction_rows.append(rarefy(values, condition, qc_column, gene))

    return (
        pd.DataFrame(endpoint_rows),
        pd.concat(rarefaction_rows, ignore_index=True),
        pd.DataFrame(coverage_rows),
        pd.concat(gene_effect_rows, ignore_index=True),
    )


def write_summary(effects: pd.DataFrame, rare_summary: pd.DataFrame, qc_summary: pd.DataFrame) -> None:
    primary = effects.loc[effects["qc_definition"] == "primary_qc"].set_index(
        "endpoint_id"
    )

    def line(endpoint: str) -> str:
        row = primary.loc[endpoint]
        rare = rare_summary.loc[
            (rare_summary["qc_definition"] == "primary_qc")
            & (rare_summary["endpoint_id"] == endpoint)
        ].iloc[0]
        return (
            f"- `{endpoint}`: mean difference {row['mean_difference']:+.3f}; "
            f"AUC {row['auc_treated_greater_vehicle']:.3f}; equal-cell median "
            f"{rare['median_mean_difference']:+.3f} "
            f"({rare['p025_mean_difference']:+.3f} to {rare['p975_mean_difference']:+.3f})."
        )

    qc_header = (
        "| Condition | Input cells | Primary QC | Strict QC | Median genes | "
        "Median UMIs | Median mitochondrial % |"
    )
    qc_rule = "|---|---:|---:|---:|---:|---:|---:|"
    qc_rows = [
        (
            f"| {row.condition} | {int(row.input_cells)} | "
            f"{int(row.primary_qc_cells)} | {int(row.strict_qc_cells)} | "
            f"{row.median_features:.0f} | {row.median_library_size:.0f} | "
            f"{row.median_pct_mito:.2f} |"
        )
        for row in qc_summary.itertuples(index=False)
    ]

    text = [
        "# GSE267981 Plexin B2 functional-reference projection",
        "",
        "This analysis contains one treated and one vehicle library. Cell-level summaries and equal-cell resampling ranges are descriptive measurement/composition checks, not replicated intervention inference.",
        "",
        "## QC inventory",
        "",
        qc_header,
        qc_rule,
        *qc_rows,
        "",
        "## Frozen phenotype projections",
        "",
        line("REC_PROGRAM_TOP50"),
        line("REACTOME_TIGHT_JUNCTION_INTERACTIONS"),
        line("CORE_HRC"),
        line("PARTIAL_EMT"),
        line("ACTIN_TURNOVER"),
        line("HALLMARK_E2F_TARGETS"),
        line("HALLMARK_G2M_CHECKPOINT"),
        "",
        "## Post-inspection focused junction audit",
        "",
        line("JUNB_TIGHT_JUNCTION_TARGETS_5_POSTHOC"),
        line("EPCAM"),
        line("CLDN3"),
        line("CLDN4"),
        line("CLDN7"),
        "",
        "## Interpretation boundary",
        "",
        "Plexin B2 does not reproduce the complete REC/HRC–broad-junction–low-cycle phenotype under the frozen equal-gene scoring and threshold-defined population. It evokes a narrower epithelial/junction component: EPCAM, CLDN3, CLDN4 and the existing five-gene focused junction set increase, whereas broad tight junction, HRC and REC change little. Because there is one library per condition, this is a descriptive functional reference rather than replicated intervention inference. It cannot show that Plexin B2 is enriched in, or causes, replacement HGP.",
        "",
    ]
    (OUTPUT / "analysis_summary.md").write_text("\n".join(text), encoding="utf-8")


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for path in FILES.values():
        if not path.exists():
            raise FileNotFoundError(path)
    programmes = load_programmes()
    counts, genes, cells = combine_libraries()
    qc = qc_table(counts, genes, cells)

    all_effects = []
    all_rarefaction = []
    all_coverage = []
    all_gene_effects = []
    for qc_column in ("primary_qc", "strict_qc"):
        effects, rarefaction, coverage, gene_effects = analyse_qc(
            counts, genes, qc, qc_column, programmes
        )
        all_effects.append(effects)
        all_rarefaction.append(rarefaction)
        all_coverage.append(coverage)
        all_gene_effects.append(gene_effects)

    effects = pd.concat(all_effects, ignore_index=True)
    rarefaction = pd.concat(all_rarefaction, ignore_index=True)
    coverage = pd.concat(all_coverage, ignore_index=True)
    gene_effects = pd.concat(all_gene_effects, ignore_index=True)
    rare_summary = (
        rarefaction.groupby(["qc_definition", "endpoint_id"], as_index=False)
        .agg(
            median_mean_difference=("mean_difference", "median"),
            p025_mean_difference=("mean_difference", lambda values: np.quantile(values, 0.025)),
            p975_mean_difference=("mean_difference", lambda values: np.quantile(values, 0.975)),
            median_median_difference=("median_difference", "median"),
            p025_median_difference=("median_difference", lambda values: np.quantile(values, 0.025)),
            p975_median_difference=("median_difference", lambda values: np.quantile(values, 0.975)),
            cells_per_condition=("cells_per_condition", "first"),
        )
    )
    qc_summary = (
        qc.groupby("condition", as_index=False)
        .agg(
            input_cells=("barcode", "size"),
            primary_qc_cells=("primary_qc", "sum"),
            strict_qc_cells=("strict_qc", "sum"),
            median_features=("n_features", "median"),
            median_library_size=("library_size", "median"),
            median_pct_mito=("pct_mito", "median"),
        )
        .sort_values("condition")
    )

    qc.to_csv(OUTPUT / "gse267981_cell_qc.tsv.gz", sep="\t", index=False)
    qc_summary.to_csv(OUTPUT / "gse267981_qc_summary.tsv", sep="\t", index=False)
    effects.to_csv(OUTPUT / "gse267981_endpoint_effects.tsv", sep="\t", index=False)
    rarefaction.to_csv(
        OUTPUT / "gse267981_equal_cell_rarefaction_iterations.tsv.gz",
        sep="\t",
        index=False,
    )
    rare_summary.to_csv(
        OUTPUT / "gse267981_equal_cell_rarefaction_summary.tsv", sep="\t", index=False
    )
    coverage.to_csv(OUTPUT / "gse267981_programme_coverage.tsv", sep="\t", index=False)
    gene_effects.to_csv(
        OUTPUT / "gse267981_programme_gene_effects.tsv.gz", sep="\t", index=False
    )
    write_summary(effects, rare_summary, qc_summary)

    provenance = {
        "analysis": "GSE267981 Plexin B2 functional-reference projection",
        "specification": str(
            ROOT / "metadata" / "plexinb2_functional_projection_spec_2026-09-01.md"
        ),
        "base_seed": BASE_SEED,
        "rarefaction_iterations": N_RAREFACTION,
        "input_files": {
            condition: {"path": str(path), "sha256": sha256(path)}
            for condition, path in FILES.items()
        },
        "gene_symbol_transfer": "case-insensitive exact human-mouse symbol identity",
        "inferential_boundary": "one library per condition; no replicated treatment inference",
        "outputs": sorted(path.name for path in OUTPUT.iterdir() if path.is_file()),
    }
    (OUTPUT / "provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
