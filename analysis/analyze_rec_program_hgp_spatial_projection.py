#!/usr/bin/env python3
"""Project the fixed Ogden REC program into E-MTAB-12043 HGP spatial data."""

from __future__ import annotations

import math
import platform
import csv
import itertools
from collections import deque
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "analysis_work" / "e_mtab_12043_spatial_full"
MANIFEST = ROOT / "data_sources" / "E_MTAB_HGP_CRLM" / "E-MTAB-12043_L1_spatial_pilot_download_manifest.tsv"
PROGRAM_FILE = ROOT / "analysis_results" / "ogden_anchor_program" / "selected_anchor_program_top50.tsv"
SPEC_FILE = ROOT / "metadata" / "rec_program_hgp_spatial_projection_spec_2026-08-23.md"
OUTPUT = ROOT / "analysis_results" / "rec_program_hgp_spatial_projection"
QC_THRESHOLDS = [100, 200, 500]
PRIMARY_THRESHOLD = 200
PRIMARY_BOUNDARY_HOPS = 5
BOUNDARY_HOP_SENSITIVITY = [3, 5, 7]
MIN_PROGRAM_GENES = 40
MIN_SPOTS = 20
BOOTSTRAP_REPLICATES = 10_000
SEED = 42
REGIONS = ["near_interface_0_500um", "deep_tumour_gt500um", "tumour_side"]
RESTRICTIONS = ["all_tumour_side", "epithelial_high"]
MARKERS = {
    # Keep spatial-region construction disjoint from the frozen REC programme.
    # KRT20 was removed after an overlap check identified it as programme rank 9.
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CEACAM5"],
    "hepatocyte": ["ALB", "APOA1", "APOA2", "ASGR1", "CPS1", "TTR"],
    "endothelial": ["VWF", "PECAM1", "EMCN", "KDR", "ENG", "RAMP2", "PLVAP", "ESAM"],
    "perivascular": ["RGS5", "PDGFRB", "CSPG4", "MCAM", "NOTCH3", "DES", "ACTA2"],
    "fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"],
}


def load_metadata() -> pd.DataFrame:
    with MANIFEST.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    unique: dict[str, dict[str, str]] = {}
    for row in rows:
        sample = row["source_name"].replace(" ", "")
        unique[sample] = {"sample": sample, "patient": row["patient"], "hgp": row["hgp"]}
    result = pd.DataFrame(unique.values()).sort_values("sample").reset_index(drop=True)
    if len(result) != 6 or result["hgp"].value_counts().to_dict() != {"rHGP": 3, "dHGP": 3}:
        raise ValueError("Expected three rHGP and three dHGP spatial patients")
    return result


def decode_strings(values: np.ndarray) -> np.ndarray:
    return np.asarray([
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in values
    ])


def read_10x_h5(path: Path) -> tuple[sparse.csr_matrix, np.ndarray, np.ndarray]:
    """Read a 10x feature-by-barcode CSC matrix and return spots by genes."""
    with h5py.File(path, "r") as handle:
        group = handle["matrix"]
        shape = tuple(int(value) for value in group["shape"][:])
        feature_by_barcode = sparse.csc_matrix(
            (group["data"][:], group["indices"][:], group["indptr"][:]),
            shape=shape,
        )
        names = decode_strings(group["features"]["name"][:]).astype(str)
        barcodes = decode_strings(group["barcodes"][:]).astype(str)
    return feature_by_barcode.transpose().tocsr(), names, barcodes


def read_positions(outs: Path) -> pd.DataFrame:
    current = outs / "spatial" / "tissue_positions.csv"
    legacy = outs / "spatial" / "tissue_positions_list.csv"
    if current.is_file():
        positions = pd.read_csv(current)
        if "barcode" not in positions.columns:
            positions = positions.rename(columns={positions.columns[0]: "barcode"})
    elif legacy.is_file():
        positions = pd.read_csv(
            legacy,
            header=None,
            names=[
                "barcode", "in_tissue", "array_row", "array_col",
                "pxl_row_in_fullres", "pxl_col_in_fullres",
            ],
        )
    else:
        raise FileNotFoundError(f"No tissue positions file under {outs / 'spatial'}")
    positions["barcode"] = positions["barcode"].astype(str)
    return positions.set_index("barcode")


def gene_vector(matrix: sparse.csr_matrix, names: np.ndarray, gene: str) -> np.ndarray:
    indices = np.flatnonzero(names == gene)
    if len(indices) == 0:
        return np.zeros(matrix.shape[0], dtype=float)
    return np.asarray(matrix[:, indices].sum(axis=1)).ravel().astype(float)


def zscore(values: np.ndarray) -> np.ndarray:
    standard_deviation = float(np.std(values, ddof=0))
    if not np.isfinite(standard_deviation) or standard_deviation == 0:
        return np.zeros_like(values, dtype=float)
    return (values - float(np.mean(values))) / standard_deviation


def marker_score(log_matrix: sparse.csr_matrix, names: np.ndarray, genes: list[str]) -> tuple[np.ndarray, int]:
    present = [gene for gene in genes if np.any(names == gene)]
    if not present:
        return np.zeros(log_matrix.shape[0], dtype=float), 0
    score = np.mean([zscore(gene_vector(log_matrix, names, gene)) for gene in present], axis=0)
    return score, len(present)


def within_two_hops(rows: np.ndarray, cols: np.ndarray, origins: np.ndarray) -> np.ndarray:
    coordinate_to_index = {(int(row), int(col)): index for index, (row, col) in enumerate(zip(rows, cols))}
    offsets = ((0, -2), (0, 2), (-1, -1), (-1, 1), (1, -1), (1, 1))
    reached = set(np.flatnonzero(origins).tolist())
    frontier = set(reached)
    for _ in range(2):
        next_frontier: set[int] = set()
        for index in frontier:
            row, col = int(rows[index]), int(cols[index])
            for row_offset, col_offset in offsets:
                neighbor = coordinate_to_index.get((row + row_offset, col + col_offset))
                if neighbor is not None and neighbor not in reached:
                    next_frontier.add(neighbor)
        reached.update(next_frontier)
        frontier = next_frontier
    result = np.zeros(len(rows), dtype=bool)
    result[list(reached)] = True
    return result


def build_masks(scores: dict[str, np.ndarray], rows: np.ndarray, cols: np.ndarray) -> dict[str, np.ndarray]:
    delta = scores["epithelial"] - scores["hepatocyte"]
    lower, upper = np.quantile(delta, [1 / 3, 2 / 3])
    liver_side = delta <= lower
    tumour_side = delta >= upper
    return {
        "tumour_side": tumour_side,
        "liver_side": liver_side,
        "interface": tumour_side & within_two_hops(rows, cols, liver_side),
    }


def hedges_g(r_values: np.ndarray, d_values: np.ndarray) -> float:
    n_r, n_d = len(r_values), len(d_values)
    pooled_variance = (
        (n_r - 1) * np.var(r_values, ddof=1) + (n_d - 1) * np.var(d_values, ddof=1)
    ) / (n_r + n_d - 2)
    if pooled_variance <= 0:
        return math.nan
    correction = 1 - 3 / (4 * (n_r + n_d) - 9)
    return float(correction * (np.mean(r_values) - np.mean(d_values)) / math.sqrt(pooled_variance))


def exact_permutation(values: np.ndarray, observed: float) -> tuple[float, float]:
    differences: list[float] = []
    all_indices = set(range(6))
    for rhgp_tuple in itertools.combinations(range(6), 3):
        rhgp_indices = set(rhgp_tuple)
        dhgp_indices = sorted(all_indices - rhgp_indices)
        differences.append(float(np.mean(values[sorted(rhgp_indices)]) - np.mean(values[dhgp_indices])))
    tolerance = 1e-12
    one_sided = sum(value >= observed - tolerance for value in differences) / len(differences)
    two_sided = sum(abs(value) >= abs(observed) - tolerance for value in differences) / len(differences)
    return one_sided, two_sided


def exact_sign_flip(values: np.ndarray) -> float:
    """Return the exhaustive two-sided sign-flip P value for paired differences."""
    observed = float(np.mean(values))
    null_means = [
        float(np.mean(values * np.asarray(signs, dtype=float)))
        for signs in itertools.product((-1.0, 1.0), repeat=len(values))
    ]
    tolerance = 1e-12
    return sum(abs(value) >= abs(observed) - tolerance for value in null_means) / len(null_means)


def distance_from_origins(rows: np.ndarray, cols: np.ndarray, origins: np.ndarray) -> np.ndarray:
    """Return shortest hex-lattice graph distance from any origin spot."""
    coordinate_to_index = {(int(row), int(col)): index for index, (row, col) in enumerate(zip(rows, cols))}
    offsets = ((0, -2), (0, 2), (-1, -1), (-1, 1), (1, -1), (1, 1))
    distance = np.full(len(rows), np.inf)
    queue: deque[int] = deque()
    for index in np.flatnonzero(origins):
        distance[index] = 0.0
        queue.append(int(index))
    while queue:
        index = queue.popleft()
        row, col = int(rows[index]), int(cols[index])
        for row_offset, col_offset in offsets:
            neighbor = coordinate_to_index.get((row + row_offset, col + col_offset))
            if neighbor is not None and not np.isfinite(distance[neighbor]):
                distance[neighbor] = distance[index] + 1.0
                queue.append(neighbor)
    return distance


def load_program() -> list[str]:
    table = pd.read_csv(PROGRAM_FILE, sep="\t")
    required = {"gene", "program_rank", "direction"}
    if not required.issubset(table.columns):
        raise ValueError(f"Program file lacks columns: {sorted(required - set(table.columns))}")
    table = table.sort_values("program_rank")
    if len(table) != 50 or table["gene"].duplicated().any() or not (table["direction"] == "anchor_up").all():
        raise ValueError("Expected 50 unique REC-up program genes")
    return table["gene"].astype(str).str.upper().tolist()


def region_pseudobulk(
    matrix: sparse.csr_matrix,
    names: np.ndarray,
    library: np.ndarray,
    masks: dict[str, np.ndarray],
    sample: str,
    threshold: int,
    restriction: str,
    program_genes: list[str],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    expression_rows: list[dict[str, object]] = []
    mask_rows: list[dict[str, object]] = []
    for region in REGIONS:
        mask = masks[region]
        n_spots = int(mask.sum())
        total_library = float(library[mask].sum()) if n_spots else 0.0
        mask_rows.append({
            "sample": sample,
            "qc_min_detected_genes": threshold,
            "restriction": restriction,
            "region": region,
            "n_spots": n_spots,
            "evaluable": n_spots >= MIN_SPOTS,
            "library_size": total_library,
        })
        for gene in program_genes:
            count = float(gene_vector(matrix[mask, :], names, gene).sum()) if n_spots else 0.0
            log2_cpm = (
                math.log2((count + 0.5) / (total_library + 1.0) * 1e6)
                if n_spots else math.nan
            )
            expression_rows.append({
                "sample": sample,
                "qc_min_detected_genes": threshold,
                "restriction": restriction,
                "region": region,
                "gene": gene,
                "n_spots": n_spots,
                "evaluable": n_spots >= MIN_SPOTS,
                "raw_count": count,
                "library_size": total_library,
                "log2_cpm": log2_cpm,
            })
    return expression_rows, mask_rows


def load_sample(
    sample: str,
    threshold: int,
    requested_program: list[str],
    boundary_hops: int = PRIMARY_BOUNDARY_HOPS,
) -> tuple[dict[str, object], list[dict[str, object]], list[dict[str, object]], set[str]]:
    outs = INPUT / sample / "outs"
    matrix, names, barcodes = read_10x_h5(outs / "filtered_feature_bc_matrix.h5")
    matrix_spots = int(matrix.shape[0])
    names = np.asarray(pd.Index(names).str.upper())
    positions = read_positions(outs).reindex(barcodes)
    if positions[["in_tissue", "array_row", "array_col"]].isna().any().any():
        raise ValueError(f"Barcode/position mismatch for {sample}")
    in_tissue = positions["in_tissue"].to_numpy(dtype=int) == 1
    detected = np.diff(matrix.indptr)
    keep = in_tissue & (detected >= threshold)
    matrix = matrix[keep, :].tocsr()
    kept_positions = positions.iloc[np.flatnonzero(keep)]
    library = np.asarray(matrix.sum(axis=1)).ravel().astype(float)
    factors = np.divide(1e4, library, out=np.zeros_like(library), where=library > 0)
    log_matrix = matrix.multiply(factors[:, None]).tocsr()
    log_matrix.data = np.log1p(log_matrix.data)

    marker_scores: dict[str, np.ndarray] = {}
    marker_coverage: dict[str, int] = {}
    for marker_name, marker_genes in MARKERS.items():
        marker_scores[marker_name], marker_coverage[marker_name] = marker_score(
            log_matrix, names, marker_genes
        )
    rows = kept_positions["array_row"].to_numpy(dtype=int)
    cols = kept_positions["array_col"].to_numpy(dtype=int)
    base_masks = build_masks(marker_scores, rows, cols)
    distance = distance_from_origins(rows, cols, base_masks["liver_side"])
    connected_tumour = base_masks["tumour_side"] & np.isfinite(distance)
    near = connected_tumour & (distance >= 1) & (distance <= boundary_hops)
    deep = connected_tumour & (distance >= boundary_hops + 1)

    epithelial_high = np.zeros(len(rows), dtype=bool)
    if np.any(base_masks["tumour_side"]):
        cutoff = float(np.median(marker_scores["epithelial"][base_masks["tumour_side"]]))
        epithelial_high = base_masks["tumour_side"] & (marker_scores["epithelial"] >= cutoff)

    masks_by_restriction = {
        "all_tumour_side": {
            "near_interface_0_500um": near,
            "deep_tumour_gt500um": deep,
            "tumour_side": base_masks["tumour_side"],
        },
        "epithelial_high": {
            "near_interface_0_500um": near & epithelial_high,
            "deep_tumour_gt500um": deep & epithelial_high,
            "tumour_side": base_masks["tumour_side"] & epithelial_high,
        },
    }
    present = set(requested_program).intersection(set(names.tolist()))
    expression_rows: list[dict[str, object]] = []
    mask_rows: list[dict[str, object]] = []
    for restriction, restriction_masks in masks_by_restriction.items():
        expressions, counts = region_pseudobulk(
            matrix,
            names,
            library,
            restriction_masks,
            sample,
            threshold,
            restriction,
            sorted(present, key=requested_program.index),
        )
        expression_rows.extend(expressions)
        mask_rows.extend(counts)
    qc = {
        "sample": sample,
        "qc_min_detected_genes": threshold,
        "boundary_hops": boundary_hops,
        "matrix_spots": matrix_spots,
        "in_tissue_spots": int(in_tissue.sum()),
        "retained_spots": int(keep.sum()),
        "median_detected_genes": float(np.median(detected[keep])) if np.any(keep) else math.nan,
        "median_umi": float(np.median(library)) if len(library) else math.nan,
        "tumour_side_spots": int(base_masks["tumour_side"].sum()),
        "liver_side_spots": int(base_masks["liver_side"].sum()),
        "connected_tumour_spots": int(connected_tumour.sum()),
        "program_genes_present": len(present),
        **{f"{key}_markers_present": value for key, value in marker_coverage.items()},
    }
    return qc, mask_rows, expression_rows, present


def calculate_boundary_hop_sensitivity(
    metadata: pd.DataFrame,
    requested_program: list[str],
    common_program: list[str],
) -> pd.DataFrame:
    """Summarise boundary localisation at 3, 5 and 7 graph hops."""
    rows: list[dict[str, object]] = []
    for boundary_hops in BOUNDARY_HOP_SENSITIVITY:
        expression_rows: list[dict[str, object]] = []
        for sample in metadata["sample"]:
            _, _, sample_expression, _ = load_sample(
                sample,
                PRIMARY_THRESHOLD,
                requested_program,
                boundary_hops=boundary_hops,
            )
            expression_rows.extend(sample_expression)
        expression = pd.DataFrame(expression_rows)
        expression = expression[expression["gene"].isin(common_program)].copy()
        scores, patients, _, _ = calculate_scores_and_effects(expression, metadata, common_program)
        selected_patients = patients[
            (patients["qc_min_detected_genes"] == PRIMARY_THRESHOLD)
            & (patients["restriction"] == "all_tumour_side")
        ].copy()
        selected_scores = scores[
            (scores["qc_min_detected_genes"] == PRIMARY_THRESHOLD)
            & (scores["restriction"] == "all_tumour_side")
        ].copy()
        if len(selected_patients) != 6 or len(selected_scores) != 12:
            raise ValueError(
                f"Boundary-hop sensitivity at {boundary_hops} hops did not retain all six patients"
            )
        spot_counts = selected_scores.pivot(
            index="patient", columns="region", values="n_spots"
        )
        endpoint_columns = {
            "Equal-gene score": "localization_index",
            "Pooled expression": "pooled_localization_index",
        }
        is_rhgp = selected_patients["hgp"].to_numpy() == "rHGP"
        for endpoint_label, endpoint_column in endpoint_columns.items():
            values = selected_patients[endpoint_column].to_numpy(dtype=float)
            contrast = group_contrast(
                values,
                is_rhgp,
                SEED + 1000 + 10 * boundary_hops + len(rows),
            )
            rows.append({
                "boundary_hops": boundary_hops,
                "near_band_definition": f"1-{boundary_hops} graph hops",
                "deep_band_definition": f">={boundary_hops + 1} graph hops",
                "endpoint": endpoint_label,
                "n_patients": len(values),
                "minimum_near_spots": int(spot_counts["near_interface_0_500um"].min()),
                "minimum_deep_spots": int(spot_counts["deep_tumour_gt500um"].min()),
                "mean_near_minus_deep": float(np.mean(values)),
                "positive_patients": int((values > 0).sum()),
                "exact_sign_flip_p_two_sided": exact_sign_flip(values),
                "rhgp_mean_near_minus_deep": contrast["rhgp_mean"],
                "dhgp_mean_near_minus_deep": contrast["dhgp_mean"],
                "rhgp_minus_dhgp_localization": contrast["rhgp_minus_dhgp"],
                "bootstrap_ci_lower": contrast["bootstrap_ci_lower"],
                "bootstrap_ci_upper": contrast["bootstrap_ci_upper"],
                "exact_hgp_permutation_p_two_sided": contrast["exact_p_two_sided"],
            })
    return pd.DataFrame(rows)


def bootstrap_mean_difference(
    values: np.ndarray,
    is_rhgp: np.ndarray,
    seed: int,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    r_values = values[is_rhgp]
    d_values = values[~is_rhgp]
    r_indices = rng.integers(0, len(r_values), size=(BOOTSTRAP_REPLICATES, len(r_values)))
    d_indices = rng.integers(0, len(d_values), size=(BOOTSTRAP_REPLICATES, len(d_values)))
    differences = r_values[r_indices].mean(axis=1) - d_values[d_indices].mean(axis=1)
    lower, upper = np.quantile(differences, [0.025, 0.975])
    return float(lower), float(upper)


def group_contrast(
    values: np.ndarray,
    is_rhgp: np.ndarray,
    seed: int,
) -> dict[str, float]:
    observed = float(values[is_rhgp].mean() - values[~is_rhgp].mean())
    one_sided, two_sided = exact_permutation(values, observed)
    ci_lower, ci_upper = bootstrap_mean_difference(values, is_rhgp, seed)
    return {
        "rhgp_mean": float(values[is_rhgp].mean()),
        "dhgp_mean": float(values[~is_rhgp].mean()),
        "rhgp_minus_dhgp": observed,
        "bootstrap_ci_lower": ci_lower,
        "bootstrap_ci_upper": ci_upper,
        "hedges_g": hedges_g(values[is_rhgp], values[~is_rhgp]),
        "exact_p_one_sided": one_sided,
        "exact_p_two_sided": two_sided,
    }


def holm_adjust_two(p_values: list[float]) -> list[float]:
    order = np.argsort(p_values)
    adjusted = np.empty(len(p_values), dtype=float)
    running = 0.0
    for rank, index in enumerate(order):
        candidate = (len(p_values) - rank) * p_values[index]
        running = max(running, candidate)
        adjusted[index] = min(running, 1.0)
    return adjusted.tolist()


def calculate_scores_and_effects(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    program_genes: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    annotated = expression.merge(metadata, on="sample", validate="many_to_one")
    score_rows: list[dict[str, object]] = []
    patient_rows: list[dict[str, object]] = []
    contrast_rows: list[dict[str, object]] = []
    loo_rows: list[dict[str, object]] = []
    contrast_counter = 0

    for (threshold, restriction), frame in annotated.groupby(
        ["qc_min_detected_genes", "restriction"], sort=True
    ):
        spatial_frame = frame[frame["region"].isin(["near_interface_0_500um", "deep_tumour_gt500um"])]
        wide = spatial_frame.pivot(
            index=["sample", "patient", "hgp", "region", "evaluable", "n_spots", "library_size"],
            columns="gene",
            values="log2_cpm",
        ).reset_index()
        if len(wide) != 12 or not wide["evaluable"].astype(bool).all():
            continue
        available = [gene for gene in program_genes if gene in wide.columns and wide[gene].notna().all()]
        gene_matrix = wide[available].to_numpy(dtype=float)
        means = gene_matrix.mean(axis=0)
        standard_deviations = gene_matrix.std(axis=0, ddof=0)
        variable = standard_deviations > 0
        if int(variable.sum()) < MIN_PROGRAM_GENES:
            continue
        standardized = (gene_matrix[:, variable] - means[variable]) / standard_deviations[variable]
        wide["rec_program_score"] = standardized.mean(axis=1)

        pooled = spatial_frame.groupby(
            ["sample", "patient", "hgp", "region"], as_index=False
        ).agg(program_raw_count=("raw_count", "sum"), library_size=("library_size", "first"))
        pooled["pooled_program_log2_cpm"] = np.log2(
            (pooled["program_raw_count"] + 0.5) / (pooled["library_size"] + 1.0) * 1e6
        )
        wide = wide.merge(
            pooled[["sample", "patient", "hgp", "region", "pooled_program_log2_cpm"]],
            on=["sample", "patient", "hgp", "region"],
            validate="one_to_one",
        )
        for _, row in wide.iterrows():
            score_rows.append({
                "sample": row["sample"],
                "patient": row["patient"],
                "hgp": row["hgp"],
                "qc_min_detected_genes": threshold,
                "restriction": restriction,
                "region": row["region"],
                "n_spots": row["n_spots"],
                "n_program_genes_scored": int(variable.sum()),
                "rec_program_score": row["rec_program_score"],
                "pooled_program_log2_cpm": row["pooled_program_log2_cpm"],
            })

        program_wide = wide.pivot(
            index=["sample", "patient", "hgp"], columns="region", values="rec_program_score"
        ).reset_index()
        pooled_wide = wide.pivot(
            index=["sample", "patient", "hgp"], columns="region", values="pooled_program_log2_cpm"
        ).reset_index()
        patients = program_wide.merge(pooled_wide, on=["sample", "patient", "hgp"], suffixes=("_score", "_pooled"))
        patients["localization_index"] = (
            patients["near_interface_0_500um_score"] - patients["deep_tumour_gt500um_score"]
        )
        patients["pooled_localization_index"] = (
            patients["near_interface_0_500um_pooled"] - patients["deep_tumour_gt500um_pooled"]
        )
        for _, row in patients.iterrows():
            patient_rows.append({
                "sample": row["sample"],
                "patient": row["patient"],
                "hgp": row["hgp"],
                "qc_min_detected_genes": threshold,
                "restriction": restriction,
                "near_interface_program_score": row["near_interface_0_500um_score"],
                "deep_tumour_program_score": row["deep_tumour_gt500um_score"],
                "localization_index": row["localization_index"],
                "near_interface_pooled_program_log2_cpm": row["near_interface_0_500um_pooled"],
                "deep_tumour_pooled_program_log2_cpm": row["deep_tumour_gt500um_pooled"],
                "pooled_localization_index": row["pooled_localization_index"],
            })

        is_rhgp = patients["hgp"].to_numpy() == "rHGP"
        endpoints = {
            "near_interface_program_score": patients["near_interface_0_500um_score"].to_numpy(dtype=float),
            "localization_index": patients["localization_index"].to_numpy(dtype=float),
            "near_interface_pooled_program_log2_cpm": patients["near_interface_0_500um_pooled"].to_numpy(dtype=float),
            "pooled_localization_index": patients["pooled_localization_index"].to_numpy(dtype=float),
        }
        current_indices: list[int] = []
        for endpoint, values in endpoints.items():
            contrast = group_contrast(values, is_rhgp, SEED + contrast_counter)
            contrast_counter += 1
            contrast_rows.append({
                "qc_min_detected_genes": threshold,
                "restriction": restriction,
                "endpoint": endpoint,
                "n_rhgp": int(is_rhgp.sum()),
                "n_dhgp": int((~is_rhgp).sum()),
                **contrast,
            })
            current_indices.append(len(contrast_rows) - 1)
            for omitted_index, omitted in patients.iterrows():
                keep = np.arange(len(patients)) != omitted_index
                effect = float(values[keep & is_rhgp].mean() - values[keep & ~is_rhgp].mean())
                loo_rows.append({
                    "qc_min_detected_genes": threshold,
                    "restriction": restriction,
                    "endpoint": endpoint,
                    "omitted_patient": omitted["patient"],
                    "omitted_hgp": omitted["hgp"],
                    "rhgp_minus_dhgp": effect,
                    "positive": effect > 0,
                })
        primary_indices = current_indices[:2]
        adjusted = holm_adjust_two([contrast_rows[index]["exact_p_two_sided"] for index in primary_indices])
        for index, p_adjusted in zip(primary_indices, adjusted):
            contrast_rows[index]["holm_p_two_primary_endpoints"] = p_adjusted

        tumour_frame = frame[frame["region"] == "tumour_side"]
        tumour_wide = tumour_frame.pivot(
            index=["sample", "patient", "hgp", "evaluable", "n_spots"],
            columns="gene",
            values="log2_cpm",
        ).reset_index()
        tumour_available = [gene for gene in program_genes if gene in tumour_wide.columns and tumour_wide[gene].notna().all()]
        if len(tumour_wide) == 6 and tumour_wide["evaluable"].astype(bool).all():
            matrix = tumour_wide[tumour_available].to_numpy(dtype=float)
            means = matrix.mean(axis=0)
            standard_deviations = matrix.std(axis=0, ddof=0)
            variable = standard_deviations > 0
            tumour_wide["rec_program_score"] = (
                (matrix[:, variable] - means[variable]) / standard_deviations[variable]
            ).mean(axis=1)
            values = tumour_wide["rec_program_score"].to_numpy(dtype=float)
            tumour_is_rhgp = tumour_wide["hgp"].to_numpy() == "rHGP"
            contrast = group_contrast(values, tumour_is_rhgp, SEED + contrast_counter)
            contrast_counter += 1
            contrast_rows.append({
                "qc_min_detected_genes": threshold,
                "restriction": restriction,
                "endpoint": "whole_tumour_side_program_score",
                "n_rhgp": int(tumour_is_rhgp.sum()),
                "n_dhgp": int((~tumour_is_rhgp).sum()),
                **contrast,
            })

    contrasts = pd.DataFrame(contrast_rows)
    if "holm_p_two_primary_endpoints" not in contrasts.columns:
        contrasts["holm_p_two_primary_endpoints"] = math.nan
    return pd.DataFrame(score_rows), pd.DataFrame(patient_rows), contrasts, pd.DataFrame(loo_rows)


def effect_value(contrasts: pd.DataFrame, threshold: int, restriction: str, endpoint: str) -> float:
    selected = contrasts[
        (contrasts["qc_min_detected_genes"] == threshold)
        & (contrasts["restriction"] == restriction)
        & (contrasts["endpoint"] == endpoint)
    ]
    return float(selected.iloc[0]["rhgp_minus_dhgp"]) if len(selected) == 1 else math.nan


def write_summary(
    coverage: pd.DataFrame,
    patients: pd.DataFrame,
    contrasts: pd.DataFrame,
    loo: pd.DataFrame,
    decision: pd.DataFrame,
    boundary_sensitivity: pd.DataFrame,
) -> None:
    primary_patients = patients[
        (patients["qc_min_detected_genes"] == PRIMARY_THRESHOLD)
        & (patients["restriction"] == "all_tumour_side")
    ].copy()
    primary_contrasts = contrasts[
        (contrasts["qc_min_detected_genes"] == PRIMARY_THRESHOLD)
        & (contrasts["restriction"] == "all_tumour_side")
        & (contrasts["endpoint"].isin(["near_interface_program_score", "localization_index"]))
    ]
    endpoint_rows = {row["endpoint"]: row for _, row in primary_contrasts.iterrows()}
    near = endpoint_rows.get("near_interface_program_score")
    localization = endpoint_rows.get("localization_index")
    rhgp_localization = primary_patients.loc[primary_patients["hgp"] == "rHGP", "localization_index"]
    dhgp_localization = primary_patients.loc[primary_patients["hgp"] == "dHGP", "localization_index"]
    lines = [
        "# Fixed REC-program projection into E-MTAB-12043",
        "",
        f"Decision: **{decision.iloc[0]['decision']}**.",
        "",
        "## Coverage",
        "",
        f"All six spatial matrices contained {int(coverage.iloc[0]['n_common_present'])} of the 50 frozen REC-program genes. "
        f"The projection used the same genes in every patient.",
        "",
        "## Primary patient-level results",
        "",
    ]
    if near is not None:
        lines.append(
            f"At the tumour-side 0–500 µm band, the mean standardized program score was "
            f"{near['rhgp_mean']:.3f} in rHGP and {near['dhgp_mean']:.3f} in dHGP "
            f"(difference {near['rhgp_minus_dhgp']:+.3f}, bootstrap 95% CI "
            f"{near['bootstrap_ci_lower']:+.3f} to {near['bootstrap_ci_upper']:+.3f}, "
            f"Hedges' g {near['hedges_g']:+.3f}, exact two-sided P={near['exact_p_two_sided']:.3f})."
        )
        lines.append("")
    if localization is not None:
        lines.append(
            f"The rHGP-minus-dHGP difference in the within-patient near-minus-deep localization index was "
            f"{localization['rhgp_minus_dhgp']:+.3f} (bootstrap 95% CI "
            f"{localization['bootstrap_ci_lower']:+.3f} to {localization['bootstrap_ci_upper']:+.3f}, "
            f"Hedges' g {localization['hedges_g']:+.3f}, exact two-sided P={localization['exact_p_two_sided']:.3f})."
        )
        lines.append("")
    lines.extend([
        f"Near-minus-deep scores were positive in {int((rhgp_localization > 0).sum())}/3 rHGP patients and "
        f"{int((dhgp_localization > 0).sum())}/3 dHGP patients.",
        "",
        "## Interpretation",
        "",
        str(decision.iloc[0]["interpretation"]),
        "",
        "The spatial bands are computationally inferred from epithelial and hepatocyte expression and graph distance. "
        "They support regional localization, not direct cell contact or causal regulation.",
        "",
        "Across 3-, 5- and 7-hop boundary definitions, the equal-gene and pooled-expression endpoints are reported in "
        "`boundary_hop_sensitivity.tsv`; these analyses are sensitivity checks rather than additional primary tests.",
        "",
        "## Reproducibility",
        "",
        f"Specification: `{SPEC_FILE.relative_to(ROOT)}`. Random seed: {SEED}; patient-level bootstrap replicates: {BOOTSTRAP_REPLICATES}.",
    ])
    (OUTPUT / "analysis_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    requested_program = load_program()
    region_markers = set(MARKERS["epithelial"] + MARKERS["hepatocyte"])
    marker_program_overlap = sorted(region_markers.intersection(requested_program))
    if marker_program_overlap:
        raise ValueError(
            "Spatial-region markers overlap the frozen REC programme: "
            + ", ".join(marker_program_overlap)
        )
    metadata = load_metadata()
    qc_rows: list[dict[str, object]] = []
    mask_rows: list[dict[str, object]] = []
    expression_rows: list[dict[str, object]] = []
    present_by_sample: dict[str, set[str]] = {}
    for threshold in QC_THRESHOLDS:
        for sample in metadata["sample"]:
            qc, masks, expression, present = load_sample(sample, threshold, requested_program)
            qc_rows.append(qc)
            mask_rows.extend(masks)
            expression_rows.extend(expression)
            present_by_sample.setdefault(sample, present)
            print(f"loaded threshold={threshold} sample={sample}", flush=True)

    common_present = set(requested_program)
    for present in present_by_sample.values():
        common_present &= present
    common_program = [gene for gene in requested_program if gene in common_present]
    if len(common_program) < MIN_PROGRAM_GENES:
        raise ValueError(f"Only {len(common_program)} of 50 program genes are common across samples")

    qc = pd.DataFrame(qc_rows).merge(metadata, on="sample", validate="many_to_one")
    masks = pd.DataFrame(mask_rows).merge(metadata, on="sample", validate="many_to_one")
    expression = pd.DataFrame(expression_rows)
    expression = expression[expression["gene"].isin(common_program)].copy()
    scores, patients, contrasts, loo = calculate_scores_and_effects(expression, metadata, common_program)
    boundary_sensitivity = calculate_boundary_hop_sensitivity(
        metadata, requested_program, common_program
    )

    primary_loo = loo[
        (loo["qc_min_detected_genes"] == PRIMARY_THRESHOLD)
        & (loo["restriction"] == "all_tumour_side")
    ]
    near_loo = primary_loo[primary_loo["endpoint"] == "near_interface_program_score"]
    localization_loo = primary_loo[primary_loo["endpoint"] == "localization_index"]
    near_by_threshold = [
        effect_value(contrasts, threshold, "all_tumour_side", "near_interface_program_score")
        for threshold in QC_THRESHOLDS
    ]
    localization_by_threshold = [
        effect_value(contrasts, threshold, "all_tumour_side", "localization_index")
        for threshold in QC_THRESHOLDS
    ]
    coverage_pass = len(common_program) >= MIN_PROGRAM_GENES
    near_stable = all(np.isfinite(value) and value > 0 for value in near_by_threshold)
    localization_stable = all(np.isfinite(value) and value > 0 for value in localization_by_threshold)
    near_loo_pass = len(near_loo) == 6 and int(near_loo["positive"].sum()) >= 5
    localization_loo_pass = len(localization_loo) == 6 and int(localization_loo["positive"].sum()) >= 5
    if coverage_pass and near_stable and near_loo_pass and localization_stable and localization_loo_pass:
        label = "STRONG_SPATIAL_SUPPORT"
        interpretation = (
            "The fixed REC program is higher in rHGP near the tumour–liver boundary and shows a larger "
            "near-versus-deep gradient in rHGP. It can be used as independent regional support for the main story."
        )
    elif coverage_pass and near_stable and near_loo_pass:
        label = "REGIONAL_NOT_PREFERENTIAL_SUPPORT"
        interpretation = (
            "The fixed REC program is reproducibly higher in rHGP in the near-interface tumour band, but the "
            "data do not show stable preferential boundary localization relative to deep tumour."
        )
    else:
        label = "NO_STABLE_HGP_SPATIAL_SUPPORT"
        interpretation = (
            "The fixed REC program does not show a stable rHGP-associated near-interface direction in this "
            "six-patient spatial cohort and should not be presented as independent HGP validation."
        )
    decision = pd.DataFrame([{
        "decision": label,
        "n_program_genes_defined": len(requested_program),
        "n_program_genes_common": len(common_program),
        "near_effect_q100": near_by_threshold[0],
        "near_effect_q200": near_by_threshold[1],
        "near_effect_q500": near_by_threshold[2],
        "near_loo_positive_n": int(near_loo["positive"].sum()) if len(near_loo) else 0,
        "localization_effect_q100": localization_by_threshold[0],
        "localization_effect_q200": localization_by_threshold[1],
        "localization_effect_q500": localization_by_threshold[2],
        "localization_loo_positive_n": int(localization_loo["positive"].sum()) if len(localization_loo) else 0,
        "interpretation": interpretation,
    }])
    coverage = pd.DataFrame([{
        "n_defined": len(requested_program),
        "n_common_present": len(common_program),
        "coverage_fraction": len(common_program) / len(requested_program),
        "matched_genes": ",".join(common_program),
        "missing_genes": ",".join(gene for gene in requested_program if gene not in common_present),
    }])
    detection = (
        expression.groupby("gene", as_index=False)
        .agg(
            total_raw_count=("raw_count", "sum"),
            fraction_region_pseudobulks_detected=("raw_count", lambda values: float((values > 0).mean())),
        )
        .merge(
            pd.DataFrame({"gene": requested_program, "program_rank": range(1, 51)}),
            on="gene",
            how="right",
        )
        .sort_values("program_rank")
    )
    run_info = pd.DataFrame([{
        "analysis_date": "2026-08-24",
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "h5py_version": h5py.__version__,
        "random_seed": SEED,
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "statistical_unit": "patient_section",
        "region_marker_program_overlap_n": 0,
        "region_marker_overlap_correction": "KRT20 removed from epithelial mask",
    }])

    OUTPUT.mkdir(parents=True, exist_ok=True)
    outputs = {
        "program_coverage.tsv": coverage,
        "program_gene_detection.tsv": detection,
        "sample_qc.tsv": qc,
        "region_counts.tsv": masks,
        "patient_region_gene_pseudobulk.tsv.gz": expression.merge(
            metadata, on="sample", validate="many_to_one"
        ),
        "patient_region_scores.tsv": scores,
        "patient_localization_values.tsv": patients,
        "hgp_contrasts.tsv": contrasts,
        "leave_one_patient_out.tsv": loo,
        "boundary_hop_sensitivity.tsv": boundary_sensitivity,
        "decision.tsv": decision,
        "run_info.tsv": run_info,
    }
    for filename, frame in outputs.items():
        compression = "gzip" if filename.endswith(".gz") else None
        frame.to_csv(OUTPUT / filename, sep="\t", index=False, encoding="utf-8", compression=compression)
    write_summary(coverage, patients, contrasts, loo, decision, boundary_sensitivity)
    manifest_lines = [
        "# Analysis output manifest",
        "",
        "Generated by `scripts/analyze_rec_program_hgp_spatial_projection.py`.",
        "",
        "| File | Purpose |",
        "|---|---|",
        "| `analysis_summary.md` | Human-readable result and interpretation |",
        "| `program_coverage.tsv` | Frozen-program coverage across spatial matrices |",
        "| `program_gene_detection.tsv` | Gene-level spatial count coverage |",
        "| `sample_qc.tsv` | Sample-level sequencing and mask QC |",
        "| `region_counts.tsv` | Patient-region spot counts and evaluability |",
        "| `patient_region_gene_pseudobulk.tsv.gz` | Auditable gene-level regional pseudobulks |",
        "| `patient_region_scores.tsv` | Patient-region REC-program scores |",
        "| `patient_localization_values.tsv` | Near, deep and within-patient localization values |",
        "| `hgp_contrasts.tsv` | Patient-level HGP effects, uncertainty and exact tests |",
        "| `leave_one_patient_out.tsv` | Directional influence analysis |",
        "| `boundary_hop_sensitivity.tsv` | Three-, five- and seven-hop boundary localisation sensitivity |",
        "| `decision.tsv` | Prespecified interpretation rule |",
        "| `run_info.tsv` | Software versions and reproducibility settings |",
    ]
    (OUTPUT / "_analysis_outputs.md").write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")
    print(decision.to_string(index=False))


if __name__ == "__main__":
    main()
