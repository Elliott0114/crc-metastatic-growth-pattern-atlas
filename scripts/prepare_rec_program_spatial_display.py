#!/usr/bin/env python3
"""Prepare FFPE spatial maps for reader-facing figures.

This script does not perform a new inferential analysis. It reproduces the
primary (>=200 detected genes per spot) region construction used by
``analyze_rec_program_hgp_spatial_projection.py`` and exports spot coordinates
for visualisation in R. All six specimens are retained for the supplementary
cohort overview; one sample per HGP is additionally selected as the within-group
median for both prespecified near-interface programme summaries in Figure 4.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

import analyze_rec_program_hgp_spatial_projection as spatial


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "analysis_work" / "e_mtab_12043_spatial_full"
SCORES = ROOT / "analysis_results" / "rec_program_hgp_spatial_projection" / "patient_region_scores.tsv"
OUTPUT = ROOT / "analysis_results" / "rec_program_spatial_display"
THRESHOLD = 200


def load_primary_samples() -> pd.DataFrame:
    scores = pd.read_csv(SCORES, sep="\t")
    primary = scores.loc[
        (scores["qc_min_detected_genes"] == THRESHOLD)
        & (scores["restriction"] == "all_tumour_side")
        & (scores["region"] == "near_interface_0_500um"),
        [
            "sample",
            "patient",
            "hgp",
            "rec_program_score",
            "pooled_program_log2_cpm",
        ],
    ].copy()
    if primary.groupby("hgp").size().to_dict() != {"dHGP": 3, "rHGP": 3}:
        raise ValueError("Expected three samples in each HGP group")
    return primary.sort_values(["hgp", "patient", "sample"]).reset_index(drop=True)


def select_representative_samples(primary: pd.DataFrame) -> pd.DataFrame:

    selected_rows: list[pd.Series] = []
    for _, group in primary.groupby("hgp", sort=True):
        group = group.copy()
        group["equal_gene_rank"] = group["rec_program_score"].rank(method="first")
        group["pooled_rank"] = group["pooled_program_log2_cpm"].rank(method="first")
        candidates = group.loc[
            (group["equal_gene_rank"] == 2) & (group["pooled_rank"] == 2)
        ]
        if len(candidates) != 1:
            raise ValueError("The two prespecified endpoints do not share a unique median sample")
        selected_rows.append(candidates.iloc[0])

    selected = pd.DataFrame(selected_rows).sort_values("hgp").reset_index(drop=True)
    selected["selection_rule"] = (
        "within-HGP median for both near-interface equal-gene and pooled-expression summaries"
    )
    return selected


def prepare_sample(row: pd.Series) -> tuple[pd.DataFrame, dict[str, object]]:
    sample = str(row["sample"])
    outs = INPUT / sample / "outs"
    matrix, names, barcodes = spatial.read_10x_h5(outs / "filtered_feature_bc_matrix.h5")
    names = np.asarray(pd.Index(names).str.upper())
    positions = spatial.read_positions(outs).reindex(barcodes)
    if positions[["in_tissue", "array_row", "array_col"]].isna().any().any():
        raise ValueError(f"Barcode/position mismatch for {sample}")

    in_tissue = positions["in_tissue"].to_numpy(dtype=int) == 1
    detected = np.diff(matrix.indptr)
    keep = in_tissue & (detected >= THRESHOLD)
    matrix = matrix[keep, :].tocsr()
    positions = positions.iloc[np.flatnonzero(keep)].copy()

    library = np.asarray(matrix.sum(axis=1)).ravel().astype(float)
    factors = np.divide(1e4, library, out=np.zeros_like(library), where=library > 0)
    log_matrix = matrix.multiply(factors[:, None]).tocsr()
    log_matrix.data = np.log1p(log_matrix.data)

    marker_scores: dict[str, np.ndarray] = {}
    for marker_name, marker_genes in spatial.MARKERS.items():
        marker_scores[marker_name], _ = spatial.marker_score(log_matrix, names, marker_genes)

    rows = positions["array_row"].to_numpy(dtype=int)
    cols = positions["array_col"].to_numpy(dtype=int)
    masks = spatial.build_masks(marker_scores, rows, cols)
    distance = spatial.distance_from_origins(rows, cols, masks["liver_side"])
    connected_tumour = masks["tumour_side"] & np.isfinite(distance)
    near = connected_tumour & (distance >= 1) & (distance <= 5)
    deep = connected_tumour & (distance >= 6)

    region = np.full(len(positions), "Other retained tissue", dtype=object)
    region[masks["liver_side"]] = "Liver-like"
    region[masks["tumour_side"]] = "Other tumour-side"
    region[deep] = "Deep tumour"
    region[near] = "Near-interface tumour"

    with (outs / "spatial" / "scalefactors_json.json").open(encoding="utf-8") as handle:
        scalefactors = json.load(handle)
    scale = float(scalefactors["tissue_hires_scalef"])
    pixels_per_micrometre = float(scalefactors["spot_diameter_fullres"]) * scale / 55.0
    image_path = outs / "spatial" / "tissue_hires_image.png"

    display = pd.DataFrame(
        {
            "sample": sample,
            "patient": str(row["patient"]),
            "hgp": str(row["hgp"]),
            "barcode": positions.index.astype(str),
            "x_hires": positions["pxl_col_in_fullres"].to_numpy(dtype=float) * scale,
            "y_hires": positions["pxl_row_in_fullres"].to_numpy(dtype=float) * scale,
            "array_row": rows,
            "array_col": cols,
            "region": region,
            "detected_genes": detected[keep],
            "library_size": library,
            "image_path": str(image_path),
            "pixels_per_micrometre": pixels_per_micrometre,
        }
    )
    summary = {
        "sample": sample,
        "patient": str(row["patient"]),
        "hgp": str(row["hgp"]),
        "image_path": str(image_path),
        "pixels_per_micrometre": pixels_per_micrometre,
        "n_retained_spots": int(len(display)),
        "n_liver_like": int(masks["liver_side"].sum()),
        "n_near_interface_tumour": int(near.sum()),
        "n_deep_tumour": int(deep.sum()),
        "n_other_tumour_side": int((masks["tumour_side"] & ~near & ~deep).sum()),
    }
    return display, summary


def main() -> None:
    primary = load_primary_samples()
    selected = select_representative_samples(primary)
    displays_by_sample: dict[str, pd.DataFrame] = {}
    summaries: list[dict[str, object]] = []
    for _, row in primary.iterrows():
        display, summary = prepare_sample(row)
        displays_by_sample[str(row["sample"])] = display
        summaries.append(summary)

    OUTPUT.mkdir(parents=True, exist_ok=True)
    all_displays = pd.concat(displays_by_sample.values(), ignore_index=True)
    all_displays.to_csv(
        OUTPUT / "all_spatial_spots.tsv.gz",
        sep="\t",
        index=False,
        encoding="utf-8",
        compression="gzip",
    )
    pd.DataFrame(summaries).to_csv(
        OUTPUT / "all_spatial_summary.tsv",
        sep="\t",
        index=False,
        encoding="utf-8",
    )

    representative_displays = pd.concat(
        [displays_by_sample[str(sample)] for sample in selected["sample"]],
        ignore_index=True,
    )
    representative_displays.to_csv(
        OUTPUT / "representative_spatial_spots.tsv.gz",
        sep="\t",
        index=False,
        encoding="utf-8",
        compression="gzip",
    )
    selected.to_csv(
        OUTPUT / "representative_sample_selection.tsv",
        sep="\t",
        index=False,
        encoding="utf-8",
    )
    pd.DataFrame(summaries).loc[
        lambda frame: frame["sample"].isin(selected["sample"])
    ].to_csv(
        OUTPUT / "representative_spatial_summary.tsv",
        sep="\t",
        index=False,
        encoding="utf-8",
    )
    (OUTPUT / "README.md").write_text(
        "# FFPE spatial display data\n\n"
        "The spot classes reproduce the primary threshold and expression-defined "
        "region construction from `analyze_rec_program_hgp_spatial_projection.py`. "
        "All six specimens are exported for the supplementary cohort overview. "
        "The representative subset is used to anchor Figure 4 in measured tissue space; all "
        "inference remains patient-level and is reported in the existing analysis "
        "contracts. Samples were selected as within-HGP medians for both primary "
        "near-interface programme summaries.\n",
        encoding="utf-8",
    )
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == "__main__":
    main()
