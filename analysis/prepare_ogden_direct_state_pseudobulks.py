#!/usr/bin/env python3
"""
Analysis: Prepare patient-state pseudobulks for direct Ogden epithelial contrasts.
Date: 2026-09-01
Random seed: 42 (no stochastic computation in this script)
Python: recorded at runtime
Key packages: anndata, numpy, pandas, scipy
"""

from __future__ import annotations

import hashlib
import json
import platform
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


np.random.seed(42)

ROOT = Path(__file__).resolve().parents[1]
INPUT = (
    ROOT
    / "data_sources"
    / "Ogden_2025_CRLM_multiome"
    / "CRCLM_multiome_GEX_decontaminated.h5ad"
)
SPECIFICATION = (
    ROOT / "metadata" / "deep_biology_phase2_analysis_spec_2026-09-01.md"
)
OUT = (
    ROOT
    / "analysis_results"
    / "deep_biology_upgrade_2026-08-31"
    / "phase2_mechanistic_specificity"
)

ANCHOR_STATE = "REC"
COMPARATORS = (
    "Hypoxia",
    "UPR",
    "iREC",
    "TA1",
    "Colonocyte",
    "Stem NOTUM",
    "Intermediate",
    "Goblet",
)
MIN_CELLS = 20
MIN_PAIRED_PATIENTS = 5


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv_gz(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, compression="gzip")


def main() -> None:
    if not INPUT.exists() or not SPECIFICATION.exists():
        raise FileNotFoundError("Missing Ogden AnnData or Phase 2 specification")

    adata = ad.read_h5ad(INPUT)
    required_obs = {"Patient", "Therapy", "Cell_type", "Cell_subtype"}
    missing_obs = sorted(required_obs.difference(adata.obs.columns))
    if missing_obs:
        raise ValueError(f"Missing observation fields: {missing_obs}")
    if adata.var_names.has_duplicates:
        raise ValueError("Ogden AnnData has duplicate gene identifiers")
    if np.any(pd.isna(adata.var_names)):
        raise ValueError("Ogden AnnData contains missing gene identifiers")

    epithelial_mask = adata.obs["Cell_type"].astype("string").eq("Epithelial").to_numpy()
    epithelial = adata[epithelial_mask]
    obs = epithelial.obs.loc[:, ["Patient", "Therapy", "Cell_subtype"]].copy()
    obs["Patient"] = obs["Patient"].astype("string")
    obs["Therapy"] = obs["Therapy"].astype("string")
    obs["Cell_subtype"] = obs["Cell_subtype"].astype("string")

    selected_states = {ANCHOR_STATE, *COMPARATORS}
    keep = obs["Cell_subtype"].isin(selected_states).to_numpy()
    obs = obs.loc[keep].reset_index(drop=True)
    matrix = epithelial.X[keep]
    if not sparse.issparse(matrix):
        matrix = sparse.csr_matrix(np.asarray(matrix, dtype=np.float64))
    else:
        matrix = sparse.csr_matrix(matrix, dtype=np.float64)
    if matrix.shape[0] != len(obs):
        raise ValueError("Expression matrix and metadata row counts differ")
    if matrix.data.size and (not np.isfinite(matrix.data).all() or (matrix.data < 0).any()):
        raise ValueError("Expression values must be finite and non-negative")

    groups = (
        obs.loc[:, ["Patient", "Cell_subtype"]]
        .drop_duplicates()
        .sort_values(["Patient", "Cell_subtype"])
        .itertuples(index=False, name=None)
    )
    metadata_rows: list[dict[str, object]] = []
    count_rows: list[np.ndarray] = []
    for patient, state in groups:
        selector = (
            obs["Patient"].eq(patient).to_numpy()
            & obs["Cell_subtype"].eq(state).to_numpy()
        )
        n_cells = int(selector.sum())
        summed = np.asarray(matrix[selector].sum(axis=0)).ravel().astype(np.float64)
        therapy_values = obs.loc[selector, "Therapy"].dropna().unique()
        if len(therapy_values) != 1:
            raise ValueError(f"Patient {patient} has inconsistent therapy metadata")
        sample_id = f"{patient}__{str(state).replace(' ', '_')}"
        metadata_rows.append(
            {
                "sample_id": sample_id,
                "patient": str(patient),
                "state": str(state),
                "therapy": str(therapy_values[0]),
                "n_cells": n_cells,
                "eligible_min_cells": n_cells >= MIN_CELLS,
                "library_size": float(summed.sum()),
            }
        )
        count_rows.append(summed)

    metadata = pd.DataFrame(metadata_rows)
    counts = np.vstack(count_rows)
    if (metadata["library_size"] <= 0).any():
        raise ValueError("One or more patient-state pseudobulks has zero library size")
    genes = pd.Index(epithelial.var_names.astype(str).str.upper(), name="gene")
    count_frame = pd.DataFrame(counts.T, columns=metadata["sample_id"])
    count_frame.insert(0, "gene", genes)

    eligibility_rows: list[dict[str, object]] = []
    for comparator in COMPARATORS:
        anchor_patients = set(
            metadata.loc[
                metadata["state"].eq(ANCHOR_STATE) & metadata["eligible_min_cells"],
                "patient",
            ]
        )
        comparator_patients = set(
            metadata.loc[
                metadata["state"].eq(comparator) & metadata["eligible_min_cells"],
                "patient",
            ]
        )
        paired = sorted(anchor_patients.intersection(comparator_patients))
        eligibility_rows.append(
            {
                "anchor_state": ANCHOR_STATE,
                "comparator_state": comparator,
                "minimum_cells_per_patient_state": MIN_CELLS,
                "n_paired_patients": len(paired),
                "eligible_minimum_five_patients": len(paired) >= MIN_PAIRED_PATIENTS,
                "paired_patient_ids": ";".join(paired),
            }
        )
    eligibility = pd.DataFrame(eligibility_rows)

    inventory = pd.DataFrame(
        [
            {
                "source_cells": adata.n_obs,
                "source_genes": adata.n_vars,
                "epithelial_cells": int(epithelial_mask.sum()),
                "selected_state_cells": len(obs),
                "patients": metadata["patient"].nunique(),
                "patient_state_pseudobulks": len(metadata),
                "missing_expression_values": 0,
                "direct_identifier_fields_used": 0,
                "analysis_unit": "pseudonymous patient",
            }
        ]
    )

    OUT.mkdir(parents=True, exist_ok=True)
    write_tsv_gz(count_frame, OUT / "ogden_direct_state_pseudobulk_counts.tsv.gz")
    metadata.to_csv(
        OUT / "ogden_direct_state_pseudobulk_metadata.tsv", sep="\t", index=False
    )
    eligibility.to_csv(
        OUT / "ogden_direct_contrast_eligibility.tsv", sep="\t", index=False
    )
    inventory.to_csv(OUT / "ogden_data_inventory.tsv", sep="\t", index=False)

    provenance = {
        "analysis_date": "2026-09-01",
        "python_version": platform.python_version(),
        "anndata_version": version("anndata"),
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "scipy_version": version("scipy"),
        "minimum_cells_per_patient_state": MIN_CELLS,
        "minimum_paired_patients": MIN_PAIRED_PATIENTS,
        "input_sha256": sha256(INPUT),
        "specification_sha256": sha256(SPECIFICATION),
    }
    with (OUT / "ogden_pseudobulk_provenance.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(provenance, handle, ensure_ascii=False, indent=2)

    print(inventory.to_string(index=False))
    print(eligibility.to_string(index=False))


if __name__ == "__main__":
    main()

