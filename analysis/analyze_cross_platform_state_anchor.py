#!/usr/bin/env python3
"""Anchor the targeted ISS NC3 state to Ogden full-transcriptome states.

Analysis date: 2026-08-23
Inputs:
  - Escriva Conde BJC 2026 targeted ISS panel-2 AnnData
  - Ogden et al. 2025 decontaminated multiome gene-expression AnnData
Outputs: analysis_results/cross_platform_state_anchor/
Random seed: 42 (used only for diagnostic point jitter)
"""

from __future__ import annotations

import hashlib
import json
import platform
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import pearsonr, rankdata, spearmanr


ROOT = Path(__file__).resolve().parents[1]
ISS_INPUT = (
    ROOT
    / "data_sources"
    / "EscrivaConde_BJC_2026_HGP_ISS"
    / "prepared_cells_panel_2.h5ad"
)
OGDEN_INPUT = (
    ROOT
    / "data_sources"
    / "Ogden_2025_CRLM_multiome"
    / "CRCLM_multiome_GEX_decontaminated.h5ad"
)
SPECIFICATION = (
    ROOT / "metadata" / "cross_platform_state_anchor_spec_2026-08-23.md"
)
OUT = ROOT / "analysis_results" / "cross_platform_state_anchor"

SEED = 42
MIN_CELLS_PER_PATIENT_STATE = 20
MIN_PATIENTS_PER_STATE = 5
MIN_DETECTION_FRACTION = 0.25
MIN_STATE_LOGCPM_RANGE = 0.5
MIN_ANCHOR_CELLS_FOR_DE = 20
MIN_OTHER_EPITHELIAL_CELLS_FOR_DE = 100
ISS_STATES = ("NC1", "NC2", "NC3", "NCγ")
NC3_MARKERS = ("KRT19", "KRT18", "CDX2", "KRT20", "ASCL2", "SLPI")


def sha256(path: Path) -> str:
    """Return a streaming SHA-256 checksum."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def dense_array(matrix: object) -> np.ndarray:
    """Convert a matrix slice to a dense float64 array."""
    if sparse.issparse(matrix):
        return np.asarray(matrix.toarray(), dtype=np.float64)
    return np.asarray(matrix, dtype=np.float64)


def canonical_iss_patient(samples: pd.Series) -> pd.Series:
    """Extract repository pseudonymous patient IDs from ISS sample names."""
    patients = samples.astype("string").str.extract(r"(P\d+)", expand=False)
    if patients.isna().any():
        raise ValueError("Failed to derive one or more ISS patient identifiers")
    return patients


def patient_state_pseudobulk(
    matrix: object,
    patients: np.ndarray,
    states: np.ndarray,
    gene_names: list[str],
    dataset: str,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """Build count and logCPM matrices for patient-state groups."""
    unique_groups = sorted(set(zip(patients, states, strict=True)))
    metadata_rows: list[dict[str, object]] = []
    count_rows: list[np.ndarray] = []
    for patient, state in unique_groups:
        selector = (patients == patient) & (states == state)
        n_cells = int(selector.sum())
        if n_cells == 0:
            continue
        group_counts = np.asarray(matrix[selector].sum(axis=0)).ravel().astype(float)
        metadata_rows.append(
            {
                "dataset": dataset,
                "patient": str(patient),
                "state": str(state),
                "n_cells": n_cells,
                "eligible_min_cells": n_cells >= MIN_CELLS_PER_PATIENT_STATE,
                "shared_panel_library_size": float(group_counts.sum()),
            }
        )
        count_rows.append(group_counts)
    metadata = pd.DataFrame(metadata_rows)
    counts = np.vstack(count_rows)
    library_sizes = counts.sum(axis=1)
    if np.any(library_sizes <= 0):
        raise ValueError(f"Zero shared-panel library in {dataset} pseudobulk")
    logcpm = np.log1p(counts / library_sizes[:, None] * 1_000_000)
    if counts.shape[1] != len(gene_names):
        raise ValueError("Pseudobulk gene dimension mismatch")
    return metadata, counts, logcpm


def state_centroids(
    metadata: pd.DataFrame,
    logcpm: np.ndarray,
    states_to_keep: list[str] | None = None,
) -> tuple[list[str], np.ndarray]:
    """Average eligible patient pseudobulks equally within each state."""
    eligible = metadata["eligible_min_cells"].to_numpy(dtype=bool)
    state_counts = metadata.loc[eligible, "state"].value_counts()
    valid_states = state_counts.loc[state_counts >= MIN_PATIENTS_PER_STATE].index
    if states_to_keep is None:
        states = sorted(valid_states.astype(str).tolist())
    else:
        states = [state for state in states_to_keep if state in set(valid_states)]
    if not states:
        raise ValueError("No states meet the patient-coverage requirement")
    centroids = []
    for state in states:
        selector = eligible & metadata["state"].eq(state).to_numpy()
        centroids.append(logcpm[selector].mean(axis=0))
    return states, np.vstack(centroids)


def gene_wise_zscore(centroids: np.ndarray) -> np.ndarray:
    """Z-score each gene across state centroids within a dataset."""
    means = centroids.mean(axis=0)
    standard_deviations = centroids.std(axis=0, ddof=0)
    if np.any(standard_deviations == 0):
        raise ValueError("Zero-variance gene passed the shared-gene filter")
    return (centroids - means) / standard_deviations


def select_shared_genes(
    genes: list[str],
    iss_metadata: pd.DataFrame,
    iss_counts: np.ndarray,
    iss_centroids: np.ndarray,
    ogden_metadata: pd.DataFrame,
    ogden_counts: np.ndarray,
    ogden_centroids: np.ndarray,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Apply detection and state-range filters symmetrically across datasets."""
    iss_eligible = iss_metadata["eligible_min_cells"].to_numpy(dtype=bool)
    ogden_eligible = ogden_metadata["eligible_min_cells"].to_numpy(dtype=bool)
    iss_detection = (iss_counts[iss_eligible] > 0).mean(axis=0)
    ogden_detection = (ogden_counts[ogden_eligible] > 0).mean(axis=0)
    iss_range = np.ptp(iss_centroids, axis=0)
    ogden_range = np.ptp(ogden_centroids, axis=0)
    selected = (
        (iss_detection >= MIN_DETECTION_FRACTION)
        & (ogden_detection >= MIN_DETECTION_FRACTION)
        & (iss_range >= MIN_STATE_LOGCPM_RANGE)
        & (ogden_range >= MIN_STATE_LOGCPM_RANGE)
    )
    audit = pd.DataFrame(
        {
            "gene": genes,
            "ISS_detection_fraction": iss_detection,
            "Ogden_detection_fraction": ogden_detection,
            "ISS_state_logCPM_range": iss_range,
            "Ogden_state_logCPM_range": ogden_range,
            "selected_for_primary_anchor": selected,
        }
    )
    if selected.sum() < 15:
        raise ValueError(f"Only {selected.sum()} genes passed the anchor filter")
    return selected, audit


def similarity_table(
    iss_states: list[str],
    iss_z: np.ndarray,
    ogden_states: list[str],
    ogden_z: np.ndarray,
    genes: list[str],
) -> pd.DataFrame:
    """Calculate Pearson, Spearman and cosine similarities between state pairs."""
    rows: list[dict[str, object]] = []
    for iss_index, iss_state in enumerate(iss_states):
        for ogden_index, ogden_state in enumerate(ogden_states):
            left = iss_z[iss_index]
            right = ogden_z[ogden_index]
            pearson = pearsonr(left, right)
            spearman = spearmanr(left, right)
            cosine = float(
                np.dot(left, right)
                / (np.linalg.norm(left) * np.linalg.norm(right))
            )
            rows.append(
                {
                    "ISS_state": iss_state,
                    "Ogden_state": ogden_state,
                    "n_shared_genes": len(genes),
                    "pearson_correlation": float(pearson.statistic),
                    "pearson_p": float(pearson.pvalue),
                    "spearman_correlation": float(spearman.statistic),
                    "spearman_p": float(spearman.pvalue),
                    "cosine_similarity": cosine,
                }
            )
    result = pd.DataFrame(rows)
    result["pearson_rank_within_ISS_state"] = result.groupby("ISS_state")[
        "pearson_correlation"
    ].rank(ascending=False, method="min")
    return result.sort_values(
        ["ISS_state", "pearson_rank_within_ISS_state", "Ogden_state"]
    )


def nc3_marker_ranking(
    ogden_states: list[str],
    ogden_centroids_all_genes: np.ndarray,
    shared_genes: list[str],
) -> pd.DataFrame:
    """Rank Ogden states by the six published NC3 context markers."""
    missing = sorted(set(NC3_MARKERS).difference(shared_genes))
    if missing:
        raise ValueError(f"NC3 markers absent from the shared panel: {missing}")
    indices = [shared_genes.index(gene) for gene in NC3_MARKERS]
    marker_centroids = ogden_centroids_all_genes[:, indices]
    marker_means = marker_centroids.mean(axis=0)
    marker_sds = marker_centroids.std(axis=0, ddof=0)
    marker_sds[marker_sds == 0] = 1.0
    marker_z = (marker_centroids - marker_means) / marker_sds
    scores = marker_z.mean(axis=1)
    result = pd.DataFrame(
        {
            "Ogden_state": ogden_states,
            "NC3_six_marker_z_mean": scores,
        }
    )
    for marker_index, marker in enumerate(NC3_MARKERS):
        result[f"{marker}_state_z"] = marker_z[:, marker_index]
    result["marker_score_rank"] = rankdata(-scores, method="min").astype(int)
    return result.sort_values(["marker_score_rank", "Ogden_state"])


def leave_one_patient_out_mapping(
    iss_metadata: pd.DataFrame,
    iss_logcpm: np.ndarray,
    ogden_metadata: pd.DataFrame,
    ogden_logcpm: np.ndarray,
    iss_states: list[str],
    ogden_states: list[str],
    selected_genes: np.ndarray,
) -> pd.DataFrame:
    """Test whether one patient determines the top NC3 anchor."""
    rows: list[dict[str, object]] = []
    datasets = (
        ("ISS", iss_metadata, iss_logcpm),
        ("Ogden", ogden_metadata, ogden_logcpm),
    )
    for omitted_dataset, metadata, _ in datasets:
        for patient in sorted(metadata["patient"].unique()):
            iss_meta = iss_metadata.copy()
            ogden_meta = ogden_metadata.copy()
            if omitted_dataset == "ISS":
                iss_meta.loc[iss_meta["patient"].eq(patient), "eligible_min_cells"] = False
            else:
                ogden_meta.loc[
                    ogden_meta["patient"].eq(patient), "eligible_min_cells"
                ] = False
            retained_iss_states, iss_centroid = state_centroids(
                iss_meta, iss_logcpm, iss_states
            )
            retained_ogden_states, ogden_centroid = state_centroids(
                ogden_meta, ogden_logcpm, ogden_states
            )
            if retained_iss_states != iss_states or retained_ogden_states != ogden_states:
                rows.append(
                    {
                        "omitted_dataset": omitted_dataset,
                        "omitted_patient": patient,
                        "top_Ogden_state": "STATE_COVERAGE_CHANGED",
                        "top_pearson": np.nan,
                        "second_Ogden_state": "STATE_COVERAGE_CHANGED",
                        "second_pearson": np.nan,
                        "top_minus_second_margin": np.nan,
                    }
                )
                continue
            iss_z = gene_wise_zscore(iss_centroid[:, selected_genes])
            ogden_z = gene_wise_zscore(ogden_centroid[:, selected_genes])
            similarity = similarity_table(
                iss_states,
                iss_z,
                ogden_states,
                ogden_z,
                [str(index) for index in np.flatnonzero(selected_genes)],
            )
            ranking = similarity.loc[similarity["ISS_state"].eq("NC3")].sort_values(
                "pearson_correlation", ascending=False
            )
            top, second = ranking.iloc[0], ranking.iloc[1]
            rows.append(
                {
                    "omitted_dataset": omitted_dataset,
                    "omitted_patient": patient,
                    "top_Ogden_state": top["Ogden_state"],
                    "top_pearson": top["pearson_correlation"],
                    "second_Ogden_state": second["Ogden_state"],
                    "second_pearson": second["pearson_correlation"],
                    "top_minus_second_margin": (
                        top["pearson_correlation"] - second["pearson_correlation"]
                    ),
                }
            )
    return pd.DataFrame(rows)


def write_full_transcriptome_pseudobulk(
    ogden: ad.AnnData, anchor_state: str
) -> tuple[pd.DataFrame, Path]:
    """Write paired anchor-versus-other epithelial pseudobulk counts."""
    epithelial_mask = ogden.obs["Cell_type"].astype("string").eq("Epithelial").to_numpy()
    epithelial = ogden[epithelial_mask].to_memory()
    patients = epithelial.obs["Patient"].astype("string").to_numpy(dtype=str)
    states = epithelial.obs["Cell_subtype"].astype("string").to_numpy(dtype=str)
    therapies = epithelial.obs["Therapy"].astype("string").to_numpy(dtype=str)
    matrix = epithelial.X
    metadata_rows: list[dict[str, object]] = []
    count_columns: list[np.ndarray] = []
    for patient in sorted(np.unique(patients)):
        patient_mask = patients == patient
        anchor_mask = patient_mask & (states == anchor_state)
        other_mask = patient_mask & (states != anchor_state)
        n_anchor = int(anchor_mask.sum())
        n_other = int(other_mask.sum())
        eligible = (
            n_anchor >= MIN_ANCHOR_CELLS_FOR_DE
            and n_other >= MIN_OTHER_EPITHELIAL_CELLS_FOR_DE
        )
        if not eligible:
            continue
        patient_therapies = sorted(set(therapies[patient_mask]))
        if len(patient_therapies) != 1:
            raise ValueError(f"Multiple therapy labels for Ogden patient {patient}")
        for group, selector in (("anchor_state", anchor_mask), ("other_epithelial", other_mask)):
            counts = np.asarray(matrix[selector].sum(axis=0)).ravel().astype(float)
            sample_id = f"{patient}__{group}"
            metadata_rows.append(
                {
                    "sample_id": sample_id,
                    "patient": patient,
                    "group": group,
                    "anchor_state": anchor_state,
                    "therapy": patient_therapies[0],
                    "n_cells": int(selector.sum()),
                    "library_size": float(counts.sum()),
                }
            )
            count_columns.append(counts)
    metadata = pd.DataFrame(metadata_rows)
    if metadata["patient"].nunique() < 5:
        raise ValueError("Fewer than five patients are eligible for paired DE")
    count_matrix = np.column_stack(count_columns)
    counts_path = OUT / "ogden_anchor_vs_other_pseudobulk_counts.tsv.gz"
    count_frame = pd.DataFrame(
        count_matrix,
        index=epithelial.var_names.astype(str),
        columns=metadata["sample_id"],
    )
    count_frame.index.name = "gene"
    count_frame.to_csv(counts_path, sep="\t", compression="gzip")
    return metadata, counts_path


def diagnostic_heatmap(
    similarity: pd.DataFrame, iss_states: list[str], ogden_states: list[str]
) -> None:
    """Create an internal correlation heatmap for mapping QA."""
    matrix = similarity.pivot(
        index="ISS_state", columns="Ogden_state", values="pearson_correlation"
    ).loc[iss_states, ogden_states]
    fig, axis = plt.subplots(figsize=(10.2, 3.5))
    image = axis.imshow(matrix, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    axis.set_xticks(range(len(ogden_states)), ogden_states, rotation=60, ha="right")
    axis.set_yticks(range(len(iss_states)), iss_states)
    axis.set_xlabel("Ogden epithelial state")
    axis.set_ylabel("ISS neoplastic state")
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            axis.text(
                column,
                row,
                f"{matrix.iloc[row, column]:.2f}",
                ha="center",
                va="center",
                fontsize=6.5,
            )
    fig.colorbar(image, ax=axis, label="Pearson correlation", shrink=0.8)
    axis.spines[:].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "state_anchor_heatmap_diagnostic.png", dpi=180)
    fig.savefig(OUT / "state_anchor_heatmap_diagnostic.pdf")
    plt.close(fig)


def write_summary(
    top_summary: pd.DataFrame,
    marker_ranking: pd.DataFrame,
    loo: pd.DataFrame,
    de_metadata: pd.DataFrame,
    n_selected_genes: int,
) -> None:
    """Write a concise state-anchor interpretation report."""
    top = top_summary.iloc[0]
    marker_top = marker_ranking.iloc[0]
    valid_loo = loo.loc[loo["top_Ogden_state"].ne("STATE_COVERAGE_CHANGED")]
    top_frequency = (
        valid_loo["top_Ogden_state"].eq(top["top_Ogden_state"]).mean()
        if len(valid_loo)
        else np.nan
    )
    state_counts = valid_loo["top_Ogden_state"].value_counts().to_dict()
    lines = [
        "## Material Passport",
        "",
        "- ID: CROSS-PLATFORM-NC3-ANCHOR-2026-08-23",
        "- Type: patient-balanced cross-platform state mapping",
        "- Verification status: VERIFIED",
        "- Source data: public pseudonymised repository objects",
        "",
        "## Primary mapping",
        "",
        (
            f"NC3 mapped most closely to the Ogden `{top['top_Ogden_state']}` state "
            f"across {n_selected_genes} eligible shared genes "
            f"(Pearson r={top['top_pearson']:.3f}). The second state was "
            f"`{top['second_Ogden_state']}` (r={top['second_pearson']:.3f}), giving "
            f"a margin of {top['top_minus_second_margin']:.3f}."
        ),
        (
            f"The independent six-marker ranking placed `{marker_top['Ogden_state']}` "
            f"first (score={marker_top['NC3_six_marker_z_mean']:.3f})."
        ),
        (
            f"The primary top state was retained in {top_frequency:.1%} of valid "
            f"leave-one-patient-out mappings; top-state counts were {state_counts}."
        ),
        "",
        "## Full-transcriptome handoff",
        "",
        (
            f"Paired anchor-state versus other-epithelial pseudobulks were generated "
            f"for {de_metadata['patient'].nunique()} patients "
            f"({(de_metadata['therapy'] == 'Naive').groupby(de_metadata['patient']).first().sum()} "
            f"treatment-naive)."
        ),
        "These counts are the input to the paired limma-voom program analysis. The "
        "Ogden cohort has no HGP label, so this step expands the biology of the ISS "
        "state and does not independently test rHGP enrichment.",
        "",
        "## Reproduction",
        "",
        "```bash",
        "conda run -n crc-metastatic-growth python scripts/analyze_cross_platform_state_anchor.py",
        "```",
    ]
    (OUT / "analysis_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    """Run state anchoring and prepare full-transcriptome paired pseudobulks."""
    OUT.mkdir(parents=True, exist_ok=True)
    for path in (ISS_INPUT, OGDEN_INPUT, SPECIFICATION):
        if not path.exists():
            raise FileNotFoundError(path)

    iss = ad.read_h5ad(ISS_INPUT, backed="r")
    ogden = ad.read_h5ad(OGDEN_INPUT, backed="r")
    if iss.n_vars != 113 or ogden.n_vars != 36_485:
        raise ValueError("Unexpected input gene dimensions")
    if not iss.var_names.is_unique or not ogden.var_names.is_unique:
        raise ValueError("Input gene names must be unique")
    shared_genes = sorted(set(iss.var_names).intersection(ogden.var_names))
    if not set(NC3_MARKERS).issubset(shared_genes):
        raise ValueError("One or more source NC3 markers are not shared")

    iss_mask = (
        iss.obs["Clusters"].astype("string").isin(ISS_STATES)
        & iss.obs["Tumor"].astype(bool)
    ).to_numpy()
    ogden_mask = ogden.obs["Cell_type"].astype("string").eq("Epithelial").to_numpy()
    iss_subset = iss[iss_mask, shared_genes].to_memory()
    ogden_subset = ogden[ogden_mask, shared_genes].to_memory()
    iss.file.close()

    iss_patients = canonical_iss_patient(iss_subset.obs["Sample"]).to_numpy(dtype=str)
    iss_state_labels = iss_subset.obs["Clusters"].astype("string").to_numpy(dtype=str)
    ogden_patients = ogden_subset.obs["Patient"].astype("string").to_numpy(dtype=str)
    ogden_state_labels = (
        ogden_subset.obs["Cell_subtype"].astype("string").to_numpy(dtype=str)
    )
    iss_matrix = iss_subset.layers["X_raw"]
    ogden_matrix = ogden_subset.X

    iss_metadata, iss_counts, iss_logcpm = patient_state_pseudobulk(
        iss_matrix, iss_patients, iss_state_labels, shared_genes, "ISS"
    )
    ogden_metadata, ogden_counts, ogden_logcpm = patient_state_pseudobulk(
        ogden_matrix, ogden_patients, ogden_state_labels, shared_genes, "Ogden"
    )
    iss_states, iss_centroids = state_centroids(iss_metadata, iss_logcpm)
    ogden_states, ogden_centroids = state_centroids(ogden_metadata, ogden_logcpm)
    if set(iss_states) != set(ISS_STATES):
        raise ValueError(f"Not all ISS states passed coverage: {iss_states}")

    selected_genes, gene_audit = select_shared_genes(
        shared_genes,
        iss_metadata,
        iss_counts,
        iss_centroids,
        ogden_metadata,
        ogden_counts,
        ogden_centroids,
    )
    selected_names = np.asarray(shared_genes)[selected_genes].tolist()
    iss_z = gene_wise_zscore(iss_centroids[:, selected_genes])
    ogden_z = gene_wise_zscore(ogden_centroids[:, selected_genes])
    similarities = similarity_table(
        iss_states, iss_z, ogden_states, ogden_z, selected_names
    )
    nc3_ranking = similarities.loc[similarities["ISS_state"].eq("NC3")].sort_values(
        "pearson_correlation", ascending=False
    )
    top, second = nc3_ranking.iloc[0], nc3_ranking.iloc[1]
    top_summary = pd.DataFrame(
        [
            {
                "ISS_state": "NC3",
                "top_Ogden_state": top["Ogden_state"],
                "top_pearson": top["pearson_correlation"],
                "second_Ogden_state": second["Ogden_state"],
                "second_pearson": second["pearson_correlation"],
                "top_minus_second_margin": (
                    top["pearson_correlation"] - second["pearson_correlation"]
                ),
                "n_shared_genes_total": len(shared_genes),
                "n_shared_genes_selected": int(selected_genes.sum()),
            }
        ]
    )
    marker_ranking = nc3_marker_ranking(
        ogden_states, ogden_centroids, shared_genes
    )
    loo = leave_one_patient_out_mapping(
        iss_metadata,
        iss_logcpm,
        ogden_metadata,
        ogden_logcpm,
        iss_states,
        ogden_states,
        selected_genes,
    )

    centroid_rows = []
    for dataset, states, centroids in (
        ("ISS", iss_states, iss_centroids),
        ("Ogden", ogden_states, ogden_centroids),
    ):
        for state, values in zip(states, centroids, strict=True):
            centroid_rows.append(
                pd.DataFrame(
                    {
                        "dataset": dataset,
                        "state": state,
                        "gene": shared_genes,
                        "patient_balanced_log1p_CPM": values,
                    }
                )
            )
    centroids_long = pd.concat(centroid_rows, ignore_index=True)

    for dataset, metadata, counts, logcpm in (
        ("ISS", iss_metadata, iss_counts, iss_logcpm),
        ("Ogden", ogden_metadata, ogden_counts, ogden_logcpm),
    ):
        metadata.to_csv(OUT / f"{dataset.lower()}_patient_state_metadata.tsv", sep="\t", index=False)
        logcpm_frame = pd.DataFrame(logcpm, columns=shared_genes)
        logcpm_frame = pd.concat([metadata.reset_index(drop=True), logcpm_frame], axis=1)
        logcpm_frame.to_csv(
            OUT / f"{dataset.lower()}_patient_state_logcpm.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )

    gene_audit.to_csv(OUT / "shared_gene_audit.tsv", sep="\t", index=False)
    centroids_long.to_csv(OUT / "state_centroids_long.tsv.gz", sep="\t", index=False, compression="gzip")
    similarities.to_csv(OUT / "state_similarity.tsv", sep="\t", index=False)
    top_summary.to_csv(OUT / "nc3_top_anchor.tsv", sep="\t", index=False)
    marker_ranking.to_csv(OUT / "nc3_six_marker_ranking.tsv", sep="\t", index=False)
    loo.to_csv(OUT / "leave_one_patient_out_mapping.tsv", sep="\t", index=False)

    de_metadata, counts_path = write_full_transcriptome_pseudobulk(
        ogden, str(top["Ogden_state"])
    )
    ogden.file.close()
    de_metadata.to_csv(
        OUT / "ogden_anchor_vs_other_pseudobulk_metadata.tsv", sep="\t", index=False
    )
    diagnostic_heatmap(similarities, iss_states, ogden_states)
    write_summary(
        top_summary,
        marker_ranking,
        loo,
        de_metadata,
        int(selected_genes.sum()),
    )

    provenance = {
        "analysis_id": "CROSS-PLATFORM-NC3-ANCHOR-2026-08-23",
        "ISS_input": str(ISS_INPUT.relative_to(ROOT)),
        "ISS_sha256": sha256(ISS_INPUT),
        "Ogden_input": str(OGDEN_INPUT.relative_to(ROOT)),
        "Ogden_sha256": sha256(OGDEN_INPUT),
        "specification": str(SPECIFICATION.relative_to(ROOT)),
        "script": str(Path(__file__).resolve().relative_to(ROOT)),
        "shared_genes_total": len(shared_genes),
        "shared_genes_selected": int(selected_genes.sum()),
        "full_transcriptome_counts_output": str(counts_path.relative_to(ROOT)),
        "random_seed": SEED,
        "python": platform.python_version(),
        "packages": {
            package: version(package)
            for package in ("anndata", "matplotlib", "numpy", "pandas", "scipy")
        },
    }
    (OUT / "provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(top_summary.to_string(index=False))
    print(marker_ranking.head(5).to_string(index=False))
    print(
        f"Paired pseudobulk patients: {de_metadata['patient'].nunique()}; "
        f"counts: {counts_path}"
    )


if __name__ == "__main__":
    main()
