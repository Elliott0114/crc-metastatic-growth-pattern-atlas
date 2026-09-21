#!/usr/bin/env python3
"""Project the frozen Ogden REC program into E-MTAB-12022 epithelial cells."""

from __future__ import annotations

import itertools
import math
import platform
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
H5AD = ROOT / "analysis_results" / "e_mtab_12022_primary_reconstruction.h5ad"
PROGRAM_FILE = ROOT / "analysis_results" / "ogden_anchor_program" / "selected_anchor_program_top50.tsv"
SPEC_FILE = ROOT / "metadata" / "rec_program_hgp_scrna_projection_spec_2026-08-23.md"
OUTPUT = ROOT / "analysis_results" / "rec_program_hgp_scrna_projection"
MIN_CELLS = 20
MIN_PROGRAM_GENES = 40
BALANCED_CELLS = 38
BALANCED_ITERATIONS = 1_000
BOOTSTRAP_REPLICATES = 10_000
SEED = 42
QC_SENSITIVITY_MIN_FEATURES = [150, 200, 300, 500]
QC_SENSITIVITY_MAX_MT = [0.20, 0.30, 0.50, 0.70]

MECHANISTIC_MODULES = {
    "hypoxia_mp6_overlap": [
        "ADM", "ANGPTL4", "ANKRD37", "EGLN3", "LDHA", "NDRG1", "NDUFA4L2", "SLC16A3",
    ],
    "regenerative_overlap": ["ADM", "APOL1", "DUOXA2", "ISG15", "OAS1", "SLC16A3"],
    "experimental_ap1_target_overlap": ["CYP3A5", "FHL2", "GSN", "KRT80", "PLAUR"],
    "nfkb_regulon_overlap": ["ABCG1", "BIRC3", "DUSP5", "MXD1", "PLAUR", "SDCBP2"],
}


def load_program() -> list[str]:
    table = pd.read_csv(PROGRAM_FILE, sep="\t").sort_values("program_rank")
    if len(table) != 50 or table["gene"].duplicated().any() or not (table["direction"] == "anchor_up").all():
        raise ValueError("Expected 50 unique REC-up program genes")
    return table["gene"].astype(str).str.upper().tolist()


def exact_permutation(values: np.ndarray, is_rhgp: np.ndarray) -> tuple[float, float, int]:
    observed = float(values[is_rhgp].mean() - values[~is_rhgp].mean())
    all_indices = set(range(len(values)))
    differences: list[float] = []
    for selected in itertools.combinations(range(len(values)), int(is_rhgp.sum())):
        rhgp_indices = sorted(selected)
        dhgp_indices = sorted(all_indices - set(selected))
        differences.append(float(values[rhgp_indices].mean() - values[dhgp_indices].mean()))
    null = np.asarray(differences)
    tolerance = 1e-12
    return (
        float(np.mean(null >= observed - tolerance)),
        float(np.mean(np.abs(null) >= abs(observed) - tolerance)),
        len(null),
    )


def hedges_g(rhgp_values: np.ndarray, dhgp_values: np.ndarray) -> float:
    pooled_variance = (
        (len(rhgp_values) - 1) * np.var(rhgp_values, ddof=1)
        + (len(dhgp_values) - 1) * np.var(dhgp_values, ddof=1)
    ) / (len(rhgp_values) + len(dhgp_values) - 2)
    if not np.isfinite(pooled_variance) or pooled_variance <= 0:
        return math.nan
    correction = 1 - 3 / (4 * (len(rhgp_values) + len(dhgp_values)) - 9)
    return float(
        correction
        * (rhgp_values.mean() - dhgp_values.mean())
        / math.sqrt(pooled_variance)
    )


def bootstrap_mean_difference(
    values: np.ndarray,
    is_rhgp: np.ndarray,
    seed: int,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    rhgp_values = values[is_rhgp]
    dhgp_values = values[~is_rhgp]
    rhgp_indices = rng.integers(
        0, len(rhgp_values), size=(BOOTSTRAP_REPLICATES, len(rhgp_values))
    )
    dhgp_indices = rng.integers(
        0, len(dhgp_values), size=(BOOTSTRAP_REPLICATES, len(dhgp_values))
    )
    differences = (
        rhgp_values[rhgp_indices].mean(axis=1) - dhgp_values[dhgp_indices].mean(axis=1)
    )
    lower, upper = np.quantile(differences, [0.025, 0.975])
    return float(lower), float(upper)


def standardize_score(log2_cpm: np.ndarray, indices: list[int]) -> np.ndarray:
    selected = log2_cpm[:, indices]
    means = selected.mean(axis=0)
    standard_deviations = selected.std(axis=0, ddof=0)
    variable = standard_deviations > 0
    if not np.any(variable):
        return np.zeros(log2_cpm.shape[0], dtype=float)
    return ((selected[:, variable] - means[variable]) / standard_deviations[variable]).mean(axis=1)


def summarize_endpoint(
    patient_table: pd.DataFrame,
    values: np.ndarray,
    endpoint: str,
    seed: int,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    is_rhgp = patient_table["hgp"].to_numpy() == "rHGP"
    observed = float(values[is_rhgp].mean() - values[~is_rhgp].mean())
    p_one_sided, p_two_sided, assignments = exact_permutation(values, is_rhgp)
    ci_lower, ci_upper = bootstrap_mean_difference(values, is_rhgp, seed)
    summary = {
        "endpoint": endpoint,
        "n_patients": len(values),
        "n_rhgp": int(is_rhgp.sum()),
        "n_dhgp": int((~is_rhgp).sum()),
        "rhgp_mean": float(values[is_rhgp].mean()),
        "dhgp_mean": float(values[~is_rhgp].mean()),
        "rhgp_minus_dhgp": observed,
        "bootstrap_ci_lower": ci_lower,
        "bootstrap_ci_upper": ci_upper,
        "hedges_g": hedges_g(values[is_rhgp], values[~is_rhgp]),
        "exact_p_one_sided": p_one_sided,
        "exact_p_two_sided": p_two_sided,
        "n_exact_assignments": assignments,
    }
    loo_rows: list[dict[str, object]] = []
    for omitted_index, omitted in patient_table.reset_index(drop=True).iterrows():
        keep = np.arange(len(values)) != omitted_index
        effect = float(values[keep & is_rhgp].mean() - values[keep & ~is_rhgp].mean())
        loo_rows.append({
            "endpoint": endpoint,
            "omitted_patient": omitted["patient"],
            "omitted_hgp": omitted["hgp"],
            "rhgp_minus_dhgp": effect,
            "positive": effect > 0,
        })
    summary["loo_positive_n"] = int(sum(row["positive"] for row in loo_rows))
    return summary, loo_rows


def balanced_resampling(
    program_matrix: sparse.csr_matrix,
    libraries: np.ndarray,
    patient_array: np.ndarray,
    patient_table: pd.DataFrame,
) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    is_rhgp = patient_table["hgp"].to_numpy() == "rHGP"
    pools = {
        patient: np.flatnonzero(patient_array == patient)
        for patient in patient_table["patient"].astype(str)
    }
    rows: list[dict[str, object]] = []
    for iteration in range(1, BALANCED_ITERATIONS + 1):
        counts: list[np.ndarray] = []
        total_libraries: list[float] = []
        for patient in patient_table["patient"].astype(str):
            selected = rng.choice(pools[patient], size=BALANCED_CELLS, replace=False)
            counts.append(
                np.asarray(program_matrix[selected, :].sum(axis=0)).ravel().astype(float)
            )
            total_libraries.append(float(libraries[selected].sum()))
        count_matrix = np.vstack(counts)
        library_array = np.asarray(total_libraries)
        log2_cpm = np.log2((count_matrix + 0.5) / (library_array[:, None] + 1.0) * 1e6)
        equal_gene_score = standardize_score(log2_cpm, list(range(log2_cpm.shape[1])))
        pooled_score = np.log2(
            (count_matrix.sum(axis=1) + 0.5) / (library_array + 1.0) * 1e6
        )
        for endpoint, values in (
            ("rec_program_score", equal_gene_score),
            ("pooled_program_log2_cpm", pooled_score),
        ):
            effect = float(values[is_rhgp].mean() - values[~is_rhgp].mean())
            rows.append({
                "iteration": iteration,
                "endpoint": endpoint,
                "balanced_cells_per_patient": BALANCED_CELLS,
                "rhgp_minus_dhgp": effect,
                "positive": effect > 0,
            })
    return pd.DataFrame(rows)


def qc_threshold_sensitivity(
    program_matrix: sparse.csr_matrix,
    libraries: np.ndarray,
    patient_array: np.ndarray,
    hgp_array: np.ndarray,
    n_features: np.ndarray,
    mt_fraction: np.ndarray,
    primary_patients: set[str],
) -> pd.DataFrame:
    """Rebuild patient pseudobulks after tightening the source-study QC rules."""
    rows: list[dict[str, object]] = []
    for minimum_features in QC_SENSITIVITY_MIN_FEATURES:
        for maximum_mt_fraction in QC_SENSITIVITY_MAX_MT:
            keep_cells = (
                (n_features >= minimum_features)
                & (mt_fraction <= maximum_mt_fraction)
            )
            cell_counts = pd.Series(patient_array[keep_cells]).value_counts()
            eligible_patients = sorted(
                str(patient)
                for patient, count in cell_counts.items()
                if int(count) >= MIN_CELLS
            )
            eligible_hgp = np.asarray([
                str(hgp_array[np.flatnonzero(patient_array == patient)[0]])
                for patient in eligible_patients
            ])
            n_rhgp = int(np.sum(eligible_hgp == "rHGP"))
            n_dhgp = int(np.sum(eligible_hgp == "dHGP"))
            common = {
                "minimum_features": minimum_features,
                "maximum_mt_fraction": maximum_mt_fraction,
                "n_eligible_patients": len(eligible_patients),
                "n_rhgp": n_rhgp,
                "n_dhgp": n_dhgp,
                "eligible_patients": ",".join(eligible_patients),
                "same_primary_patient_set": set(eligible_patients) == primary_patients,
                "minimum_eligible_cell_count": (
                    int(min(cell_counts[patient] for patient in eligible_patients))
                    if eligible_patients else 0
                ),
            }
            if n_rhgp == 0 or n_dhgp == 0:
                for endpoint in ["rec_program_score", "pooled_program_log2_cpm"]:
                    rows.append({
                        **common,
                        "endpoint": endpoint,
                        "rhgp_minus_dhgp": math.nan,
                        "exact_p_two_sided": math.nan,
                    })
                continue

            count_rows: list[np.ndarray] = []
            library_rows: list[float] = []
            for patient in eligible_patients:
                patient_cells = keep_cells & (patient_array == patient)
                count_rows.append(
                    np.asarray(program_matrix[patient_cells, :].sum(axis=0))
                    .ravel()
                    .astype(float)
                )
                library_rows.append(float(libraries[patient_cells].sum()))
            count_matrix = np.vstack(count_rows)
            library_values = np.asarray(library_rows)
            log2_cpm = np.log2(
                (count_matrix + 0.5) / (library_values[:, None] + 1.0) * 1e6
            )
            endpoint_values = {
                "rec_program_score": standardize_score(
                    log2_cpm, list(range(log2_cpm.shape[1]))
                ),
                "pooled_program_log2_cpm": np.log2(
                    (count_matrix.sum(axis=1) + 0.5)
                    / (library_values + 1.0)
                    * 1e6
                ),
            }
            is_rhgp = eligible_hgp == "rHGP"
            for endpoint, values in endpoint_values.items():
                effect = float(values[is_rhgp].mean() - values[~is_rhgp].mean())
                _, p_two_sided, _ = exact_permutation(values, is_rhgp)
                rows.append({
                    **common,
                    "endpoint": endpoint,
                    "rhgp_minus_dhgp": effect,
                    "exact_p_two_sided": p_two_sided,
                })
    return pd.DataFrame(rows)


def write_summary(
    coverage: pd.DataFrame,
    endpoint_summary: pd.DataFrame,
    balanced_summary: pd.DataFrame,
    decision: pd.DataFrame,
) -> None:
    endpoints = {row["endpoint"]: row for _, row in endpoint_summary.iterrows()}
    primary = endpoints["rec_program_score"]
    pooled = endpoints["pooled_program_log2_cpm"]
    balanced = {
        row["endpoint"]: row for _, row in balanced_summary.iterrows()
    }
    lines = [
        "# Fixed REC-program projection into E-MTAB-12022 epithelial cells",
        "",
        f"Decision: **{decision.iloc[0]['decision']}**.",
        "",
        "## Coverage and cohort",
        "",
        f"The reconstruction contained {int(coverage.iloc[0]['n_matched'])}/50 program genes. "
        f"Five patients met the fixed 20-cell threshold (three rHGP and two dHGP).",
        "",
        "## Primary result",
        "",
        f"The equal-gene REC-program score differed by {primary['rhgp_minus_dhgp']:+.3f} "
        f"(rHGP minus dHGP; bootstrap 95% CI {primary['bootstrap_ci_lower']:+.3f} to "
        f"{primary['bootstrap_ci_upper']:+.3f}; Hedges' g {primary['hedges_g']:+.3f}; "
        f"exact two-sided P={primary['exact_p_two_sided']:.3f}).",
        "",
        f"The pooled-program log2 CPM difference was {pooled['rhgp_minus_dhgp']:+.3f} "
        f"(bootstrap 95% CI {pooled['bootstrap_ci_lower']:+.3f} to "
        f"{pooled['bootstrap_ci_upper']:+.3f}; Hedges' g {pooled['hedges_g']:+.3f}; "
        f"exact two-sided P={pooled['exact_p_two_sided']:.3f}).",
        "",
        "## Cell-recovery sensitivity",
        "",
        f"After repeatedly balancing every patient to 38 epithelial cells, the HGP effect was positive in "
        f"{balanced['rec_program_score']['positive_fraction']:.1%} of iterations for the equal-gene score "
        f"and {balanced['pooled_program_log2_cpm']['positive_fraction']:.1%} for pooled program CPM.",
        "",
        "## Interpretation",
        "",
        str(decision.iloc[0]["interpretation"]),
        "",
        "The reconstructed single-cell and Visium datasets share a source study and some patients; this is "
        "an epithelial-source check across modalities, not an additional independent cohort.",
        "",
        "## Reproducibility",
        "",
        f"Specification: `{SPEC_FILE.relative_to(ROOT)}`. Seed: {SEED}; balanced iterations: "
        f"{BALANCED_ITERATIONS}; patient bootstrap replicates: {BOOTSTRAP_REPLICATES}.",
    ]
    (OUTPUT / "analysis_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    requested_program = load_program()
    adata = ad.read_h5ad(H5AD)
    epithelial_mask = adata.obs["cell_type"].astype(str).to_numpy() == "CRC/epithelial"
    epithelial_obs = adata.obs.loc[epithelial_mask].copy()
    epithelial_counts = sparse.csr_matrix(adata.X[epithelial_mask, :])
    patient_counts = (
        epithelial_obs.groupby(["patient", "hgp"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .sort_values("patient")
    )
    eligible = patient_counts[patient_counts["n_cells"] >= MIN_CELLS].reset_index(drop=True)
    if len(eligible) != 5 or eligible["hgp"].value_counts().to_dict() != {"rHGP": 3, "dHGP": 2}:
        raise ValueError("Expected five eligible patients: three rHGP and two dHGP")
    if int(eligible["n_cells"].min()) != BALANCED_CELLS:
        raise ValueError(f"Expected minimum eligible epithelial count of {BALANCED_CELLS}")

    var_names = np.asarray(pd.Index(adata.var_names.astype(str)).str.upper())
    matched_program = [gene for gene in requested_program if np.any(var_names == gene)]
    if len(matched_program) < MIN_PROGRAM_GENES:
        raise ValueError(f"Only {len(matched_program)} of 50 program genes matched")
    gene_indices = [int(np.flatnonzero(var_names == gene)[0]) for gene in matched_program]
    program_matrix = epithelial_counts[:, gene_indices].tocsr()
    cell_libraries = np.asarray(epithelial_counts.sum(axis=1)).ravel().astype(float)
    patient_array = epithelial_obs["patient"].astype(str).to_numpy()

    raw_count_rows: list[np.ndarray] = []
    library_sizes: list[float] = []
    for patient in eligible["patient"].astype(str):
        mask = patient_array == patient
        raw_count_rows.append(
            np.asarray(program_matrix[mask, :].sum(axis=0)).ravel().astype(float)
        )
        library_sizes.append(float(cell_libraries[mask].sum()))
    raw_counts = np.vstack(raw_count_rows)
    library_array = np.asarray(library_sizes)
    log2_cpm = np.log2((raw_counts + 0.5) / (library_array[:, None] + 1.0) * 1e6)

    score_table = eligible.copy()
    score_table["library_size"] = library_array
    score_table["matched_program_genes"] = len(matched_program)
    score_table["rec_program_score"] = standardize_score(
        log2_cpm, list(range(len(matched_program)))
    )
    score_table["program_raw_count"] = raw_counts.sum(axis=1)
    score_table["pooled_program_log2_cpm"] = np.log2(
        (score_table["program_raw_count"].to_numpy() + 0.5)
        / (library_array + 1.0)
        * 1e6
    )

    module_coverage_rows: list[dict[str, object]] = []
    for module_name, genes in MECHANISTIC_MODULES.items():
        matched = [gene for gene in genes if gene in matched_program]
        indices = [matched_program.index(gene) for gene in matched]
        score_table[module_name] = standardize_score(log2_cpm, indices)
        module_coverage_rows.append({
            "module": module_name,
            "n_defined": len(genes),
            "n_matched": len(matched),
            "matched_genes": ",".join(matched),
            "missing_genes": ",".join(gene for gene in genes if gene not in matched),
        })

    patient_gene_rows: list[dict[str, object]] = []
    for patient_index, patient_row in eligible.iterrows():
        for gene_index, gene in enumerate(matched_program):
            patient_gene_rows.append({
                "patient": patient_row["patient"],
                "hgp": patient_row["hgp"],
                "n_cells": int(patient_row["n_cells"]),
                "gene": gene,
                "raw_count": raw_counts[patient_index, gene_index],
                "library_size": library_array[patient_index],
                "log2_cpm": log2_cpm[patient_index, gene_index],
            })

    endpoint_names = [
        "rec_program_score",
        "pooled_program_log2_cpm",
        *MECHANISTIC_MODULES.keys(),
    ]
    endpoint_rows: list[dict[str, object]] = []
    loo_rows: list[dict[str, object]] = []
    for endpoint_index, endpoint in enumerate(endpoint_names):
        summary, endpoint_loo = summarize_endpoint(
            score_table,
            score_table[endpoint].to_numpy(dtype=float),
            endpoint,
            SEED + endpoint_index,
        )
        endpoint_rows.append(summary)
        loo_rows.extend(endpoint_loo)
    endpoint_summary = pd.DataFrame(endpoint_rows)
    loo = pd.DataFrame(loo_rows)

    balanced_iterations = balanced_resampling(
        program_matrix,
        cell_libraries,
        patient_array,
        eligible,
    )
    balanced_summary = (
        balanced_iterations.groupby("endpoint", as_index=False)
        .agg(
            iterations=("iteration", "nunique"),
            median_effect=("rhgp_minus_dhgp", "median"),
            mean_effect=("rhgp_minus_dhgp", "mean"),
            q025=("rhgp_minus_dhgp", lambda values: float(values.quantile(0.025))),
            q975=("rhgp_minus_dhgp", lambda values: float(values.quantile(0.975))),
            positive_fraction=("positive", "mean"),
        )
    )
    qc_sensitivity = qc_threshold_sensitivity(
        program_matrix=program_matrix,
        libraries=cell_libraries,
        patient_array=patient_array,
        hgp_array=epithelial_obs["hgp"].astype(str).to_numpy(),
        n_features=epithelial_obs["n_features"].to_numpy(dtype=int),
        mt_fraction=epithelial_obs["mt_fraction"].to_numpy(dtype=float),
        primary_patients=set(eligible["patient"].astype(str)),
    )

    primary = endpoint_summary[endpoint_summary["endpoint"] == "rec_program_score"].iloc[0]
    pooled = endpoint_summary[endpoint_summary["endpoint"] == "pooled_program_log2_cpm"].iloc[0]
    balanced_lookup = balanced_summary.set_index("endpoint")
    primary_stable = (
        primary["rhgp_minus_dhgp"] > 0
        and int(primary["loo_positive_n"]) == len(eligible)
        and balanced_lookup.loc["rec_program_score", "positive_fraction"] >= 0.8
    )
    pooled_stable = (
        pooled["rhgp_minus_dhgp"] > 0
        and int(pooled["loo_positive_n"]) == len(eligible)
        and balanced_lookup.loc["pooled_program_log2_cpm", "positive_fraction"] >= 0.8
    )
    if primary_stable and pooled_stable:
        label = "EPITHELIAL_PROGRAM_SUPPORT"
        interpretation = (
            "Both equal-gene and pooled summaries support higher REC-program expression in rHGP epithelial "
            "pseudobulks, including after equalizing cell recovery."
        )
    elif pooled_stable:
        label = "ABUNDANCE_WEIGHTED_EPITHELIAL_SUPPORT"
        interpretation = (
            "The abundance-weighted program is higher in rHGP epithelial pseudobulks, but the equal-gene "
            "summary is not sufficiently stable; the result supports an epithelial source without implying "
            "uniform activation of all program genes."
        )
    else:
        label = "NO_STABLE_EPITHELIAL_PROGRAM_SUPPORT"
        interpretation = (
            "The fixed REC program is not stably higher in rHGP reconstructed epithelial pseudobulks and "
            "cannot resolve the mixed-spot spatial signal."
        )
    decision = pd.DataFrame([{
        "decision": label,
        "eligible_patients": len(eligible),
        "n_rhgp": int((eligible["hgp"] == "rHGP").sum()),
        "n_dhgp": int((eligible["hgp"] == "dHGP").sum()),
        "matched_program_genes": len(matched_program),
        "primary_effect": primary["rhgp_minus_dhgp"],
        "primary_loo_positive_n": int(primary["loo_positive_n"]),
        "primary_balanced_positive_fraction": balanced_lookup.loc[
            "rec_program_score", "positive_fraction"
        ],
        "pooled_effect": pooled["rhgp_minus_dhgp"],
        "pooled_loo_positive_n": int(pooled["loo_positive_n"]),
        "pooled_balanced_positive_fraction": balanced_lookup.loc[
            "pooled_program_log2_cpm", "positive_fraction"
        ],
        "interpretation": interpretation,
    }])
    coverage = pd.DataFrame([{
        "n_defined": len(requested_program),
        "n_matched": len(matched_program),
        "coverage_fraction": len(matched_program) / len(requested_program),
        "matched_genes": ",".join(matched_program),
        "missing_genes": ",".join(gene for gene in requested_program if gene not in matched_program),
    }])
    run_info = pd.DataFrame([{
        "analysis_date": "2026-08-24",
        "python_version": platform.python_version(),
        "anndata_version": ad.__version__,
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "random_seed": SEED,
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "balanced_iterations": BALANCED_ITERATIONS,
        "balanced_cells_per_patient": BALANCED_CELLS,
        "qc_sensitivity_min_features": ",".join(
            str(value) for value in QC_SENSITIVITY_MIN_FEATURES
        ),
        "qc_sensitivity_max_mt_fraction": ",".join(
            str(value) for value in QC_SENSITIVITY_MAX_MT
        ),
        "statistical_unit": "patient",
    }])

    OUTPUT.mkdir(parents=True, exist_ok=True)
    outputs = {
        "program_coverage.tsv": coverage,
        "mechanistic_module_coverage.tsv": pd.DataFrame(module_coverage_rows),
        "patient_scores.tsv": score_table,
        "patient_gene_pseudobulk.tsv": pd.DataFrame(patient_gene_rows),
        "endpoint_summary.tsv": endpoint_summary,
        "leave_one_patient_out.tsv": loo,
        "balanced_iteration_effects.tsv.gz": balanced_iterations,
        "balanced_summary.tsv": balanced_summary,
        "qc_threshold_sensitivity.tsv": qc_sensitivity,
        "decision.tsv": decision,
        "run_info.tsv": run_info,
    }
    for filename, frame in outputs.items():
        compression = "gzip" if filename.endswith(".gz") else None
        frame.to_csv(
            OUTPUT / filename,
            sep="\t",
            index=False,
            encoding="utf-8",
            compression=compression,
        )
    write_summary(coverage, endpoint_summary, balanced_summary, decision)
    manifest_lines = [
        "# Analysis output manifest",
        "",
        "Generated by `scripts/analyze_rec_program_hgp_scrna_projection.py`.",
        "",
        "| File | Purpose |",
        "|---|---|",
        "| `analysis_summary.md` | Human-readable result and claim boundary |",
        "| `program_coverage.tsv` | Frozen full-program coverage |",
        "| `mechanistic_module_coverage.tsv` | Fixed interpretability-subset coverage |",
        "| `patient_scores.tsv` | Patient-level full and mechanistic scores |",
        "| `patient_gene_pseudobulk.tsv` | Auditable epithelial pseudobulk values |",
        "| `endpoint_summary.tsv` | Patient-level effects and exact tests |",
        "| `leave_one_patient_out.tsv` | Directional influence checks |",
        "| `balanced_iteration_effects.tsv.gz` | Equal-cell resampling effects |",
        "| `balanced_summary.tsv` | Equal-cell sensitivity summary |",
        "| `qc_threshold_sensitivity.tsv` | Patient-pseudobulk effects after tightening cell QC |",
        "| `decision.tsv` | Prespecified interpretation rule |",
        "| `run_info.tsv` | Software and reproducibility settings |",
    ]
    (OUTPUT / "_analysis_outputs.md").write_text(
        "\n".join(manifest_lines) + "\n", encoding="utf-8"
    )
    print(decision.to_string(index=False))
    print(endpoint_summary.to_string(index=False))


if __name__ == "__main__":
    main()
