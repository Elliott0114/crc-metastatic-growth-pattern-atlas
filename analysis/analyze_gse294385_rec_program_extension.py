#!/usr/bin/env python3
"""Test the fixed REC program across paired liver micro- and macrometastases."""

from __future__ import annotations

import gzip
import hashlib
import itertools
import math
import platform
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data_sources" / "Liu_2026_GSE294385"
EXTRACTED = SOURCE / "extracted"
ARCHIVES = SOURCE / "archives"
MANIFEST = ROOT / "metadata" / "gse294385_liver_sample_manifest.tsv"
ANNOTATION = SOURCE / "visium_liver_meta_after_qc.tsv.gz"
REC_PROGRAM = ROOT / "analysis_results" / "ogden_anchor_program" / "selected_anchor_program_top50.tsv"
TOP15_CORE = ROOT / "analysis_results" / "rec_program_cross_modal_concordance" / "top15_interpretable_core.tsv"
SPEC = ROOT / "metadata" / "gse294385_rec_program_extension_spec_2026-08-23.md"
OUTPUT = ROOT / "analysis_results" / "gse294385_rec_program_extension"
SEED = 42
BOOTSTRAP_REPLICATES = 10_000
BALANCED_ITERATIONS = 1_000
MIN_MICRO_SPOTS_SENSITIVITY = 20

REGION_MAP = {
    "Liver micrometastasis tumor": "micro_tumour",
    "Liver macrometastasis tumor": "macro_tumour",
    "Liver micrometastasis stroma": "micro_stroma",
    "Liver macrometastasis stroma": "macro_stroma",
}

MECHANISTIC_MODULES = {
    "hypoxia_mp6_overlap": [
        "ADM", "ANGPTL4", "ANKRD37", "EGLN3", "LDHA", "NDRG1", "NDUFA4L2", "SLC16A3",
    ],
    "regenerative_overlap": ["ADM", "APOL1", "DUOXA2", "ISG15", "OAS1", "SLC16A3"],
    "experimental_ap1_target_overlap": ["CYP3A5", "FHL2", "GSN", "KRT80", "PLAUR"],
    "nfkb_regulon_overlap": ["ABCG1", "BIRC3", "DUSP5", "MXD1", "PLAUR", "SDCBP2"],
}
MICROMETASTASIS_SIX_GENE = ["RNF40", "AEN", "WEE1", "BCL7B", "YME1L1", "COX17"]


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_programs() -> tuple[dict[str, list[str]], list[str]]:
    rec = pd.read_csv(REC_PROGRAM, sep="\t").sort_values("program_rank")
    if len(rec) != 50 or rec["gene"].duplicated().any():
        raise ValueError("Expected a unique 50-gene REC program")
    rec_genes = rec["gene"].astype(str).str.upper().tolist()
    top15 = pd.read_csv(TOP15_CORE, sep="\t").sort_values("cross_modal_support_rank")
    top15_genes = top15["gene"].astype(str).str.upper().tolist()
    if len(top15_genes) != 15 or not set(top15_genes).issubset(rec_genes):
        raise ValueError("Expected a 15-gene HGP core nested in the REC program")
    programs = {
        "rec_equal_gene": rec_genes,
        "hgp_cross_modal_top15": top15_genes,
        **MECHANISTIC_MODULES,
        "source_micrometastasis_six_gene": MICROMETASTASIS_SIX_GENE,
    }
    selected_genes = list(dict.fromkeys(gene for genes in programs.values() for gene in genes))
    return programs, selected_genes


def read_text_gzip(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle]


def read_sample(
    sample: str,
    patient: str,
    annotation: pd.DataFrame,
    selected_genes: list[str],
) -> tuple[pd.DataFrame, dict[str, object], set[str]]:
    matrix_root = EXTRACTED / sample / "filtered_feature_bc_matrix"
    barcodes = np.asarray(read_text_gzip(matrix_root / "barcodes.tsv.gz"), dtype=str)
    features = pd.read_csv(
        matrix_root / "features.tsv.gz",
        sep="\t",
        header=None,
        compression="gzip",
        dtype=str,
    )
    gene_names = features.iloc[:, 1].astype(str).str.upper().to_numpy()
    with gzip.open(matrix_root / "matrix.mtx.gz", "rb") as handle:
        matrix = sparse.csc_matrix(mmread(handle))
    if matrix.shape != (len(gene_names), len(barcodes)):
        raise ValueError(
            f"Unexpected matrix dimensions for {sample}: {matrix.shape}, "
            f"features={len(gene_names)}, barcodes={len(barcodes)}"
        )
    sample_annotation = annotation[annotation["sample"] == sample].copy()
    sample_annotation = sample_annotation.drop_duplicates("spot_barcode").set_index("spot_barcode")
    aligned_annotation = sample_annotation.reindex(barcodes)
    region = aligned_annotation["region"].to_numpy(dtype=object)
    keep = pd.notna(region)
    selected_indices = np.flatnonzero(keep)
    library = np.asarray(matrix[:, selected_indices].sum(axis=0)).ravel().astype(float)
    present = set(selected_genes).intersection(set(gene_names.tolist()))
    gene_counts = np.zeros((len(selected_indices), len(selected_genes)), dtype=np.float64)
    for gene_index, gene in enumerate(selected_genes):
        matching = np.flatnonzero(gene_names == gene)
        if len(matching):
            gene_counts[:, gene_index] = np.asarray(
                matrix[matching, :][:, selected_indices].sum(axis=0)
            ).ravel()
    frame = pd.DataFrame(gene_counts, columns=selected_genes)
    frame.insert(0, "library_size", library)
    frame.insert(0, "region", region[keep])
    frame.insert(0, "spot_barcode", barcodes[keep])
    frame.insert(0, "patient", patient)
    frame.insert(0, "sample", sample)
    qc = {
        "sample": sample,
        "patient": patient,
        "matrix_spots": len(barcodes),
        "source_annotated_spots": len(sample_annotation),
        "selected_region_spots": int(keep.sum()),
        "micro_tumour_spots": int(np.sum(region[keep] == "micro_tumour")),
        "macro_tumour_spots": int(np.sum(region[keep] == "macro_tumour")),
        "micro_stroma_spots": int(np.sum(region[keep] == "micro_stroma")),
        "macro_stroma_spots": int(np.sum(region[keep] == "macro_stroma")),
        "selected_genes_present": len(present),
    }
    return frame, qc, present


def pseudobulk_spots(
    spots: pd.DataFrame,
    group_columns: list[str],
    selected_genes: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    grouped = spots.groupby(group_columns, sort=True, observed=True)
    metadata = grouped.agg(
        n_spots=("spot_barcode", "size"),
        library_size=("library_size", "sum"),
    ).reset_index()
    counts = grouped[selected_genes].sum().reset_index()
    counts = counts.merge(metadata, on=group_columns, validate="one_to_one")
    long = counts.melt(
        id_vars=[*group_columns, "n_spots", "library_size"],
        value_vars=selected_genes,
        var_name="gene",
        value_name="raw_count",
    )
    long["log2_cpm"] = np.log2(
        (long["raw_count"] + 0.5) / (long["library_size"] + 1.0) * 1e6
    )
    return counts, long


def paired_standardized_effect(differences: np.ndarray) -> float:
    standard_deviation = float(np.std(differences, ddof=1))
    if not np.isfinite(standard_deviation) or standard_deviation == 0:
        return math.nan
    return float(np.mean(differences) / standard_deviation)


def exact_sign_flip(differences: np.ndarray) -> tuple[float, float, int]:
    observed = float(np.mean(differences))
    null = np.asarray([
        float(np.mean(differences * np.asarray(signs)))
        for signs in itertools.product((-1.0, 1.0), repeat=len(differences))
    ])
    tolerance = 1e-12
    return (
        float(np.mean(null >= observed - tolerance)),
        float(np.mean(np.abs(null) >= abs(observed) - tolerance)),
        len(null),
    )


def bootstrap_paired_mean(differences: np.ndarray, seed: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    indices = rng.integers(
        0, len(differences), size=(BOOTSTRAP_REPLICATES, len(differences))
    )
    means = differences[indices].mean(axis=1)
    lower, upper = np.quantile(means, [0.025, 0.975])
    return float(lower), float(upper)


def summarize_differences(
    values: pd.DataFrame,
    comparison: str,
    program: str,
    seed: int,
) -> dict[str, object]:
    differences = values["difference"].to_numpy(dtype=float)
    p_one_sided, p_two_sided, assignments = exact_sign_flip(differences)
    ci_lower, ci_upper = bootstrap_paired_mean(differences, seed)
    return {
        "comparison": comparison,
        "program": program,
        "n_pairs": len(differences),
        "positive_n": int((differences > 0).sum()),
        "negative_n": int((differences < 0).sum()),
        "mean_difference": float(differences.mean()),
        "median_difference": float(np.median(differences)),
        "bootstrap_ci_lower": ci_lower,
        "bootstrap_ci_upper": ci_upper,
        "paired_standardized_effect": paired_standardized_effect(differences),
        "exact_signflip_p_one_sided": p_one_sided,
        "exact_signflip_p_two_sided": p_two_sided,
        "n_exact_sign_assignments": assignments,
    }


def score_comparison(
    pseudobulk_long: pd.DataFrame,
    programs: dict[str, list[str]],
    region_a: str,
    region_b: str,
    comparison: str,
    unit_column: str,
    patient_column: str,
) -> pd.DataFrame:
    subset = pseudobulk_long[pseudobulk_long["region"].isin([region_a, region_b])].copy()
    metadata_columns = [unit_column, patient_column, "region"]
    if unit_column == patient_column:
        metadata_columns = [unit_column, "region"]
    wide = subset.pivot(
        index=metadata_columns,
        columns="gene",
        values="log2_cpm",
    ).reset_index()
    rows: list[dict[str, object]] = []
    for program_name, defined_genes in programs.items():
        genes = [gene for gene in defined_genes if gene in wide.columns and wide[gene].notna().all()]
        matrix = wide[genes].to_numpy(dtype=float)
        means = matrix.mean(axis=0)
        standard_deviations = matrix.std(axis=0, ddof=0)
        variable = standard_deviations > 0
        score = ((matrix[:, variable] - means[variable]) / standard_deviations[variable]).mean(axis=1)
        scored = wide[metadata_columns].copy()
        scored["score"] = score
        paired = scored.pivot(
            index=[column for column in metadata_columns if column != "region"],
            columns="region",
            values="score",
        ).reset_index()
        if region_a not in paired.columns or region_b not in paired.columns:
            continue
        paired = paired.dropna(subset=[region_a, region_b])
        for row in paired.itertuples(index=False):
            row_values = row._asdict()
            output = {
                "comparison": comparison,
                "program": program_name,
                "region_a": region_a,
                "region_b": region_b,
                "n_program_genes": int(variable.sum()),
                "unit": row_values[unit_column],
                "patient": row_values[patient_column],
                "region_a_score": row_values[region_a],
                "region_b_score": row_values[region_b],
                "difference": row_values[region_a] - row_values[region_b],
            }
            rows.append(output)

    if region_a == "macro_tumour" and region_b == "micro_tumour":
        rec_genes = [gene for gene in programs["rec_equal_gene"] if gene in subset["gene"].unique()]
        pooled = (
            subset[subset["gene"].isin(rec_genes)]
            .groupby(metadata_columns, as_index=False, observed=True)
            .agg(program_raw_count=("raw_count", "sum"), library_size=("library_size", "first"))
        )
        pooled["score"] = np.log2(
            (pooled["program_raw_count"] + 0.5) / (pooled["library_size"] + 1.0) * 1e6
        )
        paired = pooled.pivot(
            index=[column for column in metadata_columns if column != "region"],
            columns="region",
            values="score",
        ).reset_index().dropna(subset=[region_a, region_b])
        for row in paired.itertuples(index=False):
            row_values = row._asdict()
            rows.append({
                "comparison": comparison,
                "program": "rec_pooled_log2_cpm",
                "region_a": region_a,
                "region_b": region_b,
                "n_program_genes": len(rec_genes),
                "unit": row_values[unit_column],
                "patient": row_values[patient_column],
                "region_a_score": row_values[region_a],
                "region_b_score": row_values[region_b],
                "difference": row_values[region_a] - row_values[region_b],
            })
    return pd.DataFrame(rows)


def patient_region_counts(spots: pd.DataFrame) -> pd.DataFrame:
    return (
        spots.groupby(["patient", "region"], as_index=False, observed=True)
        .size()
        .rename(columns={"size": "n_spots"})
    )


def balanced_spot_resampling(
    spots: pd.DataFrame,
    programs: dict[str, list[str]],
    selected_genes: list[str],
) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    patients = sorted(spots["patient"].unique())
    program_subset = {
        name: genes
        for name, genes in programs.items()
        if name in {
            "rec_equal_gene",
            "hgp_cross_modal_top15",
            "experimental_ap1_target_overlap",
            "source_micrometastasis_six_gene",
        }
    }
    pools = {
        (patient, region): spots.index[
            (spots["patient"] == patient) & (spots["region"] == region)
        ].to_numpy()
        for patient in patients
        for region in ("micro_tumour", "macro_tumour")
    }
    rows: list[dict[str, object]] = []
    for iteration in range(1, BALANCED_ITERATIONS + 1):
        selected_rows: list[pd.DataFrame] = []
        for patient in patients:
            micro_indices = pools[(patient, "micro_tumour")]
            macro_indices = pools[(patient, "macro_tumour")]
            sampled_macro = rng.choice(macro_indices, size=len(micro_indices), replace=False)
            selected_rows.append(spots.loc[micro_indices])
            selected_rows.append(spots.loc[sampled_macro])
        iteration_spots = pd.concat(selected_rows, ignore_index=True)
        _, iteration_long = pseudobulk_spots(
            iteration_spots,
            ["patient", "region"],
            selected_genes,
        )
        values = score_comparison(
            iteration_long,
            program_subset,
            "macro_tumour",
            "micro_tumour",
            "macro_minus_micro_tumour",
            "patient",
            "patient",
        )
        for program, frame in values.groupby("program", sort=False):
            differences = frame["difference"].to_numpy(dtype=float)
            rows.append({
                "iteration": iteration,
                "program": program,
                "mean_paired_difference": float(differences.mean()),
                "median_paired_difference": float(np.median(differences)),
                "positive_patient_n": int((differences > 0).sum()),
                "positive": float(differences.mean()) > 0,
            })
    return pd.DataFrame(rows)


def gene_level_paired_effects(
    patient_long: pd.DataFrame,
    programs: dict[str, list[str]],
) -> pd.DataFrame:
    tumour = patient_long[patient_long["region"].isin(["micro_tumour", "macro_tumour"])]
    wide = tumour.pivot(
        index=["patient", "gene"], columns="region", values="log2_cpm"
    ).reset_index().dropna(subset=["micro_tumour", "macro_tumour"])
    wide["macro_minus_micro"] = wide["macro_tumour"] - wide["micro_tumour"]
    rows: list[dict[str, object]] = []
    for gene, frame in wide.groupby("gene", sort=False):
        differences = frame["macro_minus_micro"].to_numpy(dtype=float)
        p_one_sided, p_two_sided, assignments = exact_sign_flip(differences)
        rows.append({
            "gene": gene,
            "n_patients": len(differences),
            "mean_macro_minus_micro": float(differences.mean()),
            "median_macro_minus_micro": float(np.median(differences)),
            "positive_patient_n": int((differences > 0).sum()),
            "positive_patient_fraction": float((differences > 0).mean()),
            "paired_standardized_effect": paired_standardized_effect(differences),
            "exact_signflip_p_one_sided": p_one_sided,
            "exact_signflip_p_two_sided": p_two_sided,
            "n_exact_sign_assignments": assignments,
            **{
                f"member_{program_name}": gene in gene_set
                for program_name, gene_set in programs.items()
            },
        })
    result = pd.DataFrame(rows)
    order = np.argsort(result["exact_signflip_p_two_sided"].to_numpy())
    raw = result["exact_signflip_p_two_sided"].to_numpy()[order]
    adjusted_sorted = np.minimum.accumulate(
        (raw * len(raw) / np.arange(1, len(raw) + 1))[::-1]
    )[::-1]
    adjusted = np.empty(len(raw))
    adjusted[order] = np.minimum(adjusted_sorted, 1.0)
    result["BH_FDR_all_selected_genes"] = adjusted
    return result.sort_values(
        ["member_rec_equal_gene", "mean_macro_minus_micro"], ascending=False
    )


def main() -> None:
    defined_programs, selected_genes = load_programs()
    manifest = pd.read_csv(MANIFEST, sep="\t")
    annotation = pd.read_csv(ANNOTATION, sep="\t", compression="gzip")
    annotation = annotation[annotation["Layer3"].isin(REGION_MAP)].copy()
    annotation["region"] = annotation["Layer3"].map(REGION_MAP)

    spots_list: list[pd.DataFrame] = []
    qc_rows: list[dict[str, object]] = []
    present_by_sample: dict[str, set[str]] = {}
    for row in manifest.itertuples(index=False):
        spots, qc, present = read_sample(
            str(row.sample), str(row.patient), annotation, selected_genes
        )
        spots_list.append(spots)
        qc_rows.append(qc)
        present_by_sample[str(row.sample)] = present
        print(f"loaded {row.sample}", flush=True)
    spots = pd.concat(spots_list, ignore_index=True)
    qc = pd.DataFrame(qc_rows)
    common_genes = set(selected_genes)
    for present in present_by_sample.values():
        common_genes &= present
    missing = [gene for gene in selected_genes if gene not in common_genes]
    programs = {
        program_name: [gene for gene in genes if gene in common_genes]
        for program_name, genes in defined_programs.items()
    }
    if len(programs["rec_equal_gene"]) < 40:
        raise ValueError(
            f"Only {len(programs['rec_equal_gene'])} of 50 REC genes are common across samples"
        )
    selected_genes = [gene for gene in selected_genes if gene in common_genes]
    spots = spots[[
        "sample", "patient", "spot_barcode", "region", "library_size", *selected_genes
    ]].copy()

    patient_counts, patient_long = pseudobulk_spots(
        spots,
        ["patient", "region"],
        selected_genes,
    )
    sample_counts, sample_long = pseudobulk_spots(
        spots,
        ["sample", "patient", "region"],
        selected_genes,
    )

    comparisons = [
        ("macro_tumour", "micro_tumour", "macro_minus_micro_tumour"),
        ("micro_tumour", "micro_stroma", "micro_tumour_minus_stroma"),
        ("macro_tumour", "macro_stroma", "macro_tumour_minus_stroma"),
    ]
    score_tables: list[pd.DataFrame] = []
    for region_a, region_b, comparison in comparisons:
        score_tables.append(
            score_comparison(
                patient_long,
                programs,
                region_a,
                region_b,
                comparison,
                "patient",
                "patient",
            )
        )
    patient_scores = pd.concat(score_tables, ignore_index=True)
    contrast_rows: list[dict[str, object]] = []
    seed_counter = 0
    for (comparison, program), frame in patient_scores.groupby(
        ["comparison", "program"], sort=True
    ):
        contrast_rows.append(
            summarize_differences(frame, comparison, program, SEED + seed_counter)
        )
        seed_counter += 1
    paired_contrasts = pd.DataFrame(contrast_rows)

    region_spot_counts = patient_region_counts(spots)
    eligible_micro = region_spot_counts[
        (region_spot_counts["region"] == "micro_tumour")
        & (region_spot_counts["n_spots"] >= MIN_MICRO_SPOTS_SENSITIVITY)
    ]["patient"]
    minimum_spot_scores = patient_scores[
        (patient_scores["comparison"] == "macro_minus_micro_tumour")
        & (patient_scores["patient"].isin(eligible_micro))
    ].copy()
    minimum_spot_rows: list[dict[str, object]] = []
    for program, frame in minimum_spot_scores.groupby("program", sort=True):
        minimum_spot_rows.append(
            summarize_differences(
                frame,
                "macro_minus_micro_tumour_min20_micro_spots",
                program,
                SEED + seed_counter,
            )
        )
        seed_counter += 1
    minimum_spot_sensitivity = pd.DataFrame(minimum_spot_rows)

    paired_samples = (
        sample_counts.groupby(["sample", "patient"])["region"]
        .agg(lambda values: set(values))
        .reset_index()
    )
    paired_samples = paired_samples[
        paired_samples["region"].apply(
            lambda regions: {"micro_tumour", "macro_tumour"}.issubset(regions)
        )
    ]
    section_long = sample_long[
        sample_long["sample"].isin(paired_samples["sample"])
    ].copy()
    section_scores = score_comparison(
        section_long,
        programs,
        "macro_tumour",
        "micro_tumour",
        "within_section_macro_minus_micro",
        "sample",
        "patient",
    )
    within_section_patient = (
        section_scores.groupby(["program", "patient"], as_index=False)
        .agg(
            n_paired_sections=("unit", "nunique"),
            difference=("difference", "mean"),
        )
    )
    within_section_rows: list[dict[str, object]] = []
    for program, frame in within_section_patient.groupby("program", sort=True):
        within_section_rows.append(
            summarize_differences(
                frame,
                "within_section_macro_minus_micro_patient_mean",
                program,
                SEED + seed_counter,
            )
        )
        seed_counter += 1
    within_section_summary = pd.DataFrame(within_section_rows)

    balanced_iterations = balanced_spot_resampling(
        spots,
        programs,
        selected_genes,
    )
    balanced_summary = (
        balanced_iterations.groupby("program", as_index=False)
        .agg(
            iterations=("iteration", "nunique"),
            median_effect=("mean_paired_difference", "median"),
            mean_effect=("mean_paired_difference", "mean"),
            q025=("mean_paired_difference", lambda values: float(values.quantile(0.025))),
            q975=("mean_paired_difference", lambda values: float(values.quantile(0.975))),
            positive_fraction=("positive", "mean"),
            median_positive_patient_n=("positive_patient_n", "median"),
        )
    )
    gene_effects = gene_level_paired_effects(patient_long, programs)

    primary = paired_contrasts[
        (paired_contrasts["comparison"] == "macro_minus_micro_tumour")
        & (paired_contrasts["program"] == "rec_equal_gene")
    ].iloc[0]
    source_control = paired_contrasts[
        (paired_contrasts["comparison"] == "macro_minus_micro_tumour")
        & (paired_contrasts["program"] == "source_micrometastasis_six_gene")
    ].iloc[0]
    balanced_lookup = balanced_summary.set_index("program")
    source_control_pass = (
        source_control["mean_difference"] < 0 and int(source_control["positive_n"]) <= 3
    )
    if (
        primary["mean_difference"] > 0
        and int(primary["positive_n"]) >= 8
        and balanced_lookup.loc["rec_equal_gene", "positive_fraction"] >= 0.8
    ):
        label = "REC_PROGRAM_INCREASES_WITH_OUTGROWTH"
        interpretation = (
            "The fixed REC program is higher in liver macrometastasis tumour than micrometastasis tumour "
            "and remains positive after balancing spot counts, placing it with metastatic outgrowth."
        )
    elif (
        primary["mean_difference"] < 0
        and int(primary["positive_n"]) <= 3
        and balanced_lookup.loc["rec_equal_gene", "positive_fraction"] <= 0.2
    ):
        label = "REC_PROGRAM_HIGHER_IN_MICROMETASTATIC_PERSISTENCE"
        interpretation = (
            "The fixed REC program is higher in liver micrometastasis tumour and remains negative after "
            "spot balancing, placing it with early metastatic persistence."
        )
    else:
        label = "REC_PROGRAM_NOT_ORDERED_BY_LESION_SIZE"
        interpretation = (
            "The fixed REC program is not consistently ordered between liver micro- and macrometastasis "
            "tumour, so lesion size should not be used as its biological explanation."
        )
    decision = pd.DataFrame([{
        "decision": label,
        "primary_mean_macro_minus_micro": primary["mean_difference"],
        "primary_positive_patient_n": int(primary["positive_n"]),
        "primary_exact_signflip_p_two_sided": primary["exact_signflip_p_two_sided"],
        "primary_balanced_positive_fraction": balanced_lookup.loc[
            "rec_equal_gene", "positive_fraction"
        ],
        "source_six_gene_control_mean_macro_minus_micro": source_control["mean_difference"],
        "source_six_gene_control_positive_patient_n": int(source_control["positive_n"]),
        "source_six_gene_control_pass": source_control_pass,
        "interpretation": interpretation,
    }])

    coverage_rows = []
    for program_name, genes in defined_programs.items():
        coverage_rows.append({
            "program": program_name,
            "n_defined": len(genes),
            "n_common_present": sum(gene in common_genes for gene in genes),
            "matched_genes": ",".join(gene for gene in genes if gene in common_genes),
            "missing_genes": ",".join(gene for gene in genes if gene not in common_genes),
        })
    coverage = pd.DataFrame(coverage_rows)
    provenance_rows = []
    for row in manifest.itertuples(index=False):
        for filename in ["matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"]:
            path = EXTRACTED / row.sample / "filtered_feature_bc_matrix" / filename
            provenance_rows.append({"sample": row.sample, "patient": row.patient,
                "gsm": row.gsm, "source_url": row.url,
                "matrix_file": str(path.relative_to(ROOT)),
                "actual_size_bytes": path.stat().st_size, "sha256": sha256sum(path)})
    provenance = pd.DataFrame(provenance_rows)

    OUTPUT.mkdir(parents=True, exist_ok=True)
    outputs = {
        "program_coverage.tsv": coverage,
        "source_file_provenance.tsv": provenance,
        "sample_qc.tsv": qc,
        "patient_region_spot_counts.tsv": region_spot_counts,
        "patient_region_gene_pseudobulk.tsv.gz": patient_long,
        "sample_region_gene_pseudobulk.tsv.gz": sample_long,
        "patient_program_scores.tsv": patient_scores,
        "paired_contrasts.tsv": paired_contrasts,
        "minimum_micro_spot_sensitivity.tsv": minimum_spot_sensitivity,
        "within_section_patient_differences.tsv": within_section_patient,
        "within_section_summary.tsv": within_section_summary,
        "balanced_iteration_effects.tsv.gz": balanced_iterations,
        "balanced_summary.tsv": balanced_summary,
        "gene_level_paired_effects.tsv": gene_effects,
        "decision.tsv": decision,
        "run_info.tsv": pd.DataFrame([{
            "analysis_date": "2026-08-23",
            "python_version": platform.python_version(),
            "numpy_version": np.__version__,
            "pandas_version": pd.__version__,
            "random_seed": SEED,
            "bootstrap_replicates": BOOTSTRAP_REPLICATES,
            "balanced_iterations": BALANCED_ITERATIONS,
            "primary_unit": "patient",
            "primary_contrast": "macro_tumour_minus_micro_tumour",
        }]),
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

    primary_programs = paired_contrasts[
        (paired_contrasts["comparison"] == "macro_minus_micro_tumour")
    ].set_index("program")
    ap1 = primary_programs.loc["experimental_ap1_target_overlap"]
    top15 = primary_programs.loc["hgp_cross_modal_top15"]
    pooled = primary_programs.loc["rec_pooled_log2_cpm"]
    lines = [
        "# REC program across paired liver micro- and macrometastases",
        "",
        f"Decision: **{label}**.",
        "",
        "## Primary paired result",
        "",
        f"Across 11 patients, the macro-minus-micro equal-gene REC score averaged "
        f"{primary['mean_difference']:+.3f} (median {primary['median_difference']:+.3f}; "
        f"{int(primary['positive_n'])}/11 patients positive; bootstrap 95% CI "
        f"{primary['bootstrap_ci_lower']:+.3f} to {primary['bootstrap_ci_upper']:+.3f}; "
        f"paired standardized effect {primary['paired_standardized_effect']:+.3f}; exact sign-flip "
        f"P={primary['exact_signflip_p_two_sided']:.4f}).",
        "",
        f"The pooled 50-gene CPM difference was {pooled['mean_difference']:+.3f}, and the compact "
        f"HGP cross-modal core difference was {top15['mean_difference']:+.3f}.",
        "",
        "## Mechanistic and source controls",
        "",
        f"The experimental AP-1 target subset differed by {ap1['mean_difference']:+.3f} "
        f"({int(ap1['positive_n'])}/11 positive). The source micrometastasis six-gene control differed "
        f"by {source_control['mean_difference']:+.3f} ({int(source_control['positive_n'])}/11 positive); "
        f"control pass={source_control_pass}.",
        "",
        "## Spot-balance sensitivity",
        "",
        f"After matching macro spot counts to micro spot counts, the mean REC effect was positive in "
        f"{balanced_lookup.loc['rec_equal_gene', 'positive_fraction']:.1%} of "
        f"{BALANCED_ITERATIONS} iterations.",
        "",
        "## Interpretation",
        "",
        interpretation,
        " GSE294385 has no HGP labels, so this result positions the program along lesion outgrowth but "
        "does not assign a metastatic growth pattern.",
        "",
        "## Reproducibility",
        "",
        f"Specification: `{SPEC.relative_to(ROOT)}`. Extracted matrix inputs were size-checked and SHA-256 hashed by the public runner.",
    ]
    (OUTPUT / "analysis_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    manifest_lines = [
        "# Analysis output manifest",
        "",
        "Generated by `scripts/analyze_gse294385_rec_program_extension.py`.",
        "",
        "| File | Purpose |",
        "|---|---|",
        "| `analysis_summary.md` | Human-readable paired extension result |",
        "| `program_coverage.tsv` | Fixed program coverage |",
        "| `source_file_provenance.tsv` | GEO URLs, sizes and SHA-256 hashes |",
        "| `sample_qc.tsv` | Matrix/annotation matching by sample |",
        "| `patient_region_spot_counts.tsv` | Patient-level region counts |",
        "| `patient_region_gene_pseudobulk.tsv.gz` | Auditable patient-region gene values |",
        "| `sample_region_gene_pseudobulk.tsv.gz` | Section-level sensitivity input |",
        "| `patient_program_scores.tsv` | Paired scores for tumour and stroma comparisons |",
        "| `paired_contrasts.tsv` | Primary and compartment paired effects |",
        "| `minimum_micro_spot_sensitivity.tsv` | At least 20 micro spots sensitivity |",
        "| `within_section_patient_differences.tsv` | Within-section paired differences |",
        "| `within_section_summary.tsv` | Within-section paired summaries |",
        "| `balanced_iteration_effects.tsv.gz` | Macro spot-count downsampling effects |",
        "| `balanced_summary.tsv` | Spot-balance sensitivity summary |",
        "| `gene_level_paired_effects.tsv` | Gene-level macro-minus-micro directions |",
        "| `decision.tsv` | Prespecified interpretation |",
        "| `run_info.tsv` | Reproducibility settings |",
    ]
    (OUTPUT / "_analysis_outputs.md").write_text(
        "\n".join(manifest_lines) + "\n", encoding="utf-8"
    )
    print(decision.to_string(index=False))
    print(primary_programs.to_string())


if __name__ == "__main__":
    main()
