#!/usr/bin/env python3
"""
Analysis: Preranked GSEA for paired Ogden REC-versus-state contrasts.
Date: 2026-09-01
Random seed: 42
Python and package versions are written to the provenance JSON.

The predeclared registry and Hallmark collection are tested in all eligible
contrasts. Complete Reactome and GO Biological Process libraries are used as
exploratory discovery layers only for the three principal stress-state
comparisons (Hypoxia, UPR and iREC).
"""

from __future__ import annotations

import hashlib
import json
import platform
from importlib.metadata import version
from pathlib import Path

import gseapy as gp
import numpy as np
import pandas as pd


np.random.seed(42)

ROOT = Path(__file__).resolve().parents[1]
OUT = (
    ROOT
    / "analysis_results"
    / "deep_biology_upgrade_2026-08-31"
    / "phase2_mechanistic_specificity"
)
RANKS = OUT / "ogden_direct_state_limma_all.tsv.gz"
REGISTRY = OUT / "frozen_gene_set_registry.tsv"
FROZEN_LONG = OUT / "frozen_gene_sets_long.tsv"
LIBRARIES = {
    "targeted_predeclared": OUT / "frozen_gene_sets.gmt",
    "full_hallmark_discovery": OUT / "msigdb_h.all_2025.1.Hs.gmt",
    "full_reactome_discovery": OUT / "msigdb_c2.cp.reactome_2025.1.Hs.gmt",
    "full_gobp_discovery": OUT / "msigdb_c5.go.bp_2025.1.Hs.gmt",
}

PRINCIPAL_COMPARATORS = {"Hypoxia", "UPR", "iREC"}
MIN_SIZE = 5
MAX_SIZE = 500
TARGETED_PERMUTATIONS = 2000
DISCOVERY_PERMUTATIONS = 1000
THREADS = 4
SEED = 42


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_gmt(path: Path) -> dict[str, list[str]]:
    collection: dict[str, list[str]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(f"Malformed GMT line {line_number} in {path}")
            term = fields[0]
            genes = list(dict.fromkeys(gene.strip().upper() for gene in fields[2:] if gene.strip()))
            if term in collection:
                raise ValueError(f"Duplicate term {term} in {path}")
            collection[term] = genes
    return collection


def normalize_result(result: pd.DataFrame) -> pd.DataFrame:
    rename = {
        "Term": "term",
        "ES": "es",
        "NES": "nes",
        "NOM p-val": "nominal_p",
        "FDR q-val": "fdr_q",
        "FWER p-val": "fwer_p",
        "Tag %": "tag_fraction",
        "Gene %": "rank_fraction",
        "Lead_genes": "leading_edge_genes",
    }
    normalized = result.rename(columns=rename).copy()
    expected = set(rename.values())
    missing = sorted(expected.difference(normalized.columns))
    if missing:
        raise ValueError(f"GSEA result missing columns: {missing}")
    keep = [
        "term",
        "es",
        "nes",
        "nominal_p",
        "fdr_q",
        "fwer_p",
        "tag_fraction",
        "rank_fraction",
        "leading_edge_genes",
    ]
    normalized = normalized.loc[:, keep]
    for column in ["es", "nes", "nominal_p", "fdr_q", "fwer_p"]:
        normalized[column] = pd.to_numeric(normalized[column], errors="coerce")
    return normalized


def run_one(
    ranks: pd.DataFrame,
    gene_sets: dict[str, list[str]],
    permutations: int,
) -> pd.DataFrame:
    result = gp.prerank(
        rnk=ranks.loc[:, ["gene", "moderated_t"]],
        gene_sets=gene_sets,
        min_size=MIN_SIZE,
        max_size=MAX_SIZE,
        permutation_num=permutations,
        weight=1.0,
        ascending=False,
        threads=THREADS,
        seed=SEED,
        outdir=None,
        verbose=False,
    )
    return normalize_result(result.res2d)


def main() -> None:
    required = [RANKS, REGISTRY, FROZEN_LONG, *LIBRARIES.values()]
    missing = [path for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing GSEA inputs: {missing}")

    ranks_all = pd.read_csv(RANKS, sep="\t")
    required_columns = {"comparator_state", "gene", "moderated_t"}
    if not required_columns.issubset(ranks_all.columns):
        raise ValueError("Rank table does not contain the expected columns")
    if ranks_all.duplicated(["comparator_state", "gene"]).any():
        raise ValueError("Rank table contains duplicate comparator-gene rows")

    registry = pd.read_csv(REGISTRY, sep="\t")
    frozen_long = pd.read_csv(FROZEN_LONG, sep="\t")
    excluded_targeted = set(
        registry.loc[registry["role"].eq("overlap_audit_only"), "set_id"]
    )
    frozen_sets = {
        term: group.sort_values("gene_rank")["gene"].astype(str).str.upper().tolist()
        for term, group in frozen_long.groupby("set_id", sort=False)
        if term not in excluded_targeted
    }
    library_sets = {
        "targeted_predeclared": frozen_sets,
        "full_hallmark_discovery": parse_gmt(LIBRARIES["full_hallmark_discovery"]),
        "full_reactome_discovery": parse_gmt(LIBRARIES["full_reactome_discovery"]),
        "full_gobp_discovery": parse_gmt(LIBRARIES["full_gobp_discovery"]),
    }

    results: list[pd.DataFrame] = []
    comparators = list(dict.fromkeys(ranks_all["comparator_state"].astype(str)))
    for comparator in comparators:
        ranks = (
            ranks_all.loc[
                ranks_all["comparator_state"].eq(comparator),
                ["gene", "moderated_t"],
            ]
            .dropna()
            .sort_values("moderated_t", ascending=False)
        )
        if ranks["gene"].duplicated().any() or len(ranks) < 1000:
            raise ValueError(f"Invalid ranked list for comparator {comparator}")

        modes = ["targeted_predeclared", "full_hallmark_discovery"]
        if comparator in PRINCIPAL_COMPARATORS:
            modes.extend(["full_reactome_discovery", "full_gobp_discovery"])
        for mode in modes:
            permutations = (
                TARGETED_PERMUTATIONS
                if mode in {"targeted_predeclared", "full_hallmark_discovery"}
                else DISCOVERY_PERMUTATIONS
            )
            result = run_one(ranks, library_sets[mode], permutations)
            result.insert(0, "comparator_state", comparator)
            result.insert(1, "analysis_mode", mode)
            result.insert(2, "n_ranked_genes", len(ranks))
            result.insert(3, "permutations", permutations)
            results.append(result)
            print(
                f"Completed {comparator}: {mode} "
                f"({len(result)} tested gene sets)",
                flush=True,
            )

    combined = pd.concat(results, ignore_index=True)
    combined["direction"] = np.where(combined["nes"] >= 0, "REC_up", "REC_down")
    combined = combined.sort_values(
        ["comparator_state", "analysis_mode", "fdr_q", "nominal_p", "term"],
        na_position="last",
    )
    combined.to_csv(
        OUT / "ogden_direct_state_gsea_all.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    targeted = combined[combined["analysis_mode"].eq("targeted_predeclared")].copy()
    targeted = targeted.merge(
        registry.loc[:, ["set_id", "family", "role", "citation"]],
        left_on="term",
        right_on="set_id",
        how="left",
        validate="many_to_one",
    ).drop(columns="set_id")
    targeted.to_csv(
        OUT / "ogden_direct_state_gsea_targeted.tsv",
        sep="\t",
        index=False,
    )

    discovery_key = combined[
        ~combined["analysis_mode"].eq("targeted_predeclared")
        & combined["fdr_q"].le(0.05)
    ].copy()
    discovery_key.to_csv(
        OUT / "ogden_direct_state_gsea_discovery_fdr05.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    principal = combined[
        combined["comparator_state"].isin(PRINCIPAL_COMPARATORS)
        & combined["fdr_q"].le(0.05)
    ].copy()
    recurrence = (
        principal.groupby(["analysis_mode", "term", "direction"], as_index=False)
        .agg(
            n_principal_contrasts=("comparator_state", "nunique"),
            principal_comparators=(
                "comparator_state",
                lambda values: ";".join(sorted(set(values))),
            ),
            median_nes=("nes", "median"),
            maximum_fdr_q=("fdr_q", "max"),
        )
        .sort_values(
            ["n_principal_contrasts", "analysis_mode", "median_nes"],
            ascending=[False, True, False],
        )
    )
    recurrence.to_csv(
        OUT / "ogden_direct_state_gsea_principal_recurrence.tsv",
        sep="\t",
        index=False,
    )

    leading_rows: list[dict[str, object]] = []
    for row in principal.itertuples(index=False):
        genes = str(row.leading_edge_genes).split(";")
        for gene in genes:
            if gene and gene.lower() != "nan":
                leading_rows.append(
                    {
                        "comparator_state": row.comparator_state,
                        "analysis_mode": row.analysis_mode,
                        "term": row.term,
                        "direction": row.direction,
                        "gene": gene,
                    }
                )
    leading = pd.DataFrame(leading_rows)
    if not leading.empty:
        leading_recurrence = (
            leading.groupby(["direction", "gene"], as_index=False)
            .agg(
                n_principal_contrasts=("comparator_state", "nunique"),
                n_significant_terms=("term", "nunique"),
                principal_comparators=(
                    "comparator_state",
                    lambda values: ";".join(sorted(set(values))),
                ),
            )
            .sort_values(
                ["n_principal_contrasts", "n_significant_terms", "gene"],
                ascending=[False, False, True],
            )
        )
    else:
        leading_recurrence = pd.DataFrame(
            columns=[
                "direction",
                "gene",
                "n_principal_contrasts",
                "n_significant_terms",
                "principal_comparators",
            ]
        )
    leading_recurrence.to_csv(
        OUT / "ogden_direct_state_gsea_leading_edge_recurrence.tsv",
        sep="\t",
        index=False,
    )

    provenance = {
        "analysis_date": "2026-09-01",
        "python_version": platform.python_version(),
        "pandas_version": pd.__version__,
        "numpy_version": np.__version__,
        "gseapy_version": version("gseapy"),
        "seed": SEED,
        "threads": THREADS,
        "min_size": MIN_SIZE,
        "max_size": MAX_SIZE,
        "targeted_and_hallmark_permutations": TARGETED_PERMUTATIONS,
        "broad_discovery_permutations": DISCOVERY_PERMUTATIONS,
        "fdr_scope": "calculated by GSEA within each comparator-library run",
        "inputs": {str(path.relative_to(ROOT)): sha256(path) for path in required},
    }
    with (OUT / "ogden_direct_state_gsea_provenance.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(provenance, handle, ensure_ascii=False, indent=2)

    print(
        targeted.loc[
            targeted["term"].isin(
                [
                    "FROZEN_REC_TOP50",
                    "CANELLAS_CORE_HRC",
                    "OGDEN_EPIHR",
                    "HALLMARK_HYPOXIA",
                    "HALLMARK_E2F_TARGETS",
                    "HALLMARK_MYC_TARGETS_V1",
                    "HALLMARK_FATTY_ACID_METABOLISM",
                    "PENGWINKLER_PALMITATE_UP_TOP100",
                ]
            ),
            ["comparator_state", "term", "nes", "nominal_p", "fdr_q"],
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
