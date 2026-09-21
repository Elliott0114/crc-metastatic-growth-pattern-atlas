#!/usr/bin/env python3
"""Post hoc REC comparison with the published CDH17-positive C8 marker set."""

from __future__ import annotations

import hashlib
import json
import platform
from importlib.metadata import version
from pathlib import Path

import gseapy as gp
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


ROOT = Path(__file__).resolve().parents[1]
OUT = (
    ROOT
    / "analysis_results"
    / "deep_biology_upgrade_2026-08-31"
    / "phase2_mechanistic_specificity"
)
WORKBOOK = (
    ROOT
    / "data_sources"
    / "AlvarezVillanueva_2026_tight_junction"
    / "41467_2025_68169_MOESM5_ESM.xlsx"
)
C8_SOURCE = OUT / "alvarez_c8_top50_source.tsv"
RANKS = OUT / "ogden_direct_state_limma_all.tsv.gz"
RECURRENT = OUT / "ogden_rec_principal_consistent_genes.tsv"
FROZEN_LONG = OUT / "frozen_gene_sets_long.tsv"
REGISTRY = OUT / "frozen_gene_set_registry.tsv"
PRINCIPAL = ["Hypoxia", "UPR", "iREC"]
NAMED_GENES = ["CDH17", "CLDN2", "TJP1", "LGR5", "MKI67"]
SEED = 42
PERMUTATIONS = 2000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def gene_hash(genes: list[str]) -> str:
    payload = "\n".join(genes) + "\n"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def normalize_result(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.rename(
        columns={
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
    )


def main() -> None:
    required = [WORKBOOK, C8_SOURCE, RANKS, RECURRENT, FROZEN_LONG, REGISTRY]
    missing = [path for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing C8 bridge inputs: {missing}")

    c8 = pd.read_csv(C8_SOURCE, sep="\t")
    c8_genes = c8.sort_values("source_rank")["gene"].astype(str).str.upper().tolist()
    if len(c8_genes) != 50 or len(set(c8_genes)) != 50:
        raise ValueError("The extracted C8 source signature is not 50 unique genes")

    frozen = pd.read_csv(FROZEN_LONG, sep="\t")
    registry = pd.read_csv(REGISTRY, sep="\t")
    excluded = set(registry.loc[registry["role"].eq("overlap_audit_only"), "set_id"])
    gene_sets = {
        set_id: group.sort_values("gene_rank")["gene"].astype(str).str.upper().tolist()
        for set_id, group in frozen.groupby("set_id", sort=False)
        if set_id not in excluded
    }
    hallmark_hypoxia = set(gene_sets["HALLMARK_HYPOXIA"])
    c8_nonhypoxia = [gene for gene in c8_genes if gene not in hallmark_hypoxia]
    gene_sets["ALVAREZ_C8_TOP50_POSTHOC"] = c8_genes
    gene_sets["ALVAREZ_C8_TOP50_MINUS_HALLMARK_HYPOXIA_POSTINSPECTION"] = c8_nonhypoxia

    ranks_all = pd.read_csv(RANKS, sep="\t")
    if ranks_all.duplicated(["comparator_state", "gene"]).any():
        raise ValueError("Direct-state rank table contains duplicated comparator-gene rows")

    gsea_rows: list[pd.DataFrame] = []
    for comparator in PRINCIPAL:
        ranks = (
            ranks_all.loc[
                ranks_all["comparator_state"].eq(comparator),
                ["gene", "moderated_t"],
            ]
            .dropna()
            .sort_values("moderated_t", ascending=False)
        )
        result = gp.prerank(
            rnk=ranks,
            gene_sets=gene_sets,
            min_size=5,
            max_size=500,
            permutation_num=PERMUTATIONS,
            weight=1.0,
            ascending=False,
            threads=4,
            seed=SEED,
            outdir=None,
            verbose=False,
        )
        normalized = normalize_result(result.res2d)
        normalized = normalized.loc[
            normalized["term"].isin(
                {
                    "ALVAREZ_C8_TOP50_POSTHOC",
                    "ALVAREZ_C8_TOP50_MINUS_HALLMARK_HYPOXIA_POSTINSPECTION",
                }
            )
        ].copy()
        normalized.insert(0, "comparator_state", comparator)
        normalized.insert(1, "n_ranked_genes", len(ranks))
        normalized.insert(2, "permutations", PERMUTATIONS)
        gsea_rows.append(normalized)

    gsea = pd.concat(gsea_rows, ignore_index=True)
    gsea["direction"] = np.where(gsea["nes"].ge(0), "REC_up", "REC_down")
    gsea["source_signature_size"] = gsea["term"].map(
        {
            "ALVAREZ_C8_TOP50_POSTHOC": len(c8_genes),
            "ALVAREZ_C8_TOP50_MINUS_HALLMARK_HYPOXIA_POSTINSPECTION": len(c8_nonhypoxia),
        }
    )
    gsea.to_csv(OUT / "rec_c8_tight_junction_gsea.tsv", sep="\t", index=False)

    tested_sets = [
        set(ranks_all.loc[ranks_all["comparator_state"].eq(name), "gene"].astype(str).str.upper())
        for name in PRINCIPAL
    ]
    universe = set.intersection(*tested_sets)
    recurrent = pd.read_csv(RECURRENT, sep="\t")
    recurrent_up = set(
        recurrent.loc[recurrent["direction"].eq("REC_up_all_three"), "gene"]
        .astype(str)
        .str.upper()
    ) & universe
    c8_tested = set(c8_genes) & universe
    overlap = c8_tested & recurrent_up
    table = np.array(
        [
            [len(overlap), len(c8_tested - recurrent_up)],
            [len(recurrent_up - c8_tested), len(universe - c8_tested - recurrent_up)],
        ]
    )
    odds_ratio, fisher_p = fisher_exact(table, alternative="greater")
    overlap_summary = pd.DataFrame(
        [
            {
                "universe_size": len(universe),
                "c8_top50_tested": len(c8_tested),
                "recurrent_rec_up": len(recurrent_up),
                "observed_overlap": len(overlap),
                "expected_overlap": len(c8_tested) * len(recurrent_up) / len(universe),
                "odds_ratio": odds_ratio,
                "fisher_exact_p_greater": fisher_p,
                "overlap_genes": ";".join(sorted(overlap)),
            }
        ]
    )
    overlap_summary.to_csv(OUT / "rec_c8_recurrent_overlap.tsv", sep="\t", index=False)

    named = ranks_all.loc[
        ranks_all["comparator_state"].isin(PRINCIPAL)
        & ranks_all["gene"].astype(str).str.upper().isin(NAMED_GENES)
    ].copy()
    named["gene"] = named["gene"].astype(str).str.upper()
    named = named.sort_values(
        ["gene", "comparator_state"],
        key=lambda values: values.map(
            {**{gene: i for i, gene in enumerate(NAMED_GENES)}, **{name: i for i, name in enumerate(PRINCIPAL)}}
        ).fillna(len(NAMED_GENES) + len(PRINCIPAL)),
    )
    named.to_csv(OUT / "rec_c8_named_gene_evidence.tsv", sep="\t", index=False)

    gsea_evidence = pd.DataFrame(
        {
            "record_type": "gene_set_gsea",
            "feature_id": gsea["term"],
            "comparator_state": gsea["comparator_state"],
            "n_source_genes": gsea["source_signature_size"],
            "effect_measure": "NES_REC_minus_comparator",
            "effect": gsea["nes"],
            "nominal_p": gsea["nominal_p"],
            "fdr": gsea["fdr_q"],
            "secondary_measure": "tag_fraction",
            "secondary_value": gsea["tag_fraction"],
            "supporting_genes": gsea["leading_edge_genes"],
            "note": gsea["term"].map(
                lambda value: "source top-50 C8 markers"
                if value == "ALVAREZ_C8_TOP50_POSTHOC"
                else "post-inspection sensitivity after removing Hallmark Hypoxia genes"
            ),
        }
    )
    overlap_evidence = pd.DataFrame(
        [
            {
                "record_type": "recurrent_gene_overlap",
                "feature_id": "C8_TOP50_vs_RECURRENT_REC_UP",
                "comparator_state": "all_three_principal_contrasts",
                "n_source_genes": len(c8_tested),
                "effect_measure": "Fisher_odds_ratio",
                "effect": odds_ratio,
                "nominal_p": fisher_p,
                "fdr": np.nan,
                "secondary_measure": "observed_vs_expected_overlap",
                "secondary_value": f"{len(overlap)} vs {len(c8_tested) * len(recurrent_up) / len(universe):.3f}",
                "supporting_genes": ";".join(sorted(overlap)),
                "note": f"common tested-gene universe n={len(universe)}",
            }
        ]
    )
    named_evidence = pd.DataFrame(
        {
            "record_type": "named_gene",
            "feature_id": named["gene"],
            "comparator_state": named["comparator_state"],
            "n_source_genes": 1,
            "effect_measure": "log2FC_REC_minus_comparator",
            "effect": named["log2fc_rec_minus_comparator"],
            "nominal_p": named["p_value"],
            "fdr": named["fdr_bh"],
            "secondary_measure": "moderated_t",
            "secondary_value": named["moderated_t"],
            "supporting_genes": named["gene"],
            "note": "named feature from the functional C8 comparison",
        }
    )
    pd.concat(
        [gsea_evidence, overlap_evidence, named_evidence], ignore_index=True
    ).to_csv(
        OUT / "rec_c8_tight_junction_bridge_evidence.tsv", sep="\t", index=False
    )

    hypoxia_overlap = [gene for gene in c8_genes if gene in hallmark_hypoxia]
    c8_rows = gsea.loc[gsea["term"].eq("ALVAREZ_C8_TOP50_POSTHOC")].set_index(
        "comparator_state"
    )
    nonhyp_rows = gsea.loc[
        gsea["term"].eq(
            "ALVAREZ_C8_TOP50_MINUS_HALLMARK_HYPOXIA_POSTINSPECTION"
        )
    ].set_index("comparator_state")
    summary_lines = [
        "# REC–C8 tight-junction metastasis bridge summary",
        "",
        "This literature-triggered post hoc analysis compares REC with the published top-50 marker set of the CDH17-positive C8 population. The source study functionally linked C8/CLDN2 to clustered migration and liver colonization, but did not annotate HGP.",
        "",
        f"- {len(hypoxia_overlap)}/50 C8 markers overlap the frozen Hallmark Hypoxia set: {', '.join(hypoxia_overlap)}.",
    ]
    for comparator in PRINCIPAL:
        row = c8_rows.loc[comparator]
        dehyp = nonhyp_rows.loc[comparator]
        summary_lines.append(
            f"- REC versus {comparator}: full C8 NES {row.nes:.2f} "
            f"(nominal P {row.nominal_p:.3g}, FDR {row.fdr_q:.3g}); "
            f"after removing Hallmark-hypoxia genes, NES {dehyp.nes:.2f} "
            f"(nominal P {dehyp.nominal_p:.3g}, FDR {dehyp.fdr_q:.3g})."
        )
    overlap_row = overlap_summary.iloc[0]
    named_counts: dict[str, tuple[int, int]] = {}
    for gene in NAMED_GENES:
        rows = named.loc[named["gene"].eq(gene)]
        expected_positive = gene not in {"LGR5", "MKI67"}
        direction_count = int(
            rows["log2fc_rec_minus_comparator"].gt(0).sum()
            if expected_positive
            else rows["log2fc_rec_minus_comparator"].lt(0).sum()
        )
        fdr_count = int(
            (
                rows["fdr_bh"].lt(0.05)
                & (
                    rows["log2fc_rec_minus_comparator"].gt(0)
                    if expected_positive
                    else rows["log2fc_rec_minus_comparator"].lt(0)
                )
            ).sum()
        )
        named_counts[gene] = (direction_count, fdr_count)
    summary_lines.extend(
        [
            f"- C8 top-50 markers overlap the 138 recurrent REC-up genes at {int(overlap_row.observed_overlap)} genes versus {overlap_row.expected_overlap:.2f} expected (odds ratio {overlap_row.odds_ratio:.2f}; one-sided Fisher P {overlap_row.fisher_exact_p_greater:.3g}): {overlap_row.overlap_genes or 'none'}.",
            f"- Named features were only partly concordant: CDH17 and TJP1 were REC-up in {named_counts['CDH17'][0]}/3 and {named_counts['TJP1'][0]}/3 contrasts (FDR<0.05 in {named_counts['CDH17'][1]}/3 and {named_counts['TJP1'][1]}/3), whereas CLDN2 was REC-up in {named_counts['CLDN2'][0]}/3. LGR5 and MKI67 were REC-down in {named_counts['LGR5'][0]}/3 and {named_counts['MKI67'][0]}/3 contrasts.",
            "",
            "The full and narrowly de-hypoxia C8 sets were strongly enriched, but their leading edges remained rich in hypoxia/stress genes outside the Hallmark set and CLDN2 itself was not increased. The result therefore supports partial convergence with a low-cycling, metastasis-competent epithelial state, not identity with C8 or evidence that IKK-alpha/CLDN2 drives REC or replacement growth.",
            "",
        ]
    )
    (OUT / "rec_c8_tight_junction_bridge_summary.md").write_text(
        "\n".join(summary_lines), encoding="utf-8"
    )

    provenance = {
        "analysis_date": "2026-09-01",
        "analysis_status": "literature-triggered post hoc; de-hypoxia sensitivity added after source-signature inspection",
        "source_doi": "10.1038/s41467-025-68169-3",
        "source_workbook_url": "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-025-68169-3/MediaObjects/41467_2025_68169_MOESM5_ESM.xlsx",
        "source_sheet": "Cluster_8",
        "source_selection": "first 50 source-ranked markers, matching the publication's top-50 C8 analysis",
        "seed": SEED,
        "permutations": PERMUTATIONS,
        "python_version": platform.python_version(),
        "gseapy_version": version("gseapy"),
        "pandas_version": pd.__version__,
        "numpy_version": np.__version__,
        "input_sha256": {str(path.relative_to(ROOT)): sha256(path) for path in required},
        "c8_top50_gene_sha256": gene_hash(c8_genes),
        "c8_minus_hallmark_hypoxia_gene_sha256": gene_hash(c8_nonhypoxia),
        "hallmark_hypoxia_overlap_count": len(hypoxia_overlap),
    }
    (OUT / "rec_c8_tight_junction_bridge_provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print("\n".join(summary_lines))


if __name__ == "__main__":
    main()
