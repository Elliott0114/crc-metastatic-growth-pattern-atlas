"""Member-level evidence tables and figures for the focused evidence bridge."""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from rec_junction_context_common import ROOT, bh, definitions, digest, mean_summary
from analyze_rec_evidence_bridge_spatial import OUT, group_summary, write

PREVIOUS = ROOT / "analysis_results/rec_junction_context_2026-09-16"
COLORS = {"rHGP": "#B35806", "dHGP": "#2166AC", "Claudin": "#0072B2",
          "Polarity": "#D55E00", "HRC": "#009E73", "Junction": "#7755A5"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                     "pdf.fonttype": 42, "ps.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False})


def save(fig, name):
    (OUT / "figures").mkdir(exist_ok=True)
    fig.savefig(OUT / f"figures/{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"figures/{name}.png", dpi=320, bbox_inches="tight")
    plt.close(fig)


def main():
    sets = definitions()
    cor = pd.read_csv(OUT / "spatial/patient_correlations.tsv", sep="\t")
    direct = []
    for restriction in ["all_tumour", "epithelial_high"]:
        d = cor[(cor.restriction == restriction) & (cor.model == "adjusted")
                & cor.component.isin(["Claudin", "Polarity"])]
        wide = d.pivot(index=["patient", "hgp"], columns="component", values="rho")
        delta = np.arctanh(wide.Claudin) - np.arctanh(wide.Polarity)
        write(delta.rename("Claudin_minus_Polarity_Fisher_z").reset_index(),
              f"integration/component_difference_{restriction}.tsv")
        direct.append(dict(restriction=restriction, role="post_primary_direct_comparison",
                           **mean_summary(delta, "spatial_direct_component_difference" + restriction)))
    write(pd.DataFrame(direct), "integration/component_difference_summary.tsv")

    # Explicit overlap: these are linked assays, not additional independent patients.
    epithelial = pd.read_csv(PREVIOUS / "sensitivities/hgp_epithelial_pseudobulk_metadata.tsv", sep="\t")
    spatial_meta = pd.read_csv(OUT / "spatial/sample_qc.tsv", sep="\t")
    overlap = spatial_meta[["patient", "hgp"]].merge(epithelial[["patient", "hgp", "n_cells"]],
                                                       how="outer", on=["patient", "hgp"], indicator=True)
    write(overlap, "integration/hgp_patient_overlap.tsv")

    # Align genes and normalization before comparing model tests with patient summaries.
    count_matrix = pd.read_csv(OUT / "spatial/region_counts.tsv.gz", sep="\t", index_col=0)
    normalizations = pd.read_csv(OUT / "counts/library_normalization.tsv", sep="\t")
    eligibility = pd.read_csv(OUT / "counts/gene_eligibility.tsv.gz", sep="\t")
    spatial_defs = pd.read_csv(OUT / "spatial/definitions.tsv", sep="\t")
    matched_effects, matched_scores, coverage_effects = [], [], []
    for dataset in ["Spatial_rHGP_minus_dHGP_tumour", "Spatial_rHGP_minus_dHGP_near"]:
        retained = set(eligibility[(eligibility.dataset == dataset) & eligibility.retained].gene)
        for norm in ["TMM", "library_only"]:
            m = normalizations[(normalizations.dataset == dataset) & (normalizations.normalization == norm)].set_index("id")
            c = count_matrix.loc[count_matrix.index.isin(retained), m.index]
            lib = m.effective_library
            prior = .5 * lib / lib.mean()
            expr = np.log2(c.add(prior, axis=1).div(lib + 2 * prior, axis=1) * 1e6)
            complete_expr = np.log2(count_matrix[m.index].add(prior, axis=1).div(lib + 2 * prior, axis=1) * 1e6)
            for membership, gg in [("all_assayed", sets["Claudin"]), ("count_eligible", sorted(set(sets["Claudin"]) & retained))]:
                values = complete_expr.loc[gg].mean(axis=0)
                coverage_effects.append(dict(dataset=dataset, normalization=norm, membership=membership,
                    n_genes=len(gg), **group_summary(values, m.hgp, "coverage" + dataset + norm + membership)))
            z = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1, ddof=0).replace(0, np.nan), axis=0)
            for comp, d in spatial_defs.groupby("component"):
                gg = sorted(set(d.gene) & retained)
                for construction, matrix in [("mean_logCPM", expr), ("mean_gene_z", z)]:
                    values = matrix.loc[gg].mean(axis=0)
                    result = group_summary(values, m.hgp, dataset + norm + comp + construction)
                    matched_effects.append(dict(dataset=dataset, normalization=norm, component=comp,
                        score=construction, n_genes=len(gg), **result))
                    matched_scores.extend(m.assign(value=values, component=comp, score=construction,
                                                   n_genes=len(gg)).reset_index().to_dict("records"))
    matched_effects = pd.DataFrame(matched_effects)
    matched_effects["q"] = matched_effects.groupby(["dataset", "normalization", "score"])["p"].transform(bh)
    write(matched_effects, "integration/matched_gene_patient_HGP_effects.tsv")
    write(pd.DataFrame(matched_scores), "integration/matched_gene_patient_HGP_scores.tsv")
    coverage_effects = pd.DataFrame(coverage_effects)
    write(coverage_effects, "integration/coverage_only_claudin_HGP.tsv")


if __name__ == '__main__':
    main()
