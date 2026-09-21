"""Fixed junction components in HGP tumours and paired spatial regions."""
import itertools

import numpy as np
import pandas as pd
import statsmodels.api as sm

from rec_junction_context_common import (
    BOOT, ROOT, bh, definitions, mean_summary,
    rng_for, score_frame, write,
)


def linear_effect(frame, outcome, predictor, covariates, label, stratify=False):
    names = [predictor] + covariates
    x = np.column_stack([np.ones(len(frame)), frame[names].to_numpy(float)])
    y = frame[outcome].to_numpy(float)
    fit = sm.OLS(y, x).fit(cov_type="HC3", use_t=True)
    rng = rng_for(label)
    if stratify:
        pools = [np.flatnonzero(frame[predictor].to_numpy() == g) for g in (0, 1)]
        idx = np.concatenate([rng.choice(p, (BOOT, len(p))) for p in pools], axis=1)
    else:
        idx = rng.integers(len(frame), size=(BOOT, len(frame)))
    xb, yb = x[idx], y[idx]
    ranks = np.linalg.matrix_rank(xb)
    ok = ranks == x.shape[1]
    beta = (np.linalg.pinv(xb[ok]) @ yb[ok, :, None])[:, 1, 0]
    low, high = np.quantile(beta, [.025, .975])
    leave = [np.linalg.lstsq(np.delete(x, i, 0), np.delete(y, i), rcond=None)[0][1]
             for i in range(len(frame))]
    return {"n": len(y), "effect": fit.params[1], "low": low, "high": high,
            "p": fit.pvalues[1], "p_method": "OLS_HC3_t_approximate",
            "bootstrap_valid": int(ok.sum()), "loo_min": min(leave), "loo_max": max(leave),
            "design_condition_number": np.linalg.cond(x)}


def paired_adjusted(delta, component, covariates, label):
    # Delta = region coefficient + changes in the covariates; do not center changes.
    x = np.column_stack([np.ones(len(delta)), delta[covariates].to_numpy(float)])
    y = delta[component].to_numpy(float)
    fit = sm.OLS(y, x).fit(cov_type="HC3", use_t=True)
    idx = rng_for(label).integers(len(y), size=(BOOT, len(y)))
    xb, yb = x[idx], y[idx]
    ok = np.linalg.matrix_rank(xb) == x.shape[1]
    beta = (np.linalg.pinv(xb[ok]) @ yb[ok, :, None])[:, 0, 0]
    low, high = np.quantile(beta, [.025, .975])
    leave = [np.linalg.lstsq(np.delete(x, i, 0), np.delete(y, i), rcond=None)[0][0]
             for i in range(len(y))]
    return {"n": len(y), "effect": fit.params[0], "low": low, "high": high,
            "p": fit.pvalues[0], "p_method": "paired_change_OLS_HC3_t_approximate",
            "bootstrap_valid": int(ok.sum()), "loo_min": min(leave), "loo_max": max(leave),
            "zero_outside_covariate_range": ";".join(c for c in covariates
                if delta[c].min() > 0 or delta[c].max() < 0),
            "design_condition_number": np.linalg.cond(x)}


def main():
    sets = definitions()
    base = ROOT / "analysis_results/rec_public_upgrade_2026-09-07/specificity"
    hgp_path = base / "gse151165_logcpm.tsv"
    expr = pd.read_csv(hgp_path, sep="\t", index_col=0)
    meta = pd.read_csv(base / "gse151165_program_scores.tsv", sep="\t").set_index("sample_id")
    scores, cov_hgp, z = score_frame(expr, sets, "GSE151165", ddof=1)
    scores = scores.join(meta[["hgp", "r"]])
    assert len(scores) == 15 and scores.r.sum() == 6 and scores.index.is_unique
    write(scores.rename_axis("sample_id").reset_index(), "contexts/hgp_scores.tsv")
    allocations = np.zeros((5005, 15))
    for i, group in enumerate(itertools.combinations(range(15), 6)):
        allocations[i, list(group)] = 1 / 6
        allocations[i, [j for j in range(15) if j not in group]] = -1 / 9
    rows = []
    for comp in ["Claudin", "Polarity", "Junction", "HRC"]:
        for model, covars in [("unadjusted", []), ("HRC", ["HRC"]),
                              ("HRC_epithelial", ["HRC", "Epithelial"]),
                              ("HRC_liver", ["HRC", "Liver"])]:
            if comp == "HRC" and covars:
                continue
            result = linear_effect(scores, comp, "r", covars, "hgp" + comp + model, True)
            if not covars:
                result["p"] = float(np.mean(np.abs(allocations @ scores[comp]) >= abs(result["effect"]) - 1e-12))
                result["p_method"] = "exhaustive_5005_patient_allocations"
            rows.append(dict(dataset="GSE151165", component=comp, model=model, **result))
    gene_rows = []
    for gene in sets["Junction"]:
        if gene not in z.index or not np.isfinite(z.loc[gene]).all():
            continue
        f = scores.assign(target=z.loc[gene])
        result = linear_effect(f, "target", "r", [], "hgpgene" + gene, True)
        result["p"] = float(np.mean(np.abs(allocations @ f.target) >= abs(result["effect"]) - 1e-12))
        result["p_method"] = "exhaustive_5005_patient_allocations"
        gene_rows.append(dict(dataset="GSE151165", gene=gene, **result))

    base = ROOT / "analysis_results/rec_program_reader_facing_specificity"
    liu_path = base / "gse294385_patient_region_all_gene_counts.tsv.gz"
    counts = pd.read_csv(liu_path, sep="\t").pivot(index="gene", columns=["patient", "region"], values="raw_count")
    meta = pd.read_csv(base / "gse294385_patient_region_metadata.tsv", sep="\t").set_index(["patient", "region"])
    assert len(counts.columns) == 22
    write(counts.isna().sum().rename("n_unmeasured_genes").reset_index(),
          "contexts/regional_panel_missingness.tsv")
    assert np.allclose(counts.sum(), meta.loc[counts.columns, "library_size"])
    manifest = pd.read_csv(ROOT / "metadata/gse294385_liver_sample_manifest.tsv", sep="\t")
    manifest = manifest[manifest.selected_for_paired_extension.astype(str).str.lower() == "true"]
    common = None
    for sample in manifest["sample"]:
        path = ROOT / "data_sources/Liu_2026_GSE294385/extracted" / sample / "filtered_feature_bc_matrix/features.tsv.gz"
        measured = set(pd.read_csv(path, sep="\t", header=None, usecols=[1])[1].str.upper())
        common = measured if common is None else common & measured
    counts = counts.loc[counts.index.intersection(sorted(common))]
    assert not counts.isna().any().any()
    expr_liu = np.log2((counts + .5).div(meta.loc[counts.columns, "library_size"] + 1) * 1e6)
    scores_liu, cov_liu, z_liu = score_frame(expr_liu, sets, "GSE294385", ddof=0)
    write(scores_liu.reset_index(), "contexts/paired_region_scores.tsv")
    delta = scores_liu.xs("macro_tumour", level="region") - scores_liu.xs("micro_tumour", level="region")
    assert len(delta) == 11 and not delta.isna().any().any()
    write(delta.reset_index(), "contexts/paired_region_changes.tsv")
    for comp in ["Claudin", "Polarity", "Junction", "HRC"]:
        for model, covars in [("unadjusted", []), ("HRC", ["HRC"]),
                              ("epithelial", ["Epithelial"]), ("liver", ["Liver"]),
                              ("HRC_epithelial", ["HRC", "Epithelial"]),
                              ("HRC_liver", ["HRC", "Liver"])]:
            if comp == "HRC" and "HRC" in covars:
                continue
            if covars:
                result = paired_adjusted(delta, comp, covars, "liu" + comp + model)
            else:
                result = mean_summary(delta[comp], "liu" + comp)
                result["p_method"] = "exhaustive_2048_patient_sign_flips"
            rows.append(dict(dataset="GSE294385", component=comp, model=model, **result))
    gene_delta = z_liu.xs("macro_tumour", level="region", axis=1) - z_liu.xs("micro_tumour", level="region", axis=1)
    for gene in sets["Junction"]:
        if gene in gene_delta.index and np.isfinite(gene_delta.loc[gene]).all():
            gene_rows.append(dict(dataset="GSE294385", gene=gene,
                p_method="exhaustive_2048_patient_sign_flips",
                **mean_summary(gene_delta.loc[gene], "liugene" + gene)))

    summary = pd.DataFrame(rows)
    summary["q_within_two_components"] = np.nan
    for _, sub in summary[summary.component.isin(["Claudin", "Polarity"])].groupby(["dataset", "model"]):
        summary.loc[sub.index, "q_within_two_components"] = bh(sub.p)
    write(summary, "contexts/component_effects.tsv")
    genes = pd.DataFrame(gene_rows)
    genes["q_within_source_genes"] = genes.groupby("dataset").p.transform(bh)
    write(genes, "contexts/gene_effects.tsv")
    write(pd.concat([cov_hgp, cov_liu]), "contexts/gene_coverage.tsv")

    # A historical audit uses the original lists without the new alias correction.
    original = pd.read_csv(ROOT / "metadata/rec_public_upgrade_2026-09-07/analysis_gene_sets.tsv", sep="\t")
    audits = []
    for name in ["CANELLAS_CORE_HRC", "REACTOME_TIGHT_JUNCTION_INTERACTIONS"]:
        gg = original.loc[original.set_id == name, "gene"]
        raw_z = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1, ddof=1), axis=0)
        measured = raw_z.index.intersection(gg)
        values = raw_z.loc[measured].mean()
        # Re-read HGP metadata; the regional metadata intentionally has a different index.
        old = pd.read_csv(ROOT / "analysis_results/rec_public_upgrade_2026-09-07/specificity/gse151165_program_scores.tsv", sep="\t").set_index("sample_id")
        error = float((values - old[name]).abs().max())
        assert error < 1e-7, (name, error)
        d = float(values[old.r == 1].mean() - values[old.r == 0].mean())
        audits.append(dict(dataset="GSE151165", set_id=name, n_genes=len(measured), effect=d, max_score_error=error))
        raw_liu_z = expr_liu.sub(expr_liu.mean(axis=1), axis=0).div(expr_liu.std(axis=1, ddof=0), axis=0)
        measured = raw_liu_z.index.intersection(gg)
        values = raw_liu_z.loc[measured].mean()
        d = values.xs("macro_tumour", level="region") - values.xs("micro_tumour", level="region")
        audits.append(dict(dataset="GSE294385", set_id=name, n_genes=len(measured), effect=d.mean(), max_score_error=np.nan))
    write(pd.DataFrame(audits), "contexts/historical_reproduction.tsv")
    print(summary.loc[summary.component.isin(["Claudin", "Polarity"]), ["dataset", "component", "model", "effect", "low", "high", "p", "q_within_two_components"]].to_string(index=False))


if __name__ == "__main__":
    main()
