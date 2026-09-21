"""Resolve the existing junction programme within source-labelled epithelial cells."""
from __future__ import annotations

import json
from importlib.metadata import version

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from rec_junction_context_common import (ROOT, OUT, SPEC, canonical, definitions, write,
    digest, rng_for, bh, mean_summary, rank_design, rank_residual, correlations)


def main():
    sets = definitions()
    write(pd.DataFrame([{"component": k, "gene": g} for k, genes in sets.items() for g in genes]), "definitions.tsv")
    path = ROOT / "data_sources/Ogden_2025_CRLM_multiome/CRCLM_multiome_GEX_decontaminated.h5ad"
    a = ad.read_h5ad(path)
    names = pd.Index([canonical(g) for g in a.var_names])
    assert not names.has_duplicates
    x = sparse.csr_matrix(a.X, dtype=np.float64)
    assert np.isfinite(x.data).all() and x.data.min() >= 0
    libraries = np.asarray(x.sum(axis=1)).ravel()
    detected = np.asarray((x > 0).sum(axis=1)).ravel()
    valid = libraries > 0
    obs = a.obs.loc[valid, ["Patient", "Cell_type", "Cell_subtype", "Therapy"]].copy()
    obs.columns = ["patient", "compartment", "state", "therapy"]
    x = x[valid]
    libraries = libraries[valid]
    detected = detected[valid]
    norm = x.multiply(10_000 / libraries[:, None]).tocsr()
    norm.data = np.log1p(norm.data)
    jpos = names.get_indexer(sets["Junction"])
    present = jpos >= 0
    jgenes = np.asarray(sets["Junction"])[present]
    junction = norm[:, jpos[present]].toarray()
    origin = []
    for (patient, compartment), indices in obs.reset_index(drop=True).groupby(["patient", "compartment"], observed=True).groups.items():
        values = junction[np.asarray(indices)]
        for g, mean, detection in zip(jgenes, values.mean(axis=0), (values > 0).mean(axis=0)):
            origin.append(dict(patient=patient, compartment=compartment, n_cells=len(values), gene=g,
                               mean_log_expression=mean, detection_fraction=detection))
    write(pd.DataFrame(origin), "cells/patient_compartment_genes.tsv")
    ep = obs.compartment.eq("Epithelial").to_numpy()
    m = norm[ep]
    cell = obs.loc[ep].reset_index(names="cell_id")
    cell["log_library"] = np.log1p(libraries[ep])
    cell["log_genes"] = np.log1p(detected[ep])
    gene_means = np.asarray(m.mean(axis=0)).ravel()
    ordered = np.flatnonzero(gene_means > 0)
    ordered = ordered[np.argsort(gene_means[ordered], kind="stable")]
    bins = np.full(len(names), -1)
    bins[ordered] = np.minimum(np.arange(len(ordered)) * 25 // len(ordered), 24)
    excluded = set().union(*map(set, sets.values()))
    available = np.asarray([i for i in ordered if names[i] not in excluded])
    used_controls = set()
    locks = []
    for component, genes in sets.items():
        positions = [names.get_loc(g) for g in genes if g in names and gene_means[names.get_loc(g)] > 0]
        assert positions, component
        controls = []
        for bin_id in sorted(set(bins[positions])):
            pool = [i for i in available if bins[i] == bin_id and i not in used_controls]
            selected = rng_for(component + str(bin_id)).choice(pool, min(20, len(pool)), replace=False)
            controls.extend(selected)
            used_controls.update(selected)
        cell[component + "_raw"] = np.asarray(m[:, positions].mean(axis=1)).ravel()
        cell[component] = cell[component + "_raw"] - np.asarray(m[:, controls].mean(axis=1)).ravel()
        for g in genes:
            idx = names.get_loc(g) if g in names else -1
            locks.append(dict(component=component, gene=g, role="target", source_gene=str(a.var_names[idx]) if idx >= 0 else "",
                              measured=idx in positions))
        locks.extend(dict(component=component, gene=names[i], role="control", source_gene=str(a.var_names[i]), measured=True) for i in controls)
    cell["iCMS"] = cell.iCMS3 - cell.iCMS2
    cell["iCMS_raw"] = cell.iCMS3_raw - cell.iCMS2_raw
    write(pd.DataFrame(locks), "cells/gene_control_lock.tsv")
    write(cell, "cells/cell_scores.tsv.gz")
    write(cell.groupby(["patient", "state"], observed=True).size().reset_index(name="n_cells"), "cells/patient_state_inventory.tsv")
    jep = junction[ep]
    rows, gene_rows = [], []
    for scope in ["All_epithelial", "REC"]:
        eligible = cell if scope == "All_epithelial" else cell.loc[cell.state.eq("REC")]
        for patient, group in eligible.groupby("patient", observed=True):
            if len(group) < (100 if scope == "All_epithelial" else 30):
                continue
            for scoring in ["matched", "raw"]:
                suffix = "" if scoring == "matched" else "_raw"
                for model in ["technical", "conditional"]:
                    covars = ["log_library", "log_genes"]
                    if model == "conditional":
                        covars += [v + suffix for v in ["RSC", "iCMS", "E2F"]]
                    design = rank_design(group, covars, state=scope == "All_epithelial")
                    values = group[[v + suffix for v in ["HRC", "Claudin", "Polarity"]]].to_numpy()
                    residuals = rank_residual(values, design)
                    rho = correlations(residuals[:, 1:], residuals[:, 0])
                    for component, value in zip(["Claudin", "Polarity"], rho):
                        rows.append(dict(scope=scope, patient=patient, n_cells=len(group), scoring=scoring,
                                         model=model, component=component, rho=value, design_rank=np.linalg.matrix_rank(design)))
                    if scoring == "matched" and model == "conditional":
                        gene_values = jep[group.index]
                        r = correlations(rank_residual(gene_values, design), residuals[:, 0])
                        for gene, value, n_positive in zip(jgenes, r, (gene_values > 0).sum(axis=0)):
                            gene_rows.append(dict(scope=scope, patient=patient, gene=gene, n_cells=len(group),
                                                  n_detected=n_positive, rho=value if n_positive >= 5 else np.nan))
    patient = pd.DataFrame(rows)
    write(patient, "cells/patient_partial_correlations.tsv")
    summaries = []
    for keys, group in patient.groupby(["scope", "scoring", "model", "component"], observed=True):
        summaries.append(dict(zip(["scope", "scoring", "model", "component"], keys)) |
                         mean_summary(group.rho, str(keys), fisher=True))
    for keys, group in patient.groupby(["scope", "scoring", "model"], observed=True):
        pivot = group.pivot(index="patient", columns="component", values="rho").dropna()
        delta = np.arctanh(pivot.Polarity.clip(-.999999, .999999)) - np.arctanh(pivot.Claudin.clip(-.999999, .999999))
        summaries.append(dict(zip(["scope", "scoring", "model"], keys)) |
                         {"component": "Polarity_minus_Claudin_Fisher_z"} | mean_summary(delta, str(keys)))
    summary = pd.DataFrame(summaries)
    summary["q_within_three_tests"] = summary.groupby(["scope", "scoring", "model"], observed=True).p.transform(bh)
    write(summary, "cells/component_summary.tsv")
    genes = pd.DataFrame(gene_rows)
    write(genes, "cells/patient_gene_partial_correlations.tsv")
    summary_genes = []
    for (scope, gene), group in genes.groupby(["scope", "gene"], observed=True):
        summary_genes.append(dict(scope=scope, gene=gene) | mean_summary(group.rho, scope + gene, fisher=True))
    gs = pd.DataFrame(summary_genes)
    gs["q_genes"] = gs.groupby("scope", observed=True).p.transform(bh)
    write(gs, "cells/gene_summary.tsv")
    provenance = {"source_sha256": digest(path), "spec_sha256": digest(SPEC),
                  "n_source_cells": a.n_obs, "n_epithelial_cells": len(cell),
                  "versions": {p: version(p) for p in ["anndata", "numpy", "pandas", "scipy"]}}
    (OUT / "cells/provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(summary.loc[summary.scoring.eq("matched") & summary.model.eq("conditional")].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
