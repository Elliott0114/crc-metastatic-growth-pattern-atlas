"""Bounded E-MTAB-12022 epithelial co-expression extension; seed 42.

Reuse the existing reference-marker reconstruction. Patients, not cells, are
inferential units. No REC gate or new cell annotation is fitted.
"""
from pathlib import Path
from importlib.metadata import version
import hashlib
import itertools
import json
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_results/rec_hgp_focused_revision_2026-09-07/epithelial"
SOURCE = ROOT / "analysis_results/e_mtab_12022_primary_reconstruction.h5ad"
SETS = ROOT / "metadata/rec_public_upgrade_2026-09-07/analysis_gene_sets.tsv"
OUT.mkdir(parents=True, exist_ok=True)
sets = {k: set(g.gene) for k, g in pd.read_csv(SETS, sep="\t").groupby("set_id")}
a = ad.read_h5ad(SOURCE, backed="r")
mask = a.obs.cell_type.astype(str).eq("CRC/epithelial").to_numpy()
obs = a.obs.loc[mask, ["patient", "hgp", "total_umi", "n_features", "mt_fraction"]].copy()
x = sparse.csr_matrix(a.X[np.flatnonzero(mask), :], dtype=np.float64)
genes = pd.Index(a.var_names.str.upper())
assert genes.is_unique and obs.index.is_unique and x.shape[0] == len(obs)
a.file.close()
x.eliminate_zeros()
lib = np.asarray(x.sum(axis=1)).ravel()
ng = np.diff(x.indptr)
assert np.all(lib > 0) and np.allclose(lib, obs.total_umi) and np.array_equal(ng, obs.n_features)
log = x.multiply((10000/lib)[:, None]).tocsr()
np.log1p(log.data, out=log.data)
means = np.asarray(log.mean(axis=0)).ravel()
det = np.asarray((x > 0).sum(axis=0)).ravel()
pool = np.flatnonzero(det >= 5)
bins = np.full(len(genes), -1, dtype=int)
bins[pool] = pd.qcut(pd.Series(means[pool]).rank(method="first"), 25, labels=False).to_numpy()
positions = {g: i for i, g in enumerate(genes)}
overlap = sets["CANELLAS_CORE_HRC"] & sets["REACTOME_TIGHT_JUNCTION_INTERACTIONS"]
targets = {"repair": sets["CANELLAS_CORE_HRC"]-overlap-sets["WHITE_RSC"],
           "junction": sets["REACTOME_TIGHT_JUNCTION_INTERACTIONS"]-overlap-sets["WHITE_RSC"],
           "RSC": sets["WHITE_RSC"]}
defined = {k: sorted(g) for k, g in targets.items()}
targets = {k: sorted(g for g in gs if g in positions and det[positions[g]] >= 5)
           for k, gs in targets.items()}
assert all(len(gs) >= 5 for gs in targets.values())
assert not any(set(targets[k]) & set(targets[l]) for k, l in itertools.combinations(targets, 2))
forbidden = {positions[g] for gs in targets.values() for g in gs}
rng = np.random.default_rng(42)
controls, used, lock = {}, set(), []
for component, gs in targets.items():
    ctrl = []
    for b in sorted({bins[positions[g]] for g in gs}):
        choices = np.array([i for i in pool[bins[pool] == b] if i not in forbidden and i not in used])
        take = rng.choice(choices, min(50, len(choices)), replace=False)
        ctrl.extend(take.tolist())
        used.update(take.tolist())
    assert ctrl
    controls[component] = ctrl
    lock.extend(dict(component=component, role="target", gene=g) for g in gs)
    lock.extend(dict(component=component, role="control", gene=genes[i]) for i in ctrl)
pd.DataFrame(lock).to_csv(OUT / "gene_control_lock.tsv", sep="\t", index=False, encoding="utf-8")
pd.DataFrame([dict(component=k, n_defined_disjoint=len(defined[k]), n_measured=len(gs),
    missing_or_low_detection=";".join(sorted(set(defined[k])-set(gs)))) for k, gs in targets.items()]).to_csv(
    OUT / "component_coverage.tsv", sep="\t", index=False, encoding="utf-8")


def score(gs, ctrl, method):
    t = log[:, [positions[g] for g in gs]]
    if method == "matched":
        return np.asarray(t.mean(axis=1)).ravel()-np.asarray(log[:, ctrl].mean(axis=1)).ravel()
    mu = np.asarray(t.mean(axis=0)).ravel()
    sd = np.sqrt(np.maximum(np.asarray(t.power(2).mean(axis=0)).ravel()-mu**2, 0))
    assert np.all(sd > 0)
    return np.asarray(t.multiply(1/sd).mean(axis=1)).ravel()-np.mean(mu/sd)


def residual(v, cov):
    sd = np.std(cov, axis=0)
    if np.any(sd == 0):
        return None
    design = np.column_stack([np.ones(len(v)), (cov-cov.mean(axis=0))/sd])
    if np.linalg.matrix_rank(design) != design.shape[1]:
        return None
    return v-design@np.linalg.lstsq(design, v, rcond=None)[0]


def inference(rhos):
    zs = np.arctanh(np.asarray(rhos))
    n = len(zs)
    rng = np.random.default_rng(42)
    draws = np.tanh(rng.choice(zs, (10000, n)).mean(axis=1))
    lo, hi = np.quantile(draws, [.025, .975])
    signs = np.array(list(itertools.product((-1, 1), repeat=n)))
    p = np.mean(np.abs(signs@zs/n) >= abs(zs.mean())-1e-12)
    return dict(rho=np.tanh(zs.mean()), CI_low=lo, CI_high=hi, p_signflip=p,
                n_positive=int(np.sum(np.asarray(rhos)>0)), bootstrap_B=10000)


tech = np.column_stack([np.log1p(lib), np.log1p(ng)])
strict = obs.n_features.ge(300).to_numpy() & obs.mt_fraction.le(.20).to_numpy()
patient_rows, recovery = [], []
scores = {method: {k: score(gs, controls[k], method) for k, gs in targets.items()}
          for method in ("matched", "gene_z")}
for qc, selected in (("existing_reconstruction", np.ones(len(obs), dtype=bool)), ("strict_300genes_20pct_mt", strict)):
    for patient in sorted(obs.patient.astype(str).unique()):
        ii = np.flatnonzero(obs.patient.astype(str).eq(patient).to_numpy() & selected)
        hgp = obs.loc[obs.patient.astype(str).eq(patient), "hgp"].unique()
        assert len(hgp) == 1
        variable = {}
        for k, gs in targets.items():
            t = log[ii, :][:, [positions[g] for g in gs]]
            if len(ii):
                mu = np.asarray(t.mean(axis=0)).ravel()
                variance = np.asarray(t.power(2).mean(axis=0)).ravel()-mu**2
                variable[k] = int(np.sum(variance > 1e-12))
            else:
                variable[k] = 0
        reason = "eligible" if len(ii)>=20 and min(variable.values())>=5 else "fewer_than_20_cells_or_5_variable_component_genes"
        recovery.append(dict(qc=qc, patient=patient, hgp=hgp[0], n_cells=len(ii),
            n_repair_variable=variable["repair"], n_junction_variable=variable["junction"],
            n_RSC_variable=variable["RSC"], status=reason,
            median_features=np.median(ng[ii]) if len(ii) else np.nan,
            median_mt=obs.mt_fraction.iloc[ii].median() if len(ii) else np.nan))
        for method, ss in scores.items():
            for adjustment in ("technical_only", "technical_plus_RSC"):
                status, rho = reason, np.nan
                if reason == "eligible":
                    cov = tech[ii] if adjustment=="technical_only" else np.column_stack([tech[ii], ss["RSC"][ii]])
                    xx, yy = residual(ss["repair"][ii], cov), residual(ss["junction"][ii], cov)
                    if xx is None or yy is None or min(np.ptp(xx), np.ptp(yy))<=1e-12:
                        status = "rank_deficient_or_constant_residual"
                    else:
                        rho = float(stats.spearmanr(xx, yy).statistic)
                patient_rows.append(dict(qc=qc, patient=patient, hgp=hgp[0], method=method,
                    adjustment=adjustment, n_cells=len(ii), rho=rho, status=status))
patient = pd.DataFrame(patient_rows)
patient.to_csv(OUT / "patient_correlations.tsv", sep="\t", index=False, encoding="utf-8")
pd.DataFrame(recovery).to_csv(OUT / "patient_recovery.tsv", sep="\t", index=False, encoding="utf-8")
summaries, hgp_rows = [], []
for key, group in patient.groupby(["qc", "method", "adjustment"]):
    keydict = dict(zip(["qc", "method", "adjustment"], key))
    g = group[group.status.eq("eligible") & group.rho.notna()]
    n = len(g)
    row = keydict | dict(n_patients=n, status="estimable" if n>=3 else "fewer_than_3_patients")
    if n >= 3:
        assert np.all(np.abs(g.rho) < 1)
        row.update(inference(g.rho))
    summaries.append(row)
    counts = g.groupby("hgp", observed=True).size().to_dict()
    hrow = keydict | dict(n_rHGP=counts.get("rHGP", 0), n_dHGP=counts.get("dHGP", 0),
                         status="fewer_than_2_patients_per_HGP")
    if hrow["n_rHGP"]>=2 and hrow["n_dHGP"]>=2:
        zr, zd = [np.arctanh(g.loc[g.hgp.eq(h), "rho"].to_numpy()) for h in ("rHGP", "dHGP")]
        rng = np.random.default_rng(42)
        draws = rng.choice(zr, (10000, len(zr))).mean(axis=1)-rng.choice(zd, (10000, len(zd))).mean(axis=1)
        lo, hi = np.quantile(draws, [.025, .975])
        combined = np.r_[zr, zd]
        allocations = list(itertools.combinations(range(len(combined)), len(zr)))
        null = np.array([combined[list(ii)].mean()-np.delete(combined, ii).mean() for ii in allocations])
        delta = zr.mean()-zd.mean()
        hrow.update(status="exploratory", fisher_z_difference=delta, CI_low=lo, CI_high=hi,
                    p_allocation=np.mean(np.abs(null)>=abs(delta)-1e-12), allocations=len(allocations))
    hgp_rows.append(hrow)
summary = pd.DataFrame(summaries)
summary["p_BH_two_adjustments"] = summary.groupby(["qc", "method"]).p_signflip.transform(
    lambda s: stats.false_discovery_control(s) if s.notna().all() else np.full(len(s), np.nan))
summary.to_csv(OUT / "correlation_summary.tsv", sep="\t", index=False, encoding="utf-8")
pd.DataFrame(hgp_rows).to_csv(OUT / "exploratory_HGP_correlations.tsv", sep="\t", index=False, encoding="utf-8")
with SOURCE.open("rb") as handle:
    source_hash = hashlib.file_digest(handle, "sha256").hexdigest()
(OUT / "provenance.json").write_text(json.dumps(dict(source=str(SOURCE.relative_to(ROOT)), sha256=source_hash,
    annotation="Existing reference-marker CRC/epithelial reconstruction; not a REC or malignant call",
    cell_count=len(obs), gene_control_selection="All annotated epithelial cells; >=5 detected cells per gene; 25 mean-expression bins; up to 50 disjoint controls per bin; seed 42",
    qc_sensitivity="At least 300 detected genes and mitochondrial fraction <=0.20; scores/control selection frozen from all epithelial cells",
    score="Log1p counts per 10,000; disjoint HRC/junction/RSC targets; technical and RSC models use identical targets",
    residualization="Within patient: OLS on log1p library size and log1p detected genes, with or without RSC; Spearman correlation of residuals",
    inference="Equal-patient Fisher-z average, 10000 patient bootstraps; sign flips assume symmetric/exchangeable patient correlations; BH across two adjustments",
    packages={k: version(k) for k in ("anndata", "numpy", "pandas", "scipy")}), indent=2)+"\n", encoding="utf-8")
print(pd.DataFrame(recovery).to_string(index=False), flush=True)
print(summary.to_string(index=False), flush=True)
