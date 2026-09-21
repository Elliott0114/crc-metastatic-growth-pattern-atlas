"""Fixed definitions and patient-level statistics for the junction-context study."""
from __future__ import annotations

import hashlib
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_results/rec_junction_context_2026-09-16"
DATA = ROOT / "data_sources/rec_junction_context_2026-09-16"
SPEC = ROOT / "metadata/rec_junction_context_2026-09-16.md"
SEED = 42
BOOT = 10_000


def canonical(gene):
    gene = str(gene).upper()
    return "PALS1" if gene == "MPP5" else gene


def definitions():
    source = pd.read_csv(ROOT / "metadata/rec_public_upgrade_2026-09-07/analysis_gene_sets.tsv", sep="\t")
    sets = {k: set(map(canonical, g.gene)) for k, g in source.groupby("set_id")}
    junction = sets["REACTOME_TIGHT_JUNCTION_INTERACTIONS"]
    polarity = set("CRB3 PALS1 PARD3 PARD6A PARD6B PARD6G PATJ PRKCI".split())
    claudin = {g for g in junction if g.startswith("CLDN")}
    assert len(junction) == 30 and len(claudin) == 21
    assert claudin | polarity | {"F11R"} == junction
    result = {"Claudin": claudin, "Polarity": polarity, "Junction": junction}
    mapping = {"HRC": "CANELLAS_CORE_HRC", "RSC": "WHITE_RSC",
               "iCMS2": "ICMS2_TEMPLATE_UP", "iCMS3": "ICMS3_TEMPLATE_UP",
               "E2F": "HALLMARK_E2F_TARGETS", "Hypoxia": "HALLMARK_HYPOXIA",
               "Epithelial": "PAN_EPITHELIAL_7", "Liver": "HEPATOCYTE_CONTEXT_6"}
    for name, key in mapping.items():
        result[name] = sets[key] - junction
    return {k: sorted(v) for k, v in result.items()}


def write(frame, filename):
    path = OUT / filename
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, encoding="utf-8", na_rep="NA")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def rng_for(label):
    seed = SEED + int.from_bytes(hashlib.sha256(label.encode("utf-8")).digest()[:4], "little")
    return np.random.default_rng(seed)


def bh(values):
    values = np.asarray(values, dtype=float)
    result = np.full(len(values), np.nan)
    valid = np.isfinite(values)
    if valid.any():
        result[valid] = stats.false_discovery_control(values[valid])
    return result


def sign_p(values):
    values = np.asarray(values, dtype=float)
    signs = np.asarray(list(itertools.product([-1, 1], repeat=len(values))))
    null = (signs @ values) / len(values)
    return float(np.mean(np.abs(null) >= abs(values.mean()) - 1e-12))


def mean_summary(values, label, fisher=False):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return {"n": 0, "effect": np.nan, "low": np.nan, "high": np.nan, "p": np.nan}
    y = np.arctanh(np.clip(values, -.999999, .999999)) if fisher else values
    draws = rng_for(label).choice(y, size=(BOOT, len(y)), replace=True).mean(axis=1)
    transform = np.tanh if fisher else lambda x: x
    low, high = transform(np.quantile(draws, [.025, .975]))
    leave = [(y.sum() - v) / (len(y) - 1) for v in y] if len(y) > 1 else [np.nan]
    return {"n": len(y), "effect": float(transform(y.mean())), "low": low, "high": high,
            "p": sign_p(y) if len(y) >= 5 else np.nan, "n_positive": int((y > 0).sum()),
            "loo_min": float(transform(np.min(leave))), "loo_max": float(transform(np.max(leave)))}


def rank_residual(values, design):
    ranks = stats.rankdata(values, axis=0)
    return ranks - design @ np.linalg.lstsq(design, ranks, rcond=None)[0]


def rank_design(frame, covariates, state=False):
    columns = [np.ones(len(frame))]
    for name in covariates:
        values = stats.rankdata(frame[name].to_numpy(float))
        if values.std() > 0:
            columns.append((values - values.mean()) / values.std())
    if state:
        indicators = pd.get_dummies(frame.state, drop_first=True, dtype=float).to_numpy()
        columns.extend(indicators.T)
    return np.column_stack(columns)


def correlations(residuals, x):
    x = x - x.mean()
    residuals = residuals - residuals.mean(axis=0)
    denom = np.sqrt((residuals ** 2).sum(axis=0) * (x ** 2).sum())
    return np.divide(x @ residuals, denom, out=np.full(residuals.shape[1], np.nan), where=denom > 1e-12)


def score_frame(expression, sets, dataset, ddof=1):
    """Genes x profiles log expression; one row per canonical gene."""
    expression = expression.copy()
    raw_names = dict(zip(map(canonical, expression.index), expression.index.astype(str)))
    expression.index = [canonical(g) for g in expression.index]
    if expression.index.has_duplicates:
        raise ValueError(f"Duplicate canonical expression identifiers: {dataset}")
    sd = expression.std(axis=1, ddof=ddof)
    z = expression.sub(expression.mean(axis=1), axis=0).div(sd.replace(0, np.nan), axis=0)
    scores = pd.DataFrame(index=expression.columns)
    coverage = []
    for name, genes in sets.items():
        measured = [g for g in genes if g in z.index and np.isfinite(z.loc[g]).all()]
        scores[name] = z.loc[measured].mean(axis=0) if measured else np.nan
        for g in genes:
            coverage.append({"dataset": dataset, "component": name, "gene": g,
                             "source_gene": raw_names.get(g, ""), "measured_variable": g in measured})
    scores["iCMS"] = scores.iCMS3 - scores.iCMS2
    return scores, pd.DataFrame(coverage), z
