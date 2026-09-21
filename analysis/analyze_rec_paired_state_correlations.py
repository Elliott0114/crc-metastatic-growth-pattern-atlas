"""Paired state-correlation inference, 2026-09-07, seed 42.

Frozen scores and patient correlations are reused. Versions are recorded.
"""
from pathlib import Path
from importlib.metadata import version
import hashlib
import itertools
import json
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_results/rec_hgp_focused_revision_2026-09-07/paired_states"
OUT.mkdir(parents=True, exist_ok=True)
SOURCE = ROOT / "analysis_results/rec_public_upgrade_2026-09-07/within_cell/patient_conditional_correlations.tsv"
data = pd.read_csv(SOURCE, sep="\t")
rows, patient_rows = [], []
for threshold in (20, 10, 30):
    for method in ("matched", "gene_z"):
        selected = data[data.n_cells.ge(threshold) & data.method.eq(method)
                        & data.adjustment.eq("technical_plus_rival")]
        for rival, group in selected.groupby("rival"):
            assert not group.duplicated(["patient", "state"]).any()
            wide = group.pivot(index="patient", columns="state", values="rho")
            for comparator in ("Hypoxia", "UPR", "iREC"):
                pair = wide[["REC", comparator]].dropna()
                assert np.all(np.abs(pair.to_numpy()) < 1)
                diff = np.arctanh(pair.REC) - np.arctanh(pair[comparator])
                n = len(pair)
                key = dict(minimum_cells=threshold, method=method, rival=rival, comparator=comparator)
                rng = np.random.default_rng(42)
                draws = rng.choice(diff.to_numpy(), (10000, n)).mean(axis=1)
                low, high = np.quantile(draws, [.025, .975])
                signs = np.array(list(itertools.product((-1, 1), repeat=n)))
                p = np.mean(np.abs(signs @ diff.to_numpy() / n) >= abs(diff.mean()) - 1e-12)
                rows.append(key | dict(n_pairs=n, fisher_z_difference=diff.mean(), CI_low=low, CI_high=high,
                    mean_rho_difference=(pair.REC-pair[comparator]).mean(), n_REC_greater=int(diff.gt(0).sum()),
                    p_signflip=p, allocations=len(signs), bootstrap_B=10000, seed=42,
                    analysis_role="primary_extension" if threshold==20 and method=="matched" else "sensitivity"))
                for patient in pair.index:
                    patient_rows.append(key | dict(patient=patient, REC_rho=pair.loc[patient,"REC"],
                        comparator_rho=pair.loc[patient,comparator], fisher_z_difference=diff.loc[patient]))
results = pd.DataFrame(rows)
results["p_BH_six_tests"] = results.groupby(["minimum_cells", "method"]).p_signflip.transform(false_discovery_control)
assert results.groupby(["minimum_cells", "method"]).size().eq(6).all()
results.to_csv(OUT / "paired_state_contrasts.tsv", sep="\t", index=False, encoding="utf-8")
pd.DataFrame(patient_rows).to_csv(OUT / "paired_patient_correlations.tsv", sep="\t", index=False, encoding="utf-8")
(OUT / "provenance.json").write_text(json.dumps(dict(source=str(SOURCE.relative_to(ROOT)),
    sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(), seed=42,
    estimand="Equal-patient paired Fisher-z difference; sign-flip inference assumes symmetric/exchangeable paired differences",
    packages={k:version(k) for k in ("numpy","pandas","scipy")}), indent=2)+"\n", encoding="utf-8")
print(results[results.analysis_role.eq("primary_extension")].to_string(index=False))
