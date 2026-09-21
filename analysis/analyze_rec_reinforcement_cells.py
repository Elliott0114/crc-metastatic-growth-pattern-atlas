"""Target-membership checks with fixed patient conditioning and no outcome selection."""
from pathlib import Path
import itertools
import json

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.optimize import linear_sum_assignment

from rec_junction_context_common import ROOT, canonical, definitions, digest, rank_design, rank_residual, correlations, mean_summary, bh

OUT = ROOT / 'analysis_results/rec_manuscript_reinforcement_2026-09-17/cells'
BASE = ROOT / 'analysis_results/rec_junction_context_2026-09-16/cells'
SPEC = ROOT / 'metadata/rec_manuscript_reinforcement_2026-09-17.md'


def write(frame, name):
    frame.to_csv(OUT / name, sep='\t', index=False, encoding='utf-8', na_rep='NA')


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cell = pd.read_csv(BASE / 'cell_scores.tsv.gz', sep='\t')
    eligible = cell.patient.value_counts().loc[lambda x: x >= 100].index
    cell = cell[cell.patient.isin(eligible)].reset_index(drop=True)
    assert len(cell) == 15987 and len(eligible) == 14
    lock = pd.read_csv(BASE / 'gene_control_lock.tsv', sep='\t')
    groups = {c: sorted(lock.loc[(lock.component == c) & (lock.role == 'target') & lock.measured, 'gene']) for c in ['Claudin', 'Polarity']}
    genes = groups['Claudin'] + groups['Polarity']
    a = ad.read_h5ad(ROOT / 'data_sources/Ogden_2025_CRLM_multiome/CRCLM_multiome_GEX_decontaminated.h5ad')
    ids = a.obs_names.get_indexer(cell.cell_id)
    assert (ids >= 0).all() and np.array_equal(a.obs.iloc[ids].Patient.astype(str), cell.patient)
    names = pd.Index([canonical(g) for g in a.var_names])
    x = sparse.csr_matrix(a.X[ids], dtype=float)
    library = np.asarray(x.sum(axis=1)).ravel()
    values = x[:, names.get_indexer(genes)].toarray()
    values = np.log1p(values * (10000 / library[:, None]))
    del a, x
    for c in groups:
        ix = [genes.index(g) for g in groups[c]]
        np.testing.assert_allclose(values[:, ix].mean(axis=1), cell[c + '_raw'], atol=1e-12)
    means, detection = [], []
    for _, ix in cell.groupby('patient', sort=True).groups.items():
        means.append(values[ix].mean(axis=0))
        detection.append((values[ix] > 0).mean(axis=0))
    features = pd.DataFrame({'gene': genes, 'component': ['Claudin'] * 20 + ['Polarity'] * 8,
                             'mean_log_expression': np.mean(means, axis=0), 'detection_fraction': np.mean(detection, axis=0)})
    f = np.column_stack([np.log(np.maximum(features.mean_log_expression, 1e-6)),
                         stats.logistic.ppf(np.clip(features.detection_fraction, 1e-6, 1-1e-6))])
    f = (f - f.mean(axis=0)) / f.std(axis=0)
    features['z_log_mean'], features['z_logit_detection'] = f.T
    delta = f[20:, None, :] - f[None, :20, :]
    allowed = np.max(np.abs(delta), axis=2) <= .5
    distance = (delta ** 2).sum(axis=2)
    cost = np.concatenate([np.where(allowed, distance, 10000), np.full((8, 8), 100.)], axis=1)
    ri, ci = linear_sum_assignment(cost)
    pairs = [(int(r), int(c)) for r, c in zip(ri, ci) if c < 20 and allowed[r, c]]
    pairtable = pd.DataFrame([{'polarity_gene': genes[20+r], 'claudin_gene': genes[c],
                               'squared_distance': distance[r, c], 'delta_z_log_mean': delta[r, c, 0],
                               'delta_z_logit_detection': delta[r, c, 1]} for r, c in pairs])
    features['matched'] = features.gene.isin([genes[20+r] for r, _ in pairs] + [genes[c] for _, c in pairs])
    write(features, 'matching_features.tsv'); write(pairtable, 'measurement_matching_pairs.tsv')
    rng = np.random.default_rng(20260917)
    selections = set()
    while len(selections) < 500:
        selections.add(tuple(sorted(rng.choice(20, 8, replace=False).tolist())))
    subsets = sorted(selections)
    write(pd.DataFrame([{'subset': k, 'gene': genes[j]} for k, sub in enumerate(subsets) for j in sub]), 'size_subset_members.tsv')
    # Record all memberships before computing HRC-conditioned target effects.
    (OUT / 'selection_lock.json').write_text(json.dumps({'spec_sha256': digest(SPEC), 'n_pairs': len(pairs),
        'seed': 20260917, 'n_subsets': len(subsets), 'selection_uses_outcomes': False}, indent=2)+'\n', encoding='utf-8')
    target_c = np.column_stack([values[:, :20].mean(axis=1)] + [values[:, list(sub)].mean(axis=1) for sub in subsets]
                              + ([values[:, [c for _, c in pairs]].mean(axis=1)] if pairs else []))
    target_p = np.column_stack([values[:, 20:].mean(axis=1)] * 501
                              + ([values[:, [20+r for r, _ in pairs]].mean(axis=1)] if pairs else []))
    scopes = ['full_reference'] + ['size_matched'] * 500 + (['measurement_matched'] if pairs else [])
    result, errors = [], []
    old = pd.read_csv(BASE / 'patient_partial_correlations.tsv', sep='\t')
    old = old[(old.scope == 'All_epithelial') & (old.model == 'conditional') & (old.scoring == 'matched')]
    for patient, ix in cell.groupby('patient', sort=True).groups.items():
        g = cell.loc[ix]
        design = rank_design(g, ['log_library', 'log_genes', 'RSC', 'iCMS', 'E2F'], state=True)
        hr = rank_residual(g[['HRC']].to_numpy(), design)[:, 0]
        for mode in ['fixed_background', 'target_mean']:
            c, p = target_c[ix].copy(), target_p[ix].copy()
            if mode == 'fixed_background':
                c -= (g.Claudin_raw-g.Claudin).to_numpy()[:, None]
                p -= (g.Polarity_raw-g.Polarity).to_numpy()[:, None]
            rc = correlations(rank_residual(c, design), hr)
            rp = correlations(rank_residual(p, design), hr)
            dz = np.arctanh(np.clip(rc, -.999999, .999999))-np.arctanh(np.clip(rp, -.999999, .999999))
            for j in range(len(rc)):
                result.append({'patient': patient, 'n_cells': len(g), 'scoring': mode, 'scope': scopes[j],
                               'subset': j-1 if 0 < j <= 500 else -1, 'claudin_rho': rc[j], 'polarity_rho': rp[j], 'difference_z': dz[j]})
            if mode == 'fixed_background':
                ref = old[old.patient == patient].set_index('component').rho
                errors.extend([abs(rc[0]-ref.Claudin), abs(rp[0]-ref.Polarity)])
    result = pd.DataFrame(result)
    write(result, 'patient_subset_correlations.tsv.gz')
    assert max(errors) < 1e-10, max(errors)
    averaged = result.groupby(['scoring', 'scope', 'patient'], as_index=False).difference_z.mean()
    write(averaged, 'patient_summary_differences.tsv')
    summaries = []
    for (mode, scope), d in averaged.groupby(['scoring', 'scope']):
        row = {'scoring': mode, 'scope': scope, 'n_target_pairs': len(pairs) if scope == 'measurement_matched' else 8 if scope == 'size_matched' else 0,
               **mean_summary(d.difference_z, 'reinforcement_'+mode+scope)}
        if scope == 'measurement_matched' and len(pairs) < 3:
            row['p'] = np.nan
        summaries.append(row)
    summary = pd.DataFrame(summaries)
    test = summary.scope.ne('full_reference')
    summary['q_four_secondary_tests'] = np.nan
    summary.loc[test, 'q_four_secondary_tests'] = bh(summary.loc[test, 'p'])
    write(summary, 'comparison_summary.tsv')
    between = result[result.scope == 'size_matched'].groupby(['scoring', 'subset']).difference_z.mean().reset_index()
    write(between, 'subset_mean_differences.tsv')
    write(between.groupby('scoring').difference_z.agg(n='size', mean='mean', median='median', minimum='min', maximum='max',
          perturbation_q025=lambda s:s.quantile(.025), perturbation_q975=lambda s:s.quantile(.975), fraction_positive=lambda s:(s>0).mean()).reset_index(), 'subset_perturbation_summary.tsv')
    (OUT / 'validation.json').write_text(json.dumps({'n_patients': 14, 'n_cells': len(cell), 'full_reference_max_abs_error': max(errors),
        'selected_pairs': len(pairs), 'all_pairs_within_caliper': bool(all(allowed[r,c] for r,c in pairs)),
        'spec_sha256': digest(SPEC)}, indent=2)+'\n', encoding='utf-8')
    print(pairtable.to_string(index=False), flush=True)
    print(summary.to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
