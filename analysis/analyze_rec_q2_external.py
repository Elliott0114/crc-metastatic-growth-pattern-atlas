"""Locked cross-study patient-level component validation; no outcome-selected sources."""
from pathlib import Path
import hashlib
import itertools
import json

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats

from rec_junction_context_common import canonical, definitions

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'analysis_results/rec_q2_integration_2026-09-19/external'
SPEC = ROOT / 'metadata/rec_q2_integration_2026-09-19.md'
ATLAS = ROOT / 'data_sources/CRC_Atlas_CZI_core/crc_atlas_core_czi.h5ad'
OGDEN = ROOT / 'data_sources/Ogden_2025_CRLM_multiome/CRCLM_multiome_GEX_decontaminated.h5ad'
STUDIES = ['Che_2021_Cell_Discov', 'Liu_2024_Cancer_Res', 'Wang_2023_Sci_Adv',
           'Giguelay_2022_Theranostics', 'Sathe_2023_Clin_Cancer_Res']
PROGRAMS = ['HRC', 'Claudin', 'Polarity', 'RSC', 'iCMS2', 'iCMS3', 'E2F']
SEED = 20260919


def write(d, name):
    d.to_csv(OUT / name, sep='\t', index=False, encoding='utf-8', na_rep='NA')


def read_h5(x):
    if isinstance(x, h5py.Group):
        cat = x['categories'].asstr()[:]
        code = x['codes'][:]
        return np.where(code >= 0, cat[np.maximum(code, 0)], '')
    return x.asstr()[:]


def seed_for(label):
    return SEED + int.from_bytes(hashlib.sha256(label.encode()).digest()[:4], 'little')


def summary(values, label):
    x = np.asarray(values, float)
    x = x[np.isfinite(x)]
    n = len(x)
    if not n:
        return dict(n_patients=0, estimate=np.nan, variance=np.nan, low=np.nan, high=np.nan, p_signflip=np.nan)
    if n == 1:
        return dict(n_patients=1, estimate=float(x[0]), variance=np.nan, low=np.nan,
                    high=np.nan, p_signflip=np.nan, n_positive=int(x[0] > 0))
    lo, hi = np.quantile(np.random.default_rng(seed_for(label)).choice(x, size=(10000, n)).mean(1), [.025, .975])
    if n <= 20:
        null = np.zeros(1)
        for v in x:
            null = np.concatenate([null + v, null - v])
        p = float(np.mean(np.abs(null) >= abs(x.sum()) - 1e-12))
    else:
        p = np.nan
    return dict(n_patients=n, estimate=float(x.mean()), variance=float(x.var(ddof=1) / n) if n > 1 else np.nan,
                low=float(lo), high=float(hi), p_signflip=p, n_positive=int((x > 0).sum()))


def design_matrix(frame, enhanced=False):
    cols = [np.ones(len(frame))]
    names = ['Intercept']
    for k in ['log_library', 'log_genes', 'RSC', 'iCMS', 'E2F']:
        v = stats.rankdata(frame[k].to_numpy(float))
        if v.std() > 0:
            cols.append((v - v.mean()) / v.std()); names.append(k)
    for field in ['sample'] + (['state'] if enhanced else []):
        dummy = pd.get_dummies(frame[field], prefix=field, drop_first=True, dtype=float)
        cols.extend(dummy.to_numpy().T); names.extend(dummy.columns)
    d = np.column_stack(cols)
    return d, names


def residual_rank(values, design):
    r = stats.rankdata(values, axis=0)
    residual = r - design @ np.linalg.lstsq(design, r, rcond=None)[0]
    # Exact constants (or completely explained ranks) have no correlation.
    scale = np.maximum(1., np.linalg.norm(r - r.mean(0), axis=0))
    residual[:, np.linalg.norm(residual, axis=0) < 1e-10 * scale] = 0
    return residual


def corr(x, y):
    x = x - x.mean(0); y = y - y.mean()
    den = np.sqrt(np.sum(x*x, axis=0) * np.sum(y*y))
    return np.divide(y @ x, den, out=np.full(x.shape[1], np.nan), where=den > 1e-10)


def extract_atlas(f, obs, genes, names):
    """Read contiguous CSR blocks; all features determine libraries/detection."""
    target_index = {g:i for i,g in enumerate(genes)}
    rows = np.asarray([i for i,g in enumerate(names) if g in target_index])
    columns = np.asarray([target_index[names[i]] for i in rows])
    projection = sparse.csr_matrix((np.ones(len(rows)), (rows,columns)),shape=(len(names),len(genes)))
    assert (np.asarray(projection.sum(0)).ravel()>0).all()
    ix = obs.atlas_row.to_numpy(int)
    order = np.argsort(ix); ix = ix[order]
    obs = obs.iloc[order].reset_index(drop=True)
    g = f['raw/X']; pointer = g['indptr'][:]
    result = np.empty((len(ix), len(genes)), dtype=np.float32)
    libs = np.empty(len(ix)); det = np.empty(len(ix), dtype=int)
    maxerr = 0.
    start = 0
    while start < len(ix):
        end = start + 1
        while end < len(ix) and ix[end] - ix[start] < 3000 and ix[end] - ix[end-1] < 150:
            end += 1
        first, last = int(ix[start]), int(ix[end-1]) + 1
        a, b = int(pointer[first]), int(pointer[last])
        vals = g['data'][a:b].astype(np.float64)
        assert np.isfinite(vals).all() and (vals >= 0).all()
        err = float(np.max(np.abs(vals - np.rint(vals)), initial=0)); maxerr = max(maxerr, err)
        assert err < 1e-6, 'Selected primary source does not contain raw counts'
        block = sparse.csr_matrix((vals, g['indices'][a:b], pointer[first:last+1] - a), shape=(last-first, len(names)))
        x = block[ix[start:end] - first]
        lib = np.asarray(x.sum(1)).ravel(); ng = np.asarray((x > 0).sum(1)).ravel()
        assert (lib > 0).all()
        result[start:end] = np.log1p((x @ projection).toarray() * (1e4 / lib[:, None])).astype(np.float32)
        libs[start:end] = lib; det[start:end] = ng
        start = end
        if start % 5000 < 3000:
            print('Atlas extraction', start, '/', len(ix), flush=True)
    obs['library_size'] = libs; obs['detected_genes'] = det
    return obs, result, maxerr


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sets = definitions()
    old = pd.read_csv(ROOT/'analysis_results/rec_junction_context_2026-09-16/cells/gene_control_lock.tsv', sep='\t')
    og = ad.read_h5ad(OGDEN, backed='r')
    onames = pd.Index([canonical(g) for g in og.var_names]); assert not onames.has_duplicates
    ep = og.obs.Cell_type.astype(str).eq('Epithelial')
    count = og.obs.loc[ep].Patient.value_counts()
    oi = np.flatnonzero((ep & og.obs.Patient.isin(count[count >= 100].index)).to_numpy())
    oobs = og.obs.iloc[oi][['Patient','Sample','Cell_subtype','Therapy']].astype(str).copy()
    oobs.columns = ['patient','sample','state','therapy']; oobs['study'] = 'Ogden_2025'; oobs['cell_id'] = og.obs_names[oi]
    oobs = oobs.reset_index(drop=True)
    with h5py.File(ATLAS, 'r') as f:
        names = pd.Index([canonical(g) for g in read_h5(f['var/feature_name'])])
        # Repeated symbols (including PAR copies) are summed in raw-count space.
        symbol_counts = pd.Series(names).value_counts()
        write(symbol_counts[symbol_counts>1].rename_axis('gene').reset_index(name='n_features'),'duplicated_symbol_accounting.tsv')
        mask = np.ones(f['obs/n_counts'].shape[0], dtype=bool)
        for field, value in [('tissue','liver'),('sample_type','metastasis'),('cell_type','malignant cell')]:
            cats = f['obs/'+field+'/categories'].asstr()[:]
            mask &= f['obs/'+field+'/codes'][:] == np.flatnonzero(cats == value)[0]
        ids = np.flatnonzero(mask)
        fields = dict(study='study_id',patient='donor_id',sample='sample_id',matrix_type='matrix_type',
                      original_annotation='cell_type_coarse_original_study',therapy='treatment_status_before_resection',
                      assay='assay',doi='study_doi',bioproject='NCBI_BioProject_accession')
        aobs = pd.DataFrame({k: read_h5(f['obs/'+v])[ids] for k,v in fields.items()})
        aobs['atlas_row'] = ids
        # Cell IDs are pulled only after restricting to the five sources.
        inventory = aobs.groupby(['study','patient','matrix_type'], observed=True).agg(n_cells=('atlas_row','size'), n_samples=('sample','nunique')).reset_index()
        inventory['eligible'] = inventory.study.isin(STUDIES) & inventory.matrix_type.eq('raw counts') & inventory.n_cells.ge(100)
        inventory['reason'] = np.select([~inventory.study.isin(STUDIES),~inventory.matrix_type.eq('raw counts'),inventory.n_cells.lt(100)],
                                       ['outside_locked_five_sources','not_raw_counts','fewer_than_100_malignant_cells'], default='eligible')
        write(inventory,'patient_eligibility.tsv')
        keep = inventory.loc[inventory.eligible, ['study','patient']]
        aobs = aobs.merge(keep, on=['study','patient'], how='inner', validate='many_to_one')
        assert len(keep) == 28 and set(aobs.study) == set(STUDIES)
        full_ids = read_h5(f['obs/_index'])
        aobs['cell_id'] = full_ids[aobs.atlas_row.to_numpy()]
        del full_ids
        assert aobs.cell_id.is_unique
        coverage, locked = [], {}
        for component in PROGRAMS:
            locked[component] = {}
            for role in ['target','control']:
                candidates = old.loc[(old.component==component) & (old.role==role)].copy()
                candidates['common_feature'] = candidates.gene.isin(names) & candidates.gene.isin(onames)
                candidates['retained'] = candidates.common_feature & candidates.measured.astype(bool)
                coverage.append(candidates)
                locked[component][role] = sorted(candidates.loc[candidates.retained,'gene'])
                assert locked[component][role], (component,role)
        lock = pd.concat(coverage, ignore_index=True)
        write(lock,'common_gene_control_lock.tsv')
        members = sets['Junction']
        all_genes = sorted(set(members).intersection(names).intersection(onames) | set(lock.loc[lock.retained,'gene']))
        gene_positions = {g:i for i,g in enumerate(all_genes)}
        write(pd.DataFrame({'gene':all_genes}), 'extracted_gene_order.tsv')
        write(pd.DataFrame([dict(study=s,doi=';'.join(sorted(g.doi.unique())),bioproject=';'.join(sorted(g.bioproject.unique())),
                                assays=';'.join(sorted(g.assay.unique())),n_patients=g.patient.nunique(),n_cells=len(g))
                           for s,g in aobs.groupby('study')]), 'source_registry.tsv')
        lock_record = {'spec_sha256':hashlib.sha256(SPEC.read_bytes()).hexdigest(), 'sources':STUDIES,
                       'n_candidate_patients':28,'n_common_genes':len(all_genes),'outcome_selection':False,
                       'target_sizes':{k:len(v['target']) for k,v in locked.items()},
                       'control_sizes':{k:len(v['control']) for k,v in locked.items()}}
        (OUT/'definition_lock.json').write_text(json.dumps(lock_record,indent=2)+'\n',encoding='utf-8')
        aobs, ax, integer_error = extract_atlas(f, aobs, all_genes, names)
    ox = sparse.csr_matrix(og.X[oi], dtype=float)
    olibrary = np.asarray(ox.sum(1)).ravel(); odetection = np.asarray((ox>0).sum(1)).ravel()
    omatrix = np.log1p(ox[:,onames.get_indexer(all_genes)].toarray()*(1e4/olibrary[:,None])).astype(np.float32)
    oobs['library_size'] = olibrary; oobs['detected_genes'] = odetection
    og.file.close(); del og, ox
    write(aobs,'external_cell_metadata.tsv.gz'); write(oobs,'ogden_cell_metadata.tsv.gz')
    np.savez_compressed(OUT/'selected_expression.npz',external=ax,ogden=omatrix,genes=np.asarray(all_genes))
    all_scores, patient_results, gene_results, diagnostics, detection, numerical_checks = [], [], [], [], [], []
    for obs, expr in [(aobs,ax),(oobs,omatrix)]:
        for (study,patient), idx in obs.groupby(['study','patient']).groups.items():
            for gene in members:
                values = expr[idx,gene_positions[gene]] if gene in gene_positions else None
                detection.append(dict(study=study,patient=patient,gene=gene,n_cells=len(idx),
                                      detection_fraction=float(np.mean(values>0)) if values is not None else np.nan,
                                      mean_log_expression=float(values.mean()) if values is not None else np.nan))
        for scoring in ['fixed_background','target_mean']:
            scored = obs.copy()
            for component in PROGRAMS:
                target = expr[:,[gene_positions[g] for g in locked[component]['target']]].mean(1,dtype=np.float64)
                background = expr[:,[gene_positions[g] for g in locked[component]['control']]].mean(1,dtype=np.float64)
                scored[component] = target - background if scoring=='fixed_background' else target
            scored['iCMS'] = scored.iCMS3 - scored.iCMS2
            scored['log_library'] = np.log1p(scored.library_size)
            scored['log_genes'] = np.log1p(scored.detected_genes)
            scored['scoring'] = scoring
            all_scores.append(scored)
            for (study,patient), idx in scored.groupby(['study','patient']).groups.items():
                g = scored.loc[idx]
                for model in ['common'] + (['state_enhanced'] if study=='Ogden_2025' else []):
                    design, covariates = design_matrix(g, model=='state_enhanced')
                    rank = np.linalg.matrix_rank(design)
                    residual = residual_rank(g[['HRC','Claudin','Polarity']].to_numpy(),design)
                    rho = corr(residual[:,1:],residual[:,0])
                    delta = np.arctanh(np.clip(rho[0],-.999999,.999999))-np.arctanh(np.clip(rho[1],-.999999,.999999))
                    common = dict(study=study,patient=patient,scoring=scoring,model=model,n_cells=len(g),n_samples=g['sample'].nunique())
                    patient_results.append(common | dict(rho_claudin=rho[0],rho_polarity=rho[1],delta_z=delta))
                    diagnostics.append(common | dict(design_columns=design.shape[1],design_rank=rank,residual_df=len(g)-rank,
                                                      finite_comparison=bool(np.isfinite(delta)),covariates=';'.join(covariates)))
                    for gene in members:
                        r = np.nan
                        if gene in gene_positions:
                            values = expr[idx,gene_positions[gene]]
                            if np.ptp(values) > 0:
                                rr = residual_rank(values[:,None],design)
                                r = corr(rr,residual[:,0])[0]
                        gene_results.append(common | dict(gene=gene,rho=r,z=np.arctanh(np.clip(r,-.999999,.999999))))
                    # Independent SVD projection validates the primary least-squares implementation.
                    if scoring=='fixed_background' and model=='common':
                        u,s,_ = np.linalg.svd(design,full_matrices=False)
                        q = u[:,s > s[0]*max(design.shape)*np.finfo(float).eps]
                        rr = stats.rankdata(g[['HRC','Claudin','Polarity']].to_numpy(),axis=0)
                        rr -= q @ (q.T @ rr)
                        check = np.corrcoef(rr,rowvar=False)[0,1:]
                        error = float(np.nanmax(np.abs(rho-check)))
                        assert error < 1e-9
                        numerical_checks.append(common | dict(max_rho_error=error))
        print('Scored',obs.study.unique().tolist(),len(obs),'cells',flush=True)
    patients = pd.DataFrame(patient_results); genes = pd.DataFrame(gene_results)
    write(pd.concat(all_scores,ignore_index=True),'cell_scores.tsv.gz')
    write(patients,'patient_component_effects.tsv'); write(genes,'patient_member_correlations.tsv')
    write(pd.DataFrame(detection),'patient_member_detection.tsv'); write(pd.DataFrame(diagnostics),'model_diagnostics.tsv')
    write(pd.DataFrame(numerical_checks),'independent_projection_checks.tsv')
    summaries = []
    for keys,g in patients.groupby(['study','scoring','model']):
        base = dict(zip(['study','scoring','model'],keys))
        for endpoint,column in [('Claudin_minus_Polarity','delta_z'),('Claudin','rho_claudin'),('Polarity','rho_polarity')]:
            x = g[column].to_numpy() if endpoint=='Claudin_minus_Polarity' else np.arctanh(np.clip(g[column].to_numpy(),-.999999,.999999))
            summaries.append(base|dict(endpoint=endpoint)|summary(x,str(keys)+endpoint))
    member_summaries = []
    for keys,g in genes.groupby(['study','scoring','model','gene']):
        member_summaries.append(dict(zip(['study','scoring','model','gene'],keys))|summary(g.z,str(keys)))
    write(pd.DataFrame(summaries),'study_component_effects.tsv')
    write(pd.DataFrame(member_summaries),'study_member_effects.tsv')
    (OUT/'run_info.json').write_text(json.dumps(lock_record|{'external_cells':len(aobs),'ogden_cells':len(oobs),
        'external_raw_max_integer_error':integer_error,'atlas_path':str(ATLAS.relative_to(ROOT)),
        'atlas_bytes':ATLAS.stat().st_size,'ogden_path':str(OGDEN.relative_to(ROOT)),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n',encoding='utf-8')
    print(pd.DataFrame(summaries).query("scoring=='fixed_background' and endpoint=='Claudin_minus_Polarity'").to_string(index=False),flush=True)


if __name__=='__main__':
    main()
