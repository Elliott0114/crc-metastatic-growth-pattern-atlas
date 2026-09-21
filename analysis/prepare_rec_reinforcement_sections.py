"""Build one outcome-blind same-section tumour pair per patient from source counts."""
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

from rec_junction_context_common import ROOT, canonical, digest

OUT = ROOT / 'analysis_results/rec_manuscript_reinforcement_2026-09-17/sections'
SOURCE = ROOT / 'data_sources/Liu_2026_GSE294385'


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    qc = pd.read_csv(ROOT/'analysis_results/gse294385_rec_program_extension/sample_qc.tsv', sep='\t')
    paired = qc[(qc.micro_tumour_spots > 0) & (qc.macro_tumour_spots > 0)].copy()
    selected = paired.sort_values(['patient','micro_tumour_spots','macro_tumour_spots','sample'], ascending=[True,False,False,True]).drop_duplicates('patient')
    paired['selected'] = paired['sample'].isin(selected['sample'])
    paired.to_csv(OUT/'section_selection.tsv',sep='\t',index=False,encoding='utf-8')
    assert len(paired)==12 and len(selected)==10
    old = pd.read_csv(ROOT/'analysis_results/rec_junction_context_2026-09-16/measurement_qc/regional_common_counts.tsv.gz',sep='\t',index_col=0)
    genes = pd.Index([canonical(x) for x in old.index]); assert not genes.has_duplicates
    annotation = pd.read_csv(SOURCE/'visium_liver_meta_after_qc.tsv.gz',sep='\t')
    outputs, metadata, hashes = [], [], []
    for row in selected.itertuples(index=False):
        path = SOURCE/'extracted'/row.sample/'filtered_feature_bc_matrix'
        names = pd.read_csv(path/'features.tsv.gz',sep='\t',header=None,usecols=[1])[1].map(canonical)
        barcodes = pd.read_csv(path/'barcodes.tsv.gz',sep='\t',header=None)[0]
        with gzip.open(path/'matrix.mtx.gz','rb') as f:
            mat = sparse.csc_matrix(mmread(f))
        assert mat.shape == (len(names),len(barcodes))
        labels = annotation[annotation['sample']==row.sample].drop_duplicates('spot_barcode').set_index('spot_barcode').reindex(barcodes)
        for kind in ['micro','macro']:
            layer = 'Liver '+('micrometastasis' if kind=='micro' else 'macrometastasis')+' tumor'
            ix = np.flatnonzero(labels.Layer3.eq(layer).to_numpy())
            expected = getattr(row,kind+'_tumour_spots'); assert len(ix)==expected,(row.sample,kind,len(ix),expected)
            total = np.asarray(mat[:,ix].sum(axis=1)).ravel()
            vec = pd.Series(total,index=names).groupby(level=0).sum().reindex(genes)
            assert not vec.isna().any() and (vec>=0).all()
            identity = str(row.patient)+'|'+kind+'_tumour'
            outputs.append(vec.rename(identity))
            metadata.append({'id':identity,'patient':row.patient,'sample':row.sample,'region':kind+'_tumour',
                             'macro':int(kind=='macro'),'n_spots':len(ix),'library_size':int(total.sum())})
        hashes.append({'sample':row.sample,'matrix_sha256':digest(path/'matrix.mtx.gz')})
        print('Prepared',row.patient,row.sample,flush=True)
    counts=pd.concat(outputs,axis=1); counts.index.name='gene'
    counts.to_csv(OUT/'counts.tsv.gz',sep='\t',encoding='utf-8')
    pd.DataFrame(metadata).to_csv(OUT/'metadata.tsv',sep='\t',index=False,encoding='utf-8')
    pd.DataFrame(hashes).to_csv(OUT/'source_hashes.tsv',sep='\t',index=False,encoding='utf-8')
    print('Prepared',counts.shape,'count matrix',flush=True)


if __name__=='__main__':
    main()
