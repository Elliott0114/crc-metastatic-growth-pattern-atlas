"""Broad lineage selection independent of REC/RSC scores in eight mouse libraries."""
import json
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from rec_public_upgrade_common import META, OUT

dest=OUT / "mapk";dest.mkdir(parents=True, exist_ok=True)
manifest=pd.read_csv(META / "mrtx1133_sample_manifest.tsv",sep="\t")
assert manifest.biological_unit.nunique()==8
orth=pd.read_csv(META / "human_mouse_one_to_one.tsv",sep="\t").drop_duplicates(["mouse_ensembl","human_gene"])
mapping=orth.groupby("mouse_ensembl").human_gene.agg(lambda x:next(iter(set(x))) if len(set(x))==1 else None).dropna()
markers=dict(epithelial=["Epcam","Krt8","Krt18","Krt19"],immune=["Ptprc","Lyz2","Cd3d","Cd79a"],stromal=["Col1a1","Col1a2","Dcn","Col3a1"],endothelial=["Pecam1","Kdr","Cdh5","Emcn"])
(dest / "lineage_selection_protocol.json").write_text(json.dumps(dict(markers=markers,QC=">=200 detected genes, >=500 UMI, <=20% mitochondrial UMI; author-called barcodes",epithelial="Mean log-normalized broad epithelial markers exceeds each immune/stromal/endothelial module; at least one epithelial marker detected",strict_sensitivity="At least two epithelial markers detected",all_QC_sensitivity="All QC-passing cells",unit="one source-named mouse/library, minimum20 selected cells for primary inference",normalization="Counts pseudobulk before gene-program scoring; no REC/RSC-based cell gate"),indent=2)+"\n",encoding="utf-8")
tables=[];meta=[];cellrows=[];coverage=[]
for row in manifest.itertuples(index=False):
 with h5py.File(row.path,"r") as h:
    g=h["matrix"];x=sparse.csc_matrix((g["data"][:],g["indices"][:],g["indptr"][:]),shape=tuple(g["shape"][:])).T.tocsr()
    ids=np.array([v.decode() for v in g["features"]["id"][:]])
    genes=np.array([v.decode() for v in g["features"]["name"][:]])
    barcodes=np.array([v.decode() for v in g["barcodes"][:]])
 lib=np.asarray(x.sum(axis=1)).ravel();ng=x.getnnz(axis=1);mito=np.array([g.lower().startswith("mt-") for g in genes]);mf=np.asarray(x[:,mito].sum(axis=1)).ravel()/np.maximum(lib,1)
 q=(ng>=200)&(lib>=500)&(mf<=.20)
 log=x.astype(float).multiply((10000/np.maximum(lib,1))[:,None]).tocsr();np.log1p(log.data,out=log.data)
 ms={}
 for k,gs in markers.items():
    ix=np.flatnonzero(np.isin(genes,gs));assert len(ix)>=3
    ms[k]=np.asarray(log[:,ix].mean(axis=1)).ravel()
 ei=np.flatnonzero(np.isin(genes,markers["epithelial"]));epdet=x[:,ei].getnnz(axis=1)
 epi=q&(epdet>=1)&(ms["epithelial"]>np.maximum.reduce([ms[k] for k in ["immune","stromal","endothelial"]]))
 strict=epi&(epdet>=2)
 cellrows.append(pd.DataFrame(dict(sample_id=row.sample_id,barcode=barcodes,n_counts=lib,n_genes=ng,mitochondrial_fraction=mf,QC=q,epithelial=epi,strict_epithelial=strict,**{k+"_score":v for k,v in ms.items()})))
 human=mapping.reindex(ids).to_numpy();good=pd.notna(human)
 coverage.append(dict(sample_id=row.sample_id,n_features=len(ids),n_one_to_one_features=int(good.sum())))
 for spec,mask in [("epithelial",epi),("strict_epithelial",strict),("all_QC",q)]:
    values=np.asarray(x[mask,:].sum(axis=0)).ravel()
    count=pd.Series(values[good],index=human[good]).groupby(level=0).sum()
    sid=row.sample_id+"__"+spec;count.name=sid;tables.append(count)
    meta.append(dict(sample_id=sid,source_sample=row.sample_id,mouse=row.biological_unit,treatment=row.treatment,selection=spec,n_cells=int(mask.sum()),n_all_cells=x.shape[0],library_size=float(values.sum()),eligible=int(mask.sum())>=20))
 print(row.sample_id,"total",x.shape[0],"QC",int(q.sum()),"epithelial",int(epi.sum()),"strict",int(strict.sum()),flush=True)
pd.concat(tables,axis=1).fillna(0).rename_axis("gene").to_csv(dest / "human_ortholog_pseudobulk_counts.tsv.gz",sep="\t",compression="gzip")
pd.DataFrame(meta).to_csv(dest / "pseudobulk_metadata.tsv",sep="\t",index=False)
pd.concat(cellrows,ignore_index=True).to_csv(dest / "cell_lineage_QC.tsv.gz",sep="\t",index=False,compression="gzip")
pd.DataFrame(coverage).to_csv(dest / "ortholog_feature_coverage.tsv",sep="\t",index=False)
