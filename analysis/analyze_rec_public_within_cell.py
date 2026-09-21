"""Patient-level conditional repair/junction coexistence with disjoint programs."""
import json
import hashlib
import itertools
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse, stats
from rec_public_upgrade_common import ROOT, META, OUT, sha256

dest=OUT / "within_cell"; dest.mkdir(parents=True, exist_ok=True)
sets={k:set(v.gene) for k,v in pd.read_csv(META / "analysis_gene_sets.tsv",sep="\t").groupby("set_id")}
source=ROOT / "data_sources/Ogden_2025_CRLM_multiome/CRCLM_multiome_GEX_decontaminated.h5ad"
a=ad.read_h5ad(source,backed="r")
obs=a.obs.loc[a.obs.Cell_type.astype(str).eq("Epithelial"),["Patient","Cell_subtype"]].copy()
idx=a.obs.index.get_indexer(obs.index)
x=sparse.csr_matrix(a.X[idx,:],dtype=np.float64)
genes=pd.Index(a.var_names.str.upper());assert genes.is_unique
a.file.close()
lib=np.asarray(x.sum(axis=1)).ravel();ng=np.asarray(x.getnnz(axis=1)).ravel()
log=x.multiply((10000/lib)[:,None]).tocsr();np.log1p(log.data,out=log.data)
means=np.asarray(log.mean(axis=0)).ravel();det=np.asarray((x>0).sum(axis=0)).ravel()
pool=np.flatnonzero(det>=5)
bins=np.full(len(genes),-1,dtype=int)
bins[pool]=pd.qcut(pd.Series(means[pool]).rank(method="first"),25,labels=False).to_numpy()
positions={g:i for i,g in enumerate(genes)}
groups=obs.groupby(["Patient","Cell_subtype"],observed=True).indices
tech=np.column_stack([np.log1p(lib),np.log1p(ng)])

def residual(v,cov):
    cov=np.asarray(cov);valid=np.std(cov,axis=0)>0
    m=np.column_stack([np.ones(len(v)),cov[:,valid]])
    return v-m@np.linalg.lstsq(m,v,rcond=None)[0]

def score(g,controls,method):
    p=[positions[z] for z in sorted(g)]
    if method=="matched":
        return np.asarray(log[:,p].mean(axis=1)).ravel()-np.asarray(log[:,controls].mean(axis=1)).ravel()
    t=log[:,p];mu=np.asarray(t.mean(axis=0)).ravel();sd=np.sqrt(np.maximum(np.asarray(t.power(2).mean(axis=0)).ravel()-mu**2,0));keep=sd>0
    return np.asarray(t[:,keep].multiply(1/sd[keep]).mean(axis=1)).ravel()-np.mean(mu[keep]/sd[keep])

locks=[];rows=[]
for rival in ("WHITE_RSC","ICMS3_MINUS_ICMS2"):
    components=["WHITE_RSC"] if rival=="WHITE_RSC" else ["ICMS3_TEMPLATE_UP","ICMS2_TEMPLATE_UP"]
    rivalgenes=set.union(*(sets[k] for k in components))
    overlap=sets["CANELLAS_CORE_HRC"] & sets["REACTOME_TIGHT_JUNCTION_INTERACTIONS"]
    targets={"repair":sets["CANELLAS_CORE_HRC"]-overlap-rivalgenes,"junction":sets["REACTOME_TIGHT_JUNCTION_INTERACTIONS"]-overlap-rivalgenes}
    targets.update({k:sets[k] for k in components})
    targets={k:{g for g in gs if g in positions and det[positions[g]]>=5} for k,gs in targets.items()}
    assert not (targets["repair"] & targets["junction"])
    all_targets=set.union(*targets.values())
    forbidden={positions[g] for g in all_targets};used=set();controls={}
    rng=np.random.default_rng(42+int.from_bytes(hashlib.sha256(rival.encode()).digest()[:4],"little"))
    for k,gs in targets.items():
        if len(gs)<3:raise ValueError(f"Insufficient disjoint genes for {rival}/{k}")
        ctrl=[]
        for b in sorted({bins[positions[g]] for g in gs}):
            choices=np.array([i for i in pool[bins[pool]==b] if i not in forbidden and i not in used])
            take=rng.choice(choices,size=min(50,len(choices)),replace=False)
            ctrl.extend(take.tolist());used.update(take.tolist())
        controls[k]=ctrl
        for g in sorted(gs):locks.append(dict(rival=rival,component=k,role="target",gene=g))
        for i in ctrl:locks.append(dict(rival=rival,component=k,role="control",gene=genes[i]))
    # Gene/control membership is fixed before evaluating within-patient relations.
    pd.DataFrame(locks).to_csv(dest / "conditional_gene_control_lock.tsv",sep="\t",index=False)
    for method in ("matched","gene_z"):
        scores={k:score(gs,controls[k],method) for k,gs in targets.items()}
        rival_score=scores[components[0]] if len(components)==1 else scores[components[0]]-scores[components[1]]
        for (patient,state),ii in groups.items():
            if state not in ("REC","Hypoxia","UPR","iREC") or len(ii)<10:continue
            for adjust in ("technical_only","technical_plus_rival"):
                cov=tech[ii] if adjust=="technical_only" else np.column_stack([tech[ii],rival_score[ii]])
                xx=residual(scores["repair"][ii],cov);yy=residual(scores["junction"][ii],cov)
                rho=float(stats.spearmanr(xx,yy).statistic) if np.ptp(xx)>0 and np.ptp(yy)>0 else np.nan
                rows.append(dict(rival=rival,method=method,adjustment=adjust,patient=patient,state=state,n_cells=len(ii),n_repair_genes=len(targets["repair"]),n_junction_genes=len(targets["junction"]),rho=rho))
patient=pd.DataFrame(rows);patient.to_csv(dest / "patient_conditional_correlations.tsv",sep="\t",index=False)
summ=[]
for threshold in (10,20,30):
 for key,g in patient.loc[patient.n_cells>=threshold].groupby(["rival","method","adjustment","state"]):
    vals=g.rho.dropna().to_numpy();zs=np.arctanh(np.clip(vals,-.999999,.999999));n=len(vals)
    rng=np.random.default_rng(42);b=np.tanh(rng.choice(zs,(10000,n)).mean(axis=1));lo,hi=np.quantile(b,[.025,.975])
    signs=np.array(list(itertools.product([-1,1],repeat=n)));p=np.mean(np.abs(signs@zs/n)>=abs(zs.mean())-1e-12)
    summ.append(dict(zip(["rival","method","adjustment","state"],key))|dict(minimum_cells=threshold,n_patients=n,n_positive=int(sum(vals>0)),rho=float(np.tanh(zs.mean())),CI_low=lo,CI_high=hi,p_signflip=p))
summary=pd.DataFrame(summ)
summary["p_BH"]=summary.groupby(["minimum_cells","method","adjustment"])["p_signflip"].transform(lambda s:stats.false_discovery_control(s))
summary.to_csv(dest / "conditional_summary.tsv",sep="\t",index=False)
(dest / "methods.json").write_text(json.dumps(dict(source=str(source),cells=len(obs),patients=obs.Patient.nunique(),normalization="log1p(10000 * corrected counts/library)",targets="Repair and junction genes exclude mutual overlap and ALL conditioning-rival genes",controls="25 expression bins, up to 50 genes per represented bin; excludes all targets and previously assigned controls",inference="Equal patient Fisher-z summary, 10000 patient bootstrap draws; exhaustive patient sign flips; primary min20 cells, sensitivity10/30",conditional_covariates="Within patient/state: log1p library size, detected genes, and rival score",gene_lock_sha256=sha256(dest / "conditional_gene_control_lock.tsv")),indent=2)+"\n",encoding="utf-8")
print(summary.query("minimum_cells==20 and method=='matched' and state=='REC'").to_string(index=False),flush=True)
