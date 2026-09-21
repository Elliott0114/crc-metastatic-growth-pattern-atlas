"""Explicitly exploratory checks: component contrasts, gene influence, HGP cell source."""
import itertools

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler

from rec_junction_context_common import (
    OUT, ROOT, bh, canonical, correlations, definitions, mean_summary,
    rank_design, rank_residual, score_frame, write,
)
from analyze_rec_junction_contexts import linear_effect


def main():
    sets = definitions()
    hgp = pd.read_csv(OUT / "contexts/hgp_scores.tsv", sep="\t")
    hgp["difference"] = hgp.Claudin-hgp.Polarity
    fit = linear_effect(hgp,"difference","r",[],"component_contrast_hgp",True)
    y = hgp.difference.to_numpy()
    null=[]
    for ids in itertools.combinations(range(15),6):
        mask=np.zeros(15,dtype=bool);mask[list(ids)]=True
        null.append(y[mask].mean()-y[~mask].mean())
    fit["p"] = np.mean(np.abs(null)>=abs(fit["effect"])-1e-12)
    fit["p_method"] = "exhaustive_5005_patient_allocations"
    delta=pd.read_csv(OUT/"contexts/paired_region_changes.tsv",sep="\t")
    rows=[dict(dataset="GSE151165",contrast="Claudin_minus_Polarity_HGP_effect",**fit),
          dict(dataset="GSE294385",contrast="Claudin_minus_Polarity_regional_effect",
               **mean_summary(delta.Claudin-delta.Polarity,"component_contrast_regions"))]
    direct=pd.DataFrame(rows);direct["q_two_exploratory_contrasts"]=bh(direct.p)
    write(direct,"sensitivities/direct_context_component_contrasts.tsv")

    cell=pd.read_csv(OUT/"cells/cell_scores.tsv.gz",sep="\t")
    locks=pd.read_csv(OUT/"cells/gene_control_lock.tsv",sep="\t")
    path=ROOT/"data_sources/Ogden_2025_CRLM_multiome/CRCLM_multiome_GEX_decontaminated.h5ad"
    a=ad.read_h5ad(path)
    indices=a.obs_names.get_indexer(cell.cell_id)
    assert (indices>=0).all()
    names=pd.Index(list(map(canonical,a.var_names)))
    genes=sorted(set(sets["Claudin"])|set(sets["Polarity"]))
    x=sparse.csr_matrix(a.X)[indices][:,names.get_indexer(genes)].toarray()
    x=np.log1p(x/np.expm1(cell.log_library.to_numpy())[:,None]*1e4)
    influence=[]
    for component in ("Claudin","Polarity"):
        members=locks.loc[(locks.component==component)&(locks.role=="target")&locks.measured,"gene"].tolist()
        for gene in members:
            raw=(len(members)*cell[component+"_raw"].to_numpy()-x[:,genes.index(gene)])/(len(members)-1)
            altered=raw+cell[component].to_numpy()-cell[component+"_raw"].to_numpy()
            patient_rho=[];differences=[]
            for patient,g in cell.groupby("patient"):
                if len(g)<100:continue
                design=rank_design(g,["log_library","log_genes","RSC","iCMS","E2F"],state=True)
                values=np.column_stack([g.HRC,altered[g.index],g.Polarity if component=="Claudin" else g.Claudin])
                res=rank_residual(values,design);rho=correlations(res[:,1:],res[:,0])
                patient_rho.append(rho[0])
                c,p=(rho[0],rho[1]) if component=="Claudin" else (rho[1],rho[0])
                differences.append(np.arctanh(c)-np.arctanh(p))
            influence.append(dict(component=component,removed_gene=gene,n_patients=len(patient_rho),
                altered_component_rho=np.tanh(np.mean(np.arctanh(patient_rho))),
                claudin_minus_polarity_fisher_difference=np.mean(differences)))
    write(pd.DataFrame(influence),"sensitivities/leave_one_gene_out.tsv")
    del a,x

    # Reuse the published reconstruction and its existing >=20-cell source check.
    a=ad.read_h5ad(ROOT/"analysis_results/e_mtab_12022_primary_reconstruction.h5ad")
    mask=a.obs.cell_type.astype(str).eq("CRC/epithelial").to_numpy()
    obs=a.obs.loc[mask].reset_index(drop=True)
    counts=sparse.csr_matrix(a.X)[mask]
    assert np.all(counts.data>=0) and np.allclose(counts.data,np.round(counts.data))
    patient_counts=obs.groupby(["patient","hgp"],observed=True).size().reset_index(name="n_cells")
    eligible=patient_counts[patient_counts.n_cells>=20].copy()
    assert len(eligible)==5 and (eligible.hgp=="rHGP").sum()==3
    raw=[];columns=[]
    for row in eligible.itertuples():
        idx=np.flatnonzero(obs.patient.to_numpy()==row.patient)
        raw.append(np.asarray(counts[idx].sum(axis=0)).ravel());columns.append(row.patient)
    matrix=pd.DataFrame(np.column_stack(raw),index=[canonical(g) for g in a.var_names],columns=columns).groupby(level=0).sum()
    library=matrix.sum()
    write(matrix.rename_axis("gene").reset_index(),"sensitivities/hgp_epithelial_pseudobulk_counts.tsv.gz")
    write(eligible.assign(library_size=eligible.patient.map(library)),"sensitivities/hgp_epithelial_pseudobulk_metadata.tsv")
    measured=matrix.sum(axis=1)>0
    expression=np.log2((matrix.loc[measured]+.5).div(library+1)*1e6)
    scores,coverage,_=score_frame(expression,sets,"E-MTAB-12022_source_sensitivity",ddof=1)
    scores=scores.join(eligible.set_index("patient")[["hgp","n_cells"]]);scores["r"]=(scores.hgp=="rHGP").astype(int)
    write(scores.rename_axis("patient").reset_index(),"sensitivities/hgp_epithelial_scores.tsv")
    write(coverage,"sensitivities/hgp_epithelial_coverage.tsv")
    effects=[]
    for component in ("Claudin","Polarity"):
        result=linear_effect(scores,component,"r",[],"smallhgp"+component,True)
        values=scores[component].to_numpy();null=[]
        for ids in itertools.combinations(range(5),3):
            m=np.zeros(5,dtype=bool);m[list(ids)]=True;null.append(values[m].mean()-values[~m].mean())
        result.update(p=np.mean(np.abs(null)>=abs(result["effect"])-1e-12),p_method="10_patient_allocations_source_sensitivity")
        effects.append(dict(component=component,**result))
    write(pd.DataFrame(effects),"sensitivities/hgp_epithelial_effects.tsv")

    spatial=pd.read_csv(OUT/"spatial/tumour_spot_scores.tsv.gz",sep="\t")
    spatial["macro"]=(spatial.region=="Liver macrometastasis tumor").astype(float)
    baseline=["HRC","Epithelial","Liver","log_library","log_genes","macro","platform_IX"]
    predictive=[]
    for sensitivity,select in [("macro_only",spatial.macro.eq(1)&spatial.eligible),
                               ("at_least_five_neighbours",spatial.n_neighbours>=5)]:
        frame=spatial.loc[select]
        for patient in sorted(frame.patient.unique()):
            train,test=frame[frame.patient!=patient],frame[frame.patient==patient]
            weights=1/train.patient.map(train.patient.value_counts()).to_numpy()
            for component in ("Claudin","Polarity"):
                row=dict(sensitivity=sensitivity,patient=patient,component=component,n_spots=len(test))
                for model,cols in [("baseline",baseline),("context",baseline+["neighbour_liver","neighbour_stroma"])]:
                    cols=[c for c in cols if train[c].std(ddof=0)>1e-12]
                    scaler=StandardScaler().fit(train[cols],sample_weight=weights)
                    fit=LinearRegression().fit(scaler.transform(train[cols]),train[component],sample_weight=weights)
                    pred=fit.predict(scaler.transform(test[cols]))
                    row["rmse_"+model]=np.sqrt(np.mean((pred-test[component])**2))
                row["relative_rmse_improvement"]=1-row["rmse_context"]/row["rmse_baseline"]
                predictive.append(row)
    pred=pd.DataFrame(predictive);write(pred,"sensitivities/spatial_prediction_patients.tsv")
    summary=pd.DataFrame([dict(sensitivity=s,component=c,**mean_summary(g.relative_rmse_improvement,s+c))
                         for (s,c),g in pred.groupby(["sensitivity","component"])])
    write(summary,"sensitivities/spatial_prediction_summary.tsv")
    print(direct.to_string(index=False));print(pd.DataFrame(effects).to_string(index=False));print(summary.to_string(index=False))


if __name__=="__main__":
    main()
