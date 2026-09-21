"""Measurement robustness prompted by opposite gene-level versus aggregate directions."""
import itertools
import subprocess
import numpy as np
import pandas as pd

from rec_junction_context_common import OUT, ROOT, canonical, definitions, mean_summary, score_frame, write
from analyze_rec_junction_contexts import linear_effect


def main():
    sets=definitions()
    counts=pd.read_csv(ROOT/"analysis_results/rec_program_reader_facing_specificity/gse294385_patient_region_all_gene_counts.tsv.gz",sep="\t").pivot(index="gene",columns=["patient","region"],values="raw_count")
    metadata=pd.read_csv(ROOT/"analysis_results/rec_program_reader_facing_specificity/gse294385_patient_region_metadata.tsv",sep="\t").set_index(["patient","region"])
    libraries=metadata.loc[counts.columns,"library_size"]
    manifest=pd.read_csv(ROOT/"metadata/gse294385_liver_sample_manifest.tsv",sep="\t")
    common=None
    for sample in manifest.loc[manifest.selected_for_paired_extension,"sample"]:
        p=ROOT/"data_sources/Liu_2026_GSE294385/extracted"/sample/"filtered_feature_bc_matrix/features.tsv.gz"
        measured=set(pd.read_csv(p,sep="\t",header=None)[1].str.upper())
        common=measured if common is None else common&measured
    counts=counts.loc[counts.index.intersection(sorted(common))]
    counts.index=counts.index.map(canonical)
    assert not counts.index.duplicated().any() and not counts.isna().any().any()
    count_export=counts.copy();count_export.columns=[p+"|"+r for p,r in counts.columns]
    write(count_export.rename_axis("gene").reset_index(),"measurement_qc/regional_common_counts.tsv.gz")
    cpm=counts.div(libraries)*1e6
    qc=pd.DataFrame({"gene":counts.index,"total_count":counts.sum(axis=1).to_numpy(),
        "n_profiles_nonzero":(counts>0).sum(axis=1).to_numpy(),
        "n_profiles_cpm_ge_1":(cpm>=1).sum(axis=1).to_numpy()})
    write(qc[qc.gene.isin(sets["Junction"])],"measurement_qc/regional_junction_abundance.tsv")
    rows=[];coverage=[];all_scores=[]
    transforms={"historical_fixed_count_prior":np.log2((counts+.5).div(libraries+1)*1e6),
                "library_scaled_count_prior":np.log2((counts.add(.5*libraries/libraries.mean(),axis=1)).div(libraries+libraries/libraries.mean())*1e6),
                "fixed_CPM_prior_0.1":np.log2(cpm+.1),"fixed_CPM_prior_1":np.log2(cpm+1)}
    eligible=(cpm>=1).sum(axis=1)>=11
    for transform,expression in transforms.items():
        for filtering,expr in [("all_members",expression),("CPM_ge_1_in_half_profiles",expression.loc[eligible])]:
            scores,cov,z=score_frame(expr,sets,transform+filtering,ddof=0)
            cov["transform"]=transform;cov["filter"]=filtering;coverage.append(cov)
            all_scores.append(scores.reset_index().assign(transform=transform,filter=filtering,score="mean_gene_z"))
            for construction in ["mean_gene_z","mean_log_expression"]:
                if construction=="mean_gene_z":values=scores
                else:values=pd.DataFrame({c:expr.loc[expr.index.intersection(sets[c])].mean(axis=0) for c in ["Claudin","Polarity","Junction","HRC"]})
                delta=values.xs("macro_tumour",level="region")-values.xs("micro_tumour",level="region")
                for component in ["Claudin","Polarity","Junction","HRC"]:
                    effect=mean_summary(delta[component],transform+filtering+construction+component)
                    n=int(cov.loc[cov.component==component,"measured_variable"].sum())
                    rows.append(dict(dataset="GSE294385",transform=transform,filter=filtering,score=construction,
                                     component=component,n_genes=n,**effect))
    write(pd.DataFrame(rows),"measurement_qc/regional_score_sensitivity.tsv")
    write(pd.concat(coverage),"measurement_qc/regional_coverage.tsv")
    write(pd.concat(all_scores),"measurement_qc/regional_scores.tsv")
    prior=pd.DataFrame({"patient":libraries.index.get_level_values(0),"region":libraries.index.get_level_values(1),
                        "library_size":libraries.to_numpy(),"zero_count_prior_CPM":.5/(libraries.to_numpy()+1)*1e6})
    write(prior,"measurement_qc/regional_zero_count_baseline.tsv")
    gene_rows=[]
    for transform,expression in transforms.items():
        z=expression.sub(expression.mean(axis=1),axis=0).div(expression.std(axis=1,ddof=0),axis=0)
        delta=z.xs("macro_tumour",level="region",axis=1)-z.xs("micro_tumour",level="region",axis=1)
        for gene in sets["Junction"]:
            if gene in delta.index:
                gene_rows.append(dict(transform=transform,gene=gene,**mean_summary(delta.loc[gene],transform+gene)))
    write(pd.DataFrame(gene_rows),"measurement_qc/regional_gene_sensitivity.tsv")
    original=pd.read_csv(ROOT/"metadata/rec_public_upgrade_2026-09-07/analysis_gene_sets.tsv",sep="\t")
    legacy=[]
    for transform,expression in transforms.items():
        z=expression.sub(expression.mean(axis=1),axis=0).div(expression.std(axis=1,ddof=0),axis=0)
        for identifier in ["CANELLAS_CORE_HRC","REACTOME_TIGHT_JUNCTION_INTERACTIONS"]:
            # Exact historical source-symbol coverage, before the MPP5/PALS1 resolution.
            gg=[canonical(g) for g in original.loc[original.set_id==identifier,"gene"] if g in common]
            v=z.loc[gg].mean(axis=0)
            delta=v.xs("macro_tumour",level="region")-v.xs("micro_tumour",level="region")
            legacy.append(dict(set_id=identifier,transform=transform,n_genes=len(gg),**mean_summary(delta,"legacy"+identifier+transform)))
    legacy=pd.DataFrame(legacy)
    assert np.isclose(legacy.loc[(legacy.set_id=="REACTOME_TIGHT_JUNCTION_INTERACTIONS")&(legacy["transform"]=="historical_fixed_count_prior"),"effect"].iloc[0],-.278238,atol=1e-6)
    write(legacy,"measurement_qc/exact_historical_set_normalization_sensitivity.tsv")
    spots=pd.read_csv(OUT/"spatial/tumour_spot_scores.tsv.gz",sep="\t")
    means=spots.groupby(["patient","region"])[["Claudin","Polarity","HRC"]].mean()
    delta=means.xs("Liver macrometastasis tumor",level="region")-means.xs("Liver micrometastasis tumor",level="region")
    write(pd.DataFrame([dict(component=c,**mean_summary(delta[c],"per_spot"+c)) for c in delta]),
          "measurement_qc/region_mean_of_normalized_spots.tsv")

    # HGP audit: raw counts give outcome-blind abundance eligibility; retain original logCPM.
    bulk_cache=OUT/"measurement_qc/gse151165_raw_counts.tsv.gz"
    if not bulk_cache.exists():
        # readxl is already in the project environment; do not install openpyxl.
        code='suppressPackageStartupMessages(library(readxl)); suppressPackageStartupMessages(library(data.table)); x <- as.data.table(read_excel("data_sources/GSE151165/GSE151165_RNA_seq.raw_read_count.xlsx", .name_repair="unique_quiet")); fwrite(x,"analysis_results/rec_junction_context_2026-09-16/measurement_qc/gse151165_raw_counts.tsv.gz",sep="\\t")'
        subprocess.run(["Rscript","--vanilla","-e",code],cwd=ROOT,check=True)
    bulk=pd.read_csv(bulk_cache,sep="\t",low_memory=False)
    bulk=bulk.rename(columns={bulk.columns[0]:"gene"});bulk=bulk.dropna(subset=["gene"])
    bulk.gene=bulk.gene.map(lambda x:canonical(str(x).strip()))
    meta=pd.read_csv(OUT/"contexts/hgp_scores.tsv",sep="\t").set_index("sample_id")
    cols=meta.index.tolist()
    bulk=bulk.groupby("gene")[cols].sum()
    cpm=bulk.div(bulk.sum())*1e6
    expr=pd.read_csv(ROOT/"analysis_results/rec_public_upgrade_2026-09-07/specificity/gse151165_logcpm.tsv",sep="\t",index_col=0)
    expr.index=expr.index.map(canonical)
    eligible=(cpm>=1).sum(axis=1)>=8
    bulk_qc=pd.DataFrame({"gene":bulk.index,"total_count":bulk.sum(axis=1).to_numpy(),
                         "n_profiles_cpm_ge_1":(cpm>=1).sum(axis=1).to_numpy()})
    write(bulk_qc[bulk_qc.gene.isin(sets["Junction"])],"measurement_qc/hgp_junction_abundance.tsv")
    results=[];coverage=[]
    for filtering,e in [("all_members",expr),("CPM_ge_1_in_half_profiles",expr.loc[expr.index.intersection(bulk.index[eligible])])]:
        values,cov,_=score_frame(e,sets,filtering,ddof=1);coverage.append(cov.assign(filter=filtering))
        for construction in ["mean_gene_z","mean_log_expression"]:
            if construction=="mean_log_expression":
                values=pd.DataFrame({c:e.loc[e.index.intersection(sets[c])].mean(axis=0) for c in ["Claudin","Polarity","Junction","HRC"]})
            for component in ["Claudin","Polarity","Junction","HRC"]:
                v=values.loc[cols,component].to_numpy();r=meta.r.to_numpy().astype(bool)
                effect=v[r].mean()-v[~r].mean();null=[]
                for ids in itertools.combinations(range(15),6):
                    select=np.zeros(15,dtype=bool);select[list(ids)]=True
                    null.append(v[select].mean()-v[~select].mean())
                interval=linear_effect(pd.DataFrame({"value":v,"r":r.astype(float)}),"value","r",[],"hgp_qc"+filtering+construction+component,True)
                results.append(dict(component=component,filter=filtering,score=construction,effect=effect,
                    low=interval["low"],high=interval["high"],
                    p=np.mean(np.abs(null)>=abs(effect)-1e-12),n_genes=int(cov.loc[cov.component==component,"measured_variable"].sum())))
    write(pd.DataFrame(results),"measurement_qc/hgp_score_sensitivity.tsv")
    write(pd.concat(coverage),"measurement_qc/hgp_coverage.tsv")
    # The five-patient epithelial source check also aggregates unequal cell counts.
    sc=pd.read_csv(OUT/"sensitivities/hgp_epithelial_pseudobulk_counts.tsv.gz",sep="\t",index_col=0)
    smeta=pd.read_csv(OUT/"sensitivities/hgp_epithelial_pseudobulk_metadata.tsv",sep="\t").set_index("patient")
    lib=smeta.loc[sc.columns,"library_size"]
    positive=sc.sum(axis=1)>0;sc=sc.loc[positive]
    transforms_sc={"historical_fixed_count_prior":np.log2((sc+.5).div(lib+1)*1e6),
                   "library_scaled_count_prior":np.log2(sc.add(.5*lib/lib.mean(),axis=1).div(lib+lib/lib.mean())*1e6),
                   "fixed_CPM_prior_1":np.log2(sc.div(lib)*1e6+1)}
    sc_rows=[]
    for transform,e in transforms_sc.items():
        values,cov,_=score_frame(e,sets,transform,ddof=1)
        for component in ["Claudin","Polarity"]:
            v=values[component].to_numpy();r=smeta.loc[values.index,"hgp"].eq("rHGP").to_numpy()
            effect=v[r].mean()-v[~r].mean();null=[]
            for ids in itertools.combinations(range(5),3):
                mask=np.zeros(5,dtype=bool);mask[list(ids)]=True;null.append(v[mask].mean()-v[~mask].mean())
            sc_rows.append(dict(transform=transform,component=component,effect=effect,p=np.mean(np.abs(null)>=abs(effect)-1e-12),
                                n_genes=int(cov.loc[cov.component==component,"measured_variable"].sum())))
    write(pd.DataFrame(sc_rows),"measurement_qc/hgp_epithelial_score_sensitivity.tsv")
    print(pd.DataFrame(rows).query('score=="mean_gene_z"')[['component','transform','filter','n_genes','effect','p']].to_string(index=False))
    print(pd.DataFrame(results).to_string(index=False))


if __name__=="__main__":
    main()
