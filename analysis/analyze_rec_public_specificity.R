source("scripts/rec_public_upgrade_common.R")
dest<-file.path(OUT,"specificity");dir.create(dest,recursive=TRUE,showWarnings=FALSE)
sets<-read_sets();rec<-sets$FROZEN_REC_TOP50
newsets<-c("WHITE_RSC","WHITE_CBC","WHITE_MAPK","WHITE_WNT","ICMS2_TEMPLATE_UP","ICMS3_TEMPLATE_UP")
protein18<-fread(file.path(META,"rec_protein18_genes.tsv"))$gene
sets$REC_PROTEIN18<-protein18
signed_score<-function(z,id,remove=character()) {
  if(id=="ICMS3_MINUS_ICMS2")return(score_set(z,setdiff(sets$ICMS3_TEMPLATE_UP,remove))-score_set(z,setdiff(sets$ICMS2_TEMPLATE_UP,remove)))
  score_set(z,setdiff(sets[[id]],remove))
}
extra<-c(newsets,"ICMS3_MINUS_ICMS2","FROZEN_REC_TOP50","CANELLAS_CORE_HRC","REACTOME_TIGHT_JUNCTION_INTERACTIONS","HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","REC_PROTEIN18")
gs_genes<-function(id)if(id=="ICMS3_MINUS_ICMS2")union(sets$ICMS3_TEMPLATE_UP,sets$ICMS2_TEMPLATE_UP)else sets[[id]]

raw<-as.data.table(read_excel(file.path(ROOT,"data_sources/GSE151165/GSE151165_RNA_seq.raw_read_count.xlsx")))
sample_cols<-grep("^KR-",names(raw),value=TRUE);setnames(raw,1,"gene");raw[,gene:=toupper(trimws(gene))]
raw<-raw[!is.na(gene)&gene!="",lapply(.SD,function(x)sum(as.numeric(x))),by=gene,.SDcols=sample_cols]
x<-as.matrix(raw[,..sample_cols]);rownames(x)<-raw$gene
mt<-data.frame(sample_id=sample_cols,hgp=c(rep("dHGP",18),rep("rHGP",12)),tissue=rep(c("liver","tumour","liver","tumour"),c(9,9,6,6)))
mt<-mt[mt$tissue=="tumour",];x<-x[,mt$sample_id]
dge<-DGEList(x,group=mt$hgp);keep<-filterByExpr(dge,group=mt$hgp);dge<-calcNormFactors(dge[keep,,keep.lib.sizes=FALSE]);eff<-dge$samples$lib.size*dge$samples$norm.factors
logcpm<-cpm(x,lib.size=eff,log=TRUE,prior.count=1);z<-gene_z(logcpm)
save_tsv(data.frame(gene=rownames(logcpm),logcpm,check.names=FALSE),file.path(dest,"gse151165_logcpm.tsv"))
scores<-mt
for(id in extra)scores[[id]]<-signed_score(z,id)
scores$r<-as.numeric(scores$hgp=="rHGP")
save_tsv(scores,file.path(dest,"gse151165_program_scores.tsv"))
coverage<-rbindlist(lapply(extra,function(id)data.table(dataset="GSE151165",set_id=id,n_defined=length(gs_genes(id)),n_measured=sum(gs_genes(id)%in%rownames(z)))))
baseline<-mean(scores$FROZEN_REC_TOP50[scores$r==1])-mean(scores$FROZEN_REC_TOP50[scores$r==0])
stopifnot(abs(baseline-0.505162955238302)<1e-8)
rows<-list()
for(id in c("WHITE_RSC","ICMS3_MINUS_ICMS2","WHITE_CBC"))for(deoverlap in c(FALSE,TRUE)) {
  remove<-if(deoverlap)intersect(rec,gs_genes(id))else character()
  a<-data.frame(y=score_set(z,setdiff(rec,remove)),r=scores$r,c=signed_score(z,id,remove))
  fit<-lm(y~r+c,a);est<-coef(fit)["r"]
  groups<-split(seq_len(nrow(a)),a$r)
  draws<-replicate(5000,{ix<-unlist(lapply(groups,function(g)sample(g,length(g),replace=TRUE)));coef(lm(y~r+c,a[ix,]))["r"]})
  valid<-sum(is.finite(draws));ci<-if(valid>=4750)quantile(draws,c(.025,.975),na.rm=TRUE,names=FALSE)else c(NA_real_,NA_real_)
  unc<-mean(a$y[a$r==1])-mean(a$y[a$r==0])
  rows[[length(rows)+1]]<-data.frame(rival=id,deoverlap=deoverlap,n_patients=nrow(a),n_rHGP=sum(a$r),n_dHGP=sum(a$r==0),n_shared_removed=length(remove),n_REC_measured=sum(setdiff(rec,remove)%in%rownames(z)),n_rival_measured=sum(setdiff(gs_genes(id),remove)%in%rownames(z)),unadjusted_effect=unc,adjusted_effect=est,CI_low=ci[1],CI_high=ci[2],retention=est/unc,REC_rival_spearman=cor(a$y,a$c,method="spearman"),design_VIF=1/(1-cor(a$r,a$c)^2),bootstrap_valid=valid,bootstrap_B=5000,CI_method="Within-HGP patient bootstrap; no adjusted exact permutation P")
}
models<-rbindlist(rows);save_tsv(models,file.path(dest,"new_rival_HGP_models.tsv"))

# Existing bootstrap estimates remain valid; remove unsupported adjusted-label
# permutation P values and the FDR values derived from them in the new version.
pb<-read_matrix(file.path(OLD,"ogden_direct_state_pseudobulk_counts.tsv.gz"))
pm<-as.data.frame(fread(file.path(OLD,"ogden_direct_state_pseudobulk_metadata.tsv")))
pm<-pm[pm$n_cells>=20,];pb<-pb[,pm$sample_id,drop=FALSE]
pdge<-calcNormFactors(DGEList(pb));pz<-gene_z(cpm(pdge,log=TRUE,prior.count=1))
ps<-pm
for(id in extra)ps[[id]]<-signed_score(pz,id)
save_tsv(ps,file.path(dest,"ogden_patient_state_program_scores.tsv"))
coverage<-rbind(coverage,rbindlist(lapply(extra,function(id)data.table(dataset="Ogden_patient_state",set_id=id,n_defined=length(gs_genes(id)),n_measured=sum(gs_genes(id)%in%rownames(pz))))))
contrasts<-list()
for(id in c(newsets,"ICMS3_MINUS_ICMS2"))for(state in c("Hypoxia","UPR","iREC")) {
  a<-ps[ps$state=="REC",c("patient",id)];b<-ps[ps$state==state,c("patient",id)]
  ab<-merge(a,b,by="patient",suffixes=c("_rec","_other"));v<-ab[[2]]-ab[[3]]
  n<-length(v);ci<-c(NA_real_,NA_real_);p<-NA_real_
  if(n>=3) {
    ci<-quantile(replicate(10000,mean(sample(v,n,replace=TRUE))),c(.025,.975),names=FALSE)
    signs<-as.matrix(expand.grid(rep(list(c(-1,1)),n)));p<-mean(abs(as.numeric(signs%*%v)/n)>=abs(mean(v))-1e-12)
  }
  contrasts[[length(contrasts)+1]]<-data.frame(program=id,comparator=state,n_pairs=n,effect=mean(v),CI_low=ci[1],CI_high=ci[2],n_positive=sum(v>0),p_signflip=p)
}
cc<-rbindlist(contrasts);cc[,p_BH:=p.adjust(p_signflip,"BH")]
save_tsv(cc,file.path(dest,"ogden_program_paired_contrasts.tsv"))

concordance<-list()
for(dataset in c("GSE151165","Ogden_patient_state","GSE159216")) {
  if(dataset=="GSE151165"){zz<-z;unit<-scores$sample_id}
  if(dataset=="Ogden_patient_state"){zz<-pz;unit<-ps$sample_id}
  if(dataset=="GSE159216"){
    cm<-as.data.frame(fread(file.path(OUT,"clinical/sample_patient_manifest.tsv")));cx<-read_matrix(file.path(OUT,"clinical/gene_expression_log2.tsv.gz"));ix<-which(cm$primary_eligible);zz<-gene_z(cx[,ix]);unit<-cm$patient_id[ix]
  }
  a<-score_set(zz,rec);b<-score_set(zz,protein18)
  concordance[[dataset]]<-data.frame(dataset,n_units=length(a),n_REC50_measured=sum(rec%in%rownames(zz)),n_protein18_measured=sum(protein18%in%rownames(zz)),spearman_rho=cor(a,b,method="spearman"),interpretation=if(dataset=="Ogden_patient_state")"Descriptive, correlated patient-state profiles"else "Descriptive patient-level projection agreement")
  save_tsv(data.frame(unit,REC50=a,REC18=b),file.path(dest,paste0(dataset,"_RNA18_concordance_scores.tsv")))
}
save_tsv(rbindlist(concordance),file.path(dest,"RNA18_REC50_concordance.tsv"));save_tsv(coverage,file.path(dest,"program_coverage.tsv"))
capture.output(sessionInfo(),file=file.path(dest,"sessionInfo.txt"))
print(models);print(cc[program%in%c("WHITE_RSC","ICMS3_MINUS_ICMS2")]);print(rbindlist(concordance))
