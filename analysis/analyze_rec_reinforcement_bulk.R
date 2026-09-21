# Canonical identifiers; unchanged bulk projection normalization and explicit estimands.
suppressPackageStartupMessages({library(data.table);library(limma)})
out<-'analysis_results/rec_manuscript_reinforcement_2026-09-17/bulk'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
can<-function(x){x<-toupper(x);x[x=='MPP5']<-'PALS1';x}
raw<-fread('analysis_results/rec_public_upgrade_2026-09-07/specificity/gse151165_logcpm.tsv')
source_gene<-raw[[1]];x<-as.matrix(raw[,-1]);rownames(x)<-can(source_gene)
stopifnot(!anyDuplicated(rownames(x)))
meta<-fread('analysis_results/rec_public_upgrade_2026-09-07/specificity/gse151165_program_scores.tsv')[,.(sample_id,hgp,r)]
x<-x[,meta$sample_id,drop=FALSE]
valid<-apply(x,1,sd)>0 & rowSums(!is.finite(x))==0
x<-x[valid,,drop=FALSE];z<-t(scale(t(x)))
definition<-fread('metadata/rec_public_upgrade_2026-09-07/analysis_gene_sets.tsv')
definition[,source_gene:=gene];definition[,gene:=can(gene)]
sets<-lapply(split(definition$gene,definition$set_id),unique)
ids<-c(HRC='CANELLAS_CORE_HRC',Junction='REACTOME_TIGHT_JUNCTION_INTERACTIONS',E2F='HALLMARK_E2F_TARGETS',G2M='HALLMARK_G2M_CHECKPOINT',REC50='FROZEN_REC_TOP50',Epithelial='PAN_EPITHELIAL_7',Liver='HEPATOCYTE_CONTEXT_6')
target<-sets[ids];names(target)<-names(ids)
overlap<-intersect(target$HRC,target$Junction)
target$HRC_disjoint<-setdiff(target$HRC,overlap);target$Junction_disjoint<-setdiff(target$Junction,overlap)
score<-function(g,mat=z)colMeans(mat[intersect(g,rownames(mat)),,drop=FALSE])
scores<-copy(meta)
coverage<-list()
for(k in names(target)){
 scores[,(k):=score(target[[k]])]
 coverage[[k]]<-data.table(component=k,gene=target[[k]],measured_variable=target[[k]]%in%rownames(z))
}
stopifnot(sum(coverage$Junction$measured_variable)==28,abs(mean(scores$REC50[scores$r==1])-mean(scores$REC50[scores$r==0])-.505162955238302)<1e-8)
fwrite(scores,file.path(out,'patient_scores.tsv'),sep='\t');fwrite(rbindlist(coverage),file.path(out,'gene_coverage.tsv'),sep='\t')
fwrite(data.table(gene=rownames(x),x,check.names=FALSE),file.path(out,'logcpm.tsv.gz'),sep='\t')
fwrite(definition,file.path(out,'canonical_gene_sets.tsv'),sep='\t')
alloc<-combn(15,6);groups<-split(seq_len(15),meta$r)
set.seed(20260917)
boot<-replicate(10000,unlist(lapply(groups,function(g)sample(g,length(g),replace=TRUE))))
effect<-function(y){
 est<-mean(y[meta$r==1])-mean(y[meta$r==0])
 nm<-colMeans(matrix(y[alloc],nrow=6));null<-nm-(sum(y)-6*nm)/9
 by<-matrix(y[boot],nrow=15);br<-matrix(meta$r[boot],nrow=15)
 draws<-colSums(by*br)/6-colSums(by*(1-br))/9
 ci<-quantile(draws,c(.025,.975),names=FALSE)
 list(effect=est,low=ci[1],high=ci[2],p=mean(abs(null)>=abs(est)-1e-12))
}
rows<-list();loo<-list()
for(k in c('HRC','Junction','E2F','G2M','HRC_disjoint','Junction_disjoint')){
 rows[[k]]<-as.data.table(c(list(component=k,n_patients=15L,n_genes=sum(coverage[[k]]$measured_variable)),effect(scores[[k]])))
 for(i in seq_len(15))for(method in c('fixed_scores','gene_restandardization')){
  y<-if(method=='fixed_scores')scores[[k]][-i] else score(target[[k]],t(scale(t(x[,-i,drop=FALSE]))))
  rr<-meta$r[-i]
  loo[[length(loo)+1]]<-data.table(component=k,omitted=meta$sample_id[i],method=method,effect=mean(y[rr==1])-mean(y[rr==0]))
 }
}
effects<-rbindlist(rows);effects[,family:=ifelse(grepl('disjoint',component),'two_disjoint','four_full')]
effects[,q:=p.adjust(p,'BH'),by=family]
fwrite(effects,file.path(out,'score_effects.tsv'),sep='\t');fwrite(rbindlist(loo),file.path(out,'patient_omissions.tsv'),sep='\t')
models<-list();model_omissions<-list();model_diagnostics<-list()
for(k in c('HRC','Junction'))for(covariate in c('Epithelial','Liver')){
 shared<-intersect(target[[k]],target[[covariate]])
 for(disjoint in if(length(shared))c(FALSE,TRUE)else FALSE){
  removed<-if(disjoint)shared else character()
  y<-score(setdiff(target[[k]],removed));cv<-score(setdiff(target[[covariate]],removed))
  design<-cbind(1,r=meta$r,context=cv);est<-qr.solve(design,y)[2]
  draws<-apply(boot,2,function(ix){a<-design[ix,,drop=FALSE];if(qr(a)$rank==3)qr.solve(a,y[ix])[2] else NA_real_})
  ok<-is.finite(draws);ci<-if(sum(ok)>=9500)quantile(draws[ok],c(.025,.975),names=FALSE)else c(NA_real_,NA_real_)
  models[[length(models)+1]]<-data.table(component=k,covariate=covariate,deoverlap=disjoint,shared_removed=paste(removed,collapse=';'),
     n_genes=length(intersect(setdiff(target[[k]],removed),rownames(z))),effect=est,low=ci[1],high=ci[2],valid_bootstraps=sum(ok))
  fitted_model<-lm(y~meta$r+cv)
  model_diagnostics[[length(model_diagnostics)+1]]<-data.table(component=k,covariate=covariate,deoverlap=disjoint,
    sample_id=meta$sample_id,residual=residuals(fitted_model),cook_distance=cooks.distance(fitted_model),
    design_rank=qr(design)$rank,hgp_VIF=1/(1-cor(meta$r,cv)^2))
  for(i in seq_len(15))model_omissions[[length(model_omissions)+1]]<-data.table(component=k,covariate=covariate,deoverlap=disjoint,
    omitted=meta$sample_id[i],effect=qr.solve(design[-i,,drop=FALSE],y[-i])[2])
 }
}
fwrite(rbindlist(models),file.path(out,'composition_models.tsv'),sep='\t')
fwrite(rbindlist(model_omissions),file.path(out,'composition_patient_omissions.tsv'),sep='\t')
fwrite(rbindlist(model_diagnostics),file.path(out,'composition_diagnostics.tsv'),sep='\t')
camera_rows<-list();directions<-list();design<-model.matrix(~r,as.data.frame(meta))
for(disjoint in c(FALSE,TRUE)){
 nameset<-if(disjoint)c('HRC_disjoint','Junction_disjoint')else c('HRC','Junction')
 ind<-ids2indices(target[nameset],rownames(x))
 res<-as.data.table(camera(x,ind,design,contrast=2,inter.gene.cor=NA,trend.var=TRUE,sort=FALSE),keep.rownames='component')
 res[,`:=`(disjoint=disjoint,n_background=nrow(x))];camera_rows[[length(camera_rows)+1]]<-res
 for(k in nameset){
  g<-intersect(target[[k]],rownames(z));dz<-rowMeans(z[g,meta$r==1,drop=FALSE])-rowMeans(z[g,meta$r==0,drop=FALSE])
  directions[[k]]<-data.table(component=k,gene=g,gene_z_effect=dz)
 }
}
fwrite(rbindlist(camera_rows),file.path(out,'camera_tests.tsv'),sep='\t')
fwrite(rbindlist(directions),file.path(out,'gene_directions.tsv'),sep='\t')
fwrite(data.table(gene=rownames(x)),file.path(out,'camera_background.tsv'),sep='\t')
capture.output(sessionInfo(),file=file.path(out,'session_info.txt'))
print(effects);print(rbindlist(camera_rows));print(rbindlist(models))
