source("scripts/rec_public_upgrade_common.R")
suppressPackageStartupMessages(library(limma))
dest<-file.path(OUT,"mapk");sets<-read_sets()
sets$FOCUSED_JUNCTION_5<-c("CLDN3","CLDN4","CLDN7","CRB3","F11R")
programs<-c("FROZEN_REC_TOP50","WHITE_RSC","WHITE_CBC","WHITE_MAPK","WHITE_WNT","CANELLAS_CORE_HRC","REACTOME_TIGHT_JUNCTION_INTERACTIONS","FOCUSED_JUNCTION_5","HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT")
x<-read_matrix(file.path(dest,"human_ortholog_pseudobulk_counts.tsv.gz"));meta<-as.data.frame(fread(file.path(dest,"pseudobulk_metadata.tsv")))
stopifnot(all(x>=0),all(x==round(x)))
results<-list();cover<-list();all_scores<-list();gene_result<-list()
for(spec in c("epithelial","strict_epithelial","all_QC")) {
  a<-meta[meta$selection==spec&meta$eligible,,drop=FALSE]
  a$treatment<-factor(a$treatment,levels=c("Vehicle","MRTX1133"))
  stopifnot(!anyDuplicated(a$mouse));n<-table(a$treatment)
  if(any(n<3))stop("Fewer than three eligible biological units in a MAPK arm")
  y<-x[,a$sample_id,drop=FALSE];design<-model.matrix(~treatment,a)
  d<-DGEList(y,group=a$treatment);keep<-filterByExpr(d,design=design);d<-calcNormFactors(d[keep,,keep.lib.sizes=FALSE])
  eff<-d$samples$lib.size*d$samples$norm.factors
  lc<-cpm(y,lib.size=eff,log=TRUE,prior.count=1);z<-gene_z(lc)
  v<-voom(d,design,plot=FALSE);fit<-eBayes(lmFit(v,design));tt<-topTable(fit,coef=2,number=Inf,sort.by="none")
  save_tsv(data.frame(gene=rownames(tt),tt),file.path(dest,paste0(spec,"_limma_voom_genes.tsv")))
  score<-a
  for(id in programs) {
    score[[id]]<-score_set(z,sets[[id]])
    cover[[length(cover)+1]]<-data.frame(selection=spec,program=id,n_defined=length(sets[[id]]),n_measured=sum(sets[[id]]%in%rownames(z)))
  }
  all_scores[[spec]]<-score
  isdrug<-as.integer(a$treatment=="MRTX1133")
  for(id in programs) {
    val<-score[[id]];obs<-mean(val[isdrug==1])-mean(val[isdrug==0]);groups<-split(seq_along(val),isdrug)
    draws<-replicate(5000,{i<-sample(groups[["1"]],length(groups[["1"]]),TRUE);j<-sample(groups[["0"]],length(groups[["0"]]),TRUE);mean(val[i])-mean(val[j])})
    ci<-quantile(draws,c(.025,.975),names=FALSE)
    alloc<-combn(seq_along(val),sum(isdrug));null<-apply(alloc,2,function(ix)mean(val[ix])-mean(val[-ix]))
    p<-mean(abs(null)>=abs(obs)-1e-12)
    results[[length(results)+1]]<-data.frame(selection=spec,program=id,n_drug=sum(isdrug),n_vehicle=sum(isdrug==0),effect=obs,CI_low=ci[1],CI_high=ci[2],p_permutation=p,n_allocations=ncol(alloc),bootstrap_B=5000)
  }
  for(g in intersect(c("LGR5","DUSP5","DUSP6","ETV4","ETV5","EPCAM","CLDN3","CLDN4","CLDN7","F11R","GATA6","HNF4A","MKI67"),rownames(tt)))gene_result[[length(gene_result)+1]]<-data.frame(selection=spec,gene=g,tt[g,,drop=FALSE])
}
result<-rbindlist(results);result[,p_BH:=p.adjust(p_permutation,"BH"),by=selection]
save_tsv(result,file.path(dest,"program_intervention_effects.tsv"));save_tsv(rbindlist(cover),file.path(dest,"program_coverage.tsv"));save_tsv(rbindlist(all_scores),file.path(dest,"sample_program_scores.tsv"));save_tsv(rbindlist(gene_result),file.path(dest,"key_gene_intervention_effects.tsv"))
capture.output(sessionInfo(),file=file.path(dest,"sessionInfo.txt"))
print(result[selection=="epithelial"])
