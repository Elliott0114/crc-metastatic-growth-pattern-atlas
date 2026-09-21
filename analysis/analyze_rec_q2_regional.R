#!/usr/bin/env Rscript
# Patient-paired count models with disjoint expression proxies for tissue sources.
suppressPackageStartupMessages({library(data.table);library(edgeR);library(limma);library(digest)})
root <- normalizePath('.')
out <- file.path(root,'analysis_results/rec_q2_integration_2026-09-19/regional')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
old <- 'analysis_results/rec_junction_context_2026-09-16'
same <- 'analysis_results/rec_manuscript_reinforcement_2026-09-17/sections'
sets <- fread(file.path(old,'definitions.tsv'))[component %in% c('HRC','Claudin','Polarity','Junction')]
targets <- unique(sets$gene)
proxies <- list(Epithelial=c('CDH1','EPCAM','KRT18','KRT19','KRT8','MUC1','TACSTD2'),
                Liver=c('ALB','APOA1','APOA2','ASGR1','CPS1','TTR'),
                Endothelial=c('VWF','PECAM1','EMCN','KDR','ENG','RAMP2','PLVAP','ESAM'))
proxy_lock <- rbindlist(lapply(names(proxies),function(k)data.table(proxy=k,gene=proxies[[k]])))
proxy_lock[,excluded_target_overlap:=gene %in% targets]
fwrite(proxy_lock,file.path(out,'proxy_definition_lock.tsv'),sep='\t')
proxies <- lapply(proxies,setdiff,y=targets)
stopifnot(all(lengths(proxies)>=2))
read_table <- function(path) if(grepl('\\.gz$',path))fread(cmd=paste('gzip -dc',shQuote(path))) else fread(path)
read_counts <- function(path) {
 x<-read_table(path);g<-toupper(x[[1]]);g[g=='MPP5']<-'PALS1';m<-as.matrix(x[,-1]);rownames(m)<-g
 stopifnot(!anyDuplicated(g),all(is.finite(m)),all(m>=0),all(m==round(m)));m
}
main <- read_counts(file.path(old,'measurement_qc/regional_common_counts.tsv.gz'))
meta <- fread(file.path(old,'measurement_qc/regional_zero_count_baseline.tsv'))
meta[,zero_count_prior_CPM:=NULL]
meta[,`:=`(id=paste(patient,region,sep='|'),macro=as.integer(region=='macro_tumour'))]
meta<-meta[match(colnames(main),id)]
data_list<-list(all_pairs=list(counts=main,meta=meta),same_section=list(counts=read_counts(file.path(same,'counts.tsv.gz')),meta=fread(file.path(same,'metadata.tsv'))))
gene_outputs<-list();programme_outputs<-list();score_outputs<-list();diagnostics<-list();vifs<-list();proxy_corr<-list();coverage<-list();checks<-list();fitted<-list();display<-list()

for(dataset in names(data_list)) {
 counts<-data_list[[dataset]]$counts;m<-copy(data_list[[dataset]]$meta);m<-m[match(colnames(counts),id)]
 stopifnot(identical(colnames(counts),m$id),all(table(m$patient)==2),all(colSums(counts)<=m$library_size))
 keep<-rowSums(cpm(counts,lib.size=m$library_size)>=1)>=ceiling(ncol(counts)/2) & rowSums(counts)>=15
 fwrite(data.table(gene=rownames(counts),retained=keep),file.path(out,paste0(dataset,'_gene_eligibility.tsv')),sep='\t')
 for(norm in c('TMM','library_only')) {
  d<-DGEList(counts,lib.size=m$library_size);d<-d[keep,,keep.lib.sizes=TRUE]
  d<-calcNormFactors(d,method=if(norm=='TMM')'TMM' else 'none')
  # Proxies are defined using assayed genes, independent of outcome-gene eligibility.
  full<-DGEList(counts,lib.size=m$library_size,norm.factors=d$samples$norm.factors)
  expr<-cpm(full,log=TRUE,prior.count=.5)
  zz<-t(scale(t(expr)));mm<-copy(m)
  for(k in names(proxies)) {
   gg<-intersect(proxies[[k]],rownames(expr));gg<-gg[apply(expr[gg,,drop=FALSE],1,sd)>0]
   stopifnot(length(gg)>=2,!any(gg%in%targets))
   mm[,(k):=as.numeric(scale(colMeans(zz[gg,,drop=FALSE])))]
   coverage[[paste(dataset,norm,k)]]<-data.table(dataset=dataset,normalization=norm,proxy=k,gene=proxies[[k]],used=proxies[[k]]%in%gg)
  }
  mm[,`:=`(dataset=dataset,normalization=norm,norm_factor=d$samples$norm.factors)]
  score_outputs[[paste(dataset,norm)]]<-mm
  cc<-as.data.table(as.table(cor(as.matrix(mm[,names(proxies),with=FALSE]))))
  setnames(cc,c('proxy1','proxy2','correlation'));cc[,`:=`(dataset=dataset,normalization=norm)];proxy_corr[[paste(dataset,norm)]]<-cc
  for(model in c('unadjusted','source_adjusted')) {
   form<-if(model=='unadjusted')~factor(patient)+macro else ~factor(patient)+macro+Epithelial+Liver+Endothelial
   design<-model.matrix(form,as.data.frame(mm));co<-match('macro',colnames(design));rank<-qr(design)$rank
   qcid<-paste(dataset,norm,model)
   usable<-rank==ncol(design) & nrow(design)>rank
   diagnostics[[qcid]]<-data.table(dataset=dataset,normalization=norm,model=model,n_patients=uniqueN(m$patient),n_profiles=nrow(m),
     n_genes=sum(keep),design_rank=rank,design_columns=ncol(design),residual_df=nrow(design)-rank,condition_number=kappa(design),estimable=usable)
   if(!usable)next
   if(model=='source_adjusted') {
    bas<-model.matrix(~factor(patient),as.data.frame(mm));vv<-as.matrix(mm[,.(macro,Epithelial,Liver,Endothelial)])
    within<-vv-bas%*%qr.coef(qr(bas),vv)
    for(j in seq_len(ncol(within))) {
     xx<-within[,-j,drop=FALSE];yy<-within[,j]
     fitv<-lm.fit(cbind(1,xx),yy);vif<-sum((yy-mean(yy))^2)/sum(fitv$residuals^2)
     vifs[[paste(qcid,j)]]<-data.table(dataset=dataset,normalization=norm,variable=colnames(vv)[j],within_patient_VIF=vif)
    }
   }
   dd<-estimateDisp(d,design,robust=TRUE)
   fit<-glmQLFit(dd,design,robust=TRUE,prior.count=0)
   qlf<-glmQLFTest(fit,coef=co)
   tab<-as.data.table(qlf$table,keep.rownames='gene')
   tab[,`:=`(dataset=dataset,normalization=norm,model=model,q_genome=p.adjust(PValue,'BH'),q_junction=NA_real_,
              approximate_QL_Wald_low=NA_real_,approximate_QL_Wald_high=NA_real_,approximate_QL_SE=NA_real_)]
   ij<-which(tab$gene %in% sets[component=='Junction',gene]);tab$q_junction[ij]<-p.adjust(tab$PValue[ij],'BH')
   # Model-based information intervals; QL significance is taken from glmQLFTest.
   for(j in ij) {
    mu<-fit$fitted.values[j,];phi<-rep_len(fit$dispersion,nrow(fit))[j]
    w<-mu/(1+phi*mu)
    if(!is.null(fit$weights))w<-w*fit$weights[j,]
    info<-crossprod(design,design*w)
    iv<-tryCatch(solve(info)[co,co],error=function(e)NA_real_)
    qlvar<-rep_len(fit$s2.post,nrow(fit))[j]
    se<-sqrt(iv*qlvar)/log(2)
    df<-rep_len(qlf$df.total,nrow(fit))[j]
    crit<-qt(.975,df=df)
    tab$approximate_QL_SE[j]<-se
    tab$approximate_QL_Wald_low[j]<-tab$logFC[j]-crit*se
    tab$approximate_QL_Wald_high[j]<-tab$logFC[j]+crit*se
   }
   index<-lapply(split(sets$gene,sets$component),function(g)which(rownames(dd)%in%g))
   ft<-as.data.table(fry(dd,index=index,design=design,contrast=co,sort='none'),keep.rownames='component')
   ft[,`:=`(dataset=dataset,normalization=norm,model=model,n_patients=uniqueN(m$patient),
             q_direction=p.adjust(PValue,'BH'),q_mixed=p.adjust(PValue.Mixed,'BH'))]
   ft[,median_gene_logFC:=vapply(component,function(k)median(tab[gene%in%sets[component==k,gene],logFC]),numeric(1))]
   ft[,n_positive_logFC:=vapply(component,function(k)sum(tab[gene%in%sets[component==k,gene],logFC]>0),integer(1))]
   gene_outputs[[qcid]]<-tab;programme_outputs[[qcid]]<-ft
   fitted[[qcid]]<-list(design=design,metadata=mm,coef=co,fit=fit)
   if(model=='unadjusted') {
    old_tab<-if(dataset=='all_pairs')fread('analysis_results/rec_evidence_bridge_2026-09-16/counts/junction_gene_results.tsv')[dataset=='Liu_macro_minus_micro' & normalization==norm] else read_table(file.path(same,'gene_results.tsv.gz'))[normalization==norm]
    z<-merge(tab[,.(gene,logFC,PValue)],old_tab[,.(gene,old_logFC=logFC,old_p=PValue)],by='gene')
    err<-max(abs(z$logFC-z$old_logFC));perr<-max(abs(z$PValue-z$old_p))
    stopifnot(err<1e-7,perr<1e-7)
    checks[[qcid]]<-data.table(dataset=dataset,normalization=norm,n_shared=nrow(z),max_logFC_error=err,max_p_error=perr)
    disp<-copy(mm)
    for(k in names(index)) {
     gg<-sets[component==k,gene];gg<-intersect(gg,rownames(dd))
     disp[,(paste0(k,'_score')):=colMeans(zz[gg,,drop=FALSE])]
    }
    display[[qcid]]<-disp
   }
   print(ft[,.(dataset,normalization,model,component,NGenes,Direction,q_direction,q_mixed)])
  }
 }
}
allgenes<-rbindlist(gene_outputs)
fwrite(allgenes,file.path(out,'genome_results.tsv.gz'),sep='\t',na='NA')
junction<-merge(CJ(dataset=names(data_list),normalization=c('TMM','library_only'),model=c('unadjusted','source_adjusted'),gene=sets[component=='Junction',gene]),allgenes[gene%in%sets[component=='Junction',gene]],by=c('dataset','normalization','model','gene'),all.x=TRUE)
junction[,count_eligible:=!is.na(logFC)]
fwrite(junction,file.path(out,'all30_member_results.tsv'),sep='\t',na='NA')
fwrite(rbindlist(programme_outputs),file.path(out,'programme_tests.tsv'),sep='\t')
fwrite(rbindlist(score_outputs,fill=TRUE),file.path(out,'patient_region_proxy_scores.tsv'),sep='\t')
fwrite(rbindlist(display,fill=TRUE),file.path(out,'patient_region_display_scores.tsv'),sep='\t')
fwrite(rbindlist(diagnostics),file.path(out,'model_diagnostics.tsv'),sep='\t')
fwrite(rbindlist(vifs),file.path(out,'within_patient_VIF.tsv'),sep='\t')
fwrite(rbindlist(proxy_corr),file.path(out,'proxy_correlations.tsv'),sep='\t')
fwrite(rbindlist(coverage),file.path(out,'proxy_coverage.tsv'),sep='\t')
fwrite(rbindlist(checks),file.path(out,'legacy_reconstruction_checks.tsv'),sep='\t')
saveRDS(fitted,file.path(out,'model_fits.rds'))
capture.output(sessionInfo(),file=file.path(out,'session_info.txt'))
cat('Regional analysis complete\n')
