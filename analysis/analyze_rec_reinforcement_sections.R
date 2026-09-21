# Same-section patient pairs with the established raw-count estimand.
suppressPackageStartupMessages({library(data.table);library(edgeR);library(limma)})
out <- 'analysis_results/rec_manuscript_reinforcement_2026-09-17/sections'
x <- fread(cmd=paste('gzip -dc',shQuote(file.path(out,'counts.tsv.gz')))); genes<-x[[1]]
counts<-as.matrix(x[,-1]);rownames(counts)<-genes
meta<-fread(file.path(out,'metadata.tsv'));meta<-meta[match(colnames(counts),id)]
stopifnot(ncol(counts)==20,uniqueN(meta$patient)==10,all(table(meta$patient)==2),
          identical(colnames(counts),meta$id),all(colSums(counts)<=meta$library_size))
sets<-fread('analysis_results/rec_junction_context_2026-09-16/definitions.tsv')[component %in% c('HRC','Claudin','Polarity','Junction')]
design<-model.matrix(~factor(patient)+macro,as.data.frame(meta))
keep<-rowSums(cpm(counts,lib.size=meta$library_size)>=1)>=10 & rowSums(counts)>=15
fwrite(data.table(gene=genes,retained=keep),file.path(out,'gene_eligibility.tsv'),sep='\t')
programmes<-list();gene_results<-list();norms<-list();fits<-list();checks<-list()
for(norm in c('TMM','library_only')) {
 d<-DGEList(counts,lib.size=meta$library_size);d<-d[keep,,keep.lib.sizes=TRUE]
 d<-calcNormFactors(d,method=if(norm=='TMM')'TMM' else 'none')
 d<-estimateDisp(d,design,robust=TRUE);fit<-glmQLFit(d,design,robust=TRUE,prior.count=0)
 qlf<-glmQLFTest(fit,coef=ncol(design));tab<-as.data.table(qlf$table,keep.rownames='gene')
 tab[,`:=`(normalization=norm,q_genome=p.adjust(PValue,'BH'),q_junction=NA_real_)]
 isj<-tab$gene %in% sets[component=='Junction',gene]
 tab$q_junction[isj]<-p.adjust(tab$PValue[isj],'BH')
 indices<-lapply(split(sets$gene,sets$component),function(g)which(rownames(d)%in%g))
 f<-as.data.table(fry(d,index=indices,design=design,contrast=ncol(design),sort='none'),keep.rownames='component')
 f[,`:=`(normalization=norm,n_patients=10L,n_profiles=20L,q_direction=p.adjust(PValue,'BH'),q_mixed=p.adjust(PValue.Mixed,'BH'))]
 f[,median_gene_logFC:=vapply(component,function(k)median(tab[gene%in%sets[component==k,gene],logFC]),numeric(1))]
 f[,n_positive_logFC:=vapply(component,function(k)sum(tab[gene%in%sets[component==k,gene],logFC]>0),integer(1))]
 old<-fread('analysis_results/rec_evidence_bridge_2026-09-16/counts/junction_gene_results.tsv')[dataset=='Liu_macro_minus_micro' & normalization==norm,.(gene,original_logFC=logFC)]
 shared<-merge(tab[,.(gene,logFC)],old,by='gene')
 checks[[norm]]<-data.table(normalization=norm,shared_junction_genes=nrow(shared),same_direction=sum(sign(shared$logFC)==sign(shared$original_logFC)),
                          spearman_shared_logFC=cor(shared$logFC,shared$original_logFC,method='spearman'),residual_df=nrow(design)-ncol(design))
 fwrite(shared,file.path(out,paste0('shared_gene_comparison_',norm,'.tsv')),sep='\t')
 m<-copy(meta);m[,`:=`(normalization=norm,norm_factor=d$samples$norm.factors,effective_library=d$samples$lib.size*d$samples$norm.factors)]
 norms[[norm]]<-m;programmes[[norm]]<-f;gene_results[[norm]]<-tab;fits[[norm]]<-list(d=d,design=design,meta=meta,fit=fit)
 # Normalization-matched gene-z summaries support patient display, not set-test CI.
 expr<-cpm(d,log=TRUE,prior.count=.5);zz<-t(scale(t(expr)))
 scores<-copy(m)
 for(k in names(indices))scores[,(k):=colMeans(zz[indices[[k]],,drop=FALSE])]
 fwrite(scores,file.path(out,paste0('patient_region_scores_',norm,'.tsv')),sep='\t')
 print(f[,.(component,NGenes,Direction,q_direction,q_mixed,n_positive_logFC)])
}
fwrite(rbindlist(programmes),file.path(out,'programme_tests.tsv'),sep='\t')
fwrite(rbindlist(gene_results),file.path(out,'gene_results.tsv.gz'),sep='\t')
fwrite(rbindlist(norms),file.path(out,'library_normalization.tsv'),sep='\t')
fwrite(rbindlist(checks),file.path(out,'comparison_diagnostics.tsv'),sep='\t')
saveRDS(fits,file.path(out,'fits.rds'))
capture.output(sessionInfo(),file=file.path(out,'session_info.txt'))
