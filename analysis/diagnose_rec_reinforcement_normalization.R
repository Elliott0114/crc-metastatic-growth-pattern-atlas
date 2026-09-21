# Descriptive normalization diagnostics; never chooses models by significance.
suppressPackageStartupMessages({library(data.table);library(edgeR)})
base<-'analysis_results/rec_evidence_bridge_2026-09-16'
out<-'analysis_results/rec_manuscript_reinforcement_2026-09-17/diagnostics'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
readgz<-function(p)fread(cmd=paste('gzip -dc',shQuote(p)))
canonical<-function(g){g<-toupper(g);g[g=='MPP5']<-'PALS1';g}
fits<-readRDS(file.path(base,'counts/primary_fits.rds'))
bulk<-readgz('analysis_results/rec_junction_context_2026-09-16/measurement_qc/gse151165_raw_counts.tsv.gz')
bulk_ids<-fits[['Bulk_rHGP_minus_dHGP']]$meta$id
bulkmat<-as.matrix(bulk[,..bulk_ids]);rownames(bulkmat)<-canonical(bulk[[1]])
sp<-readgz(file.path(base,'spatial/region_counts.tsv.gz'));spmat<-as.matrix(sp[,-1]);rownames(spmat)<-sp[[1]]
def<-fread('analysis_results/rec_manuscript_reinforcement_2026-09-17/bulk/canonical_gene_sets.tsv')
markers<-list(Epithelial=def[set_id=='PAN_EPITHELIAL_7',gene],Liver=def[set_id=='HEPATOCYTE_CONTEXT_6',gene])
patient<-list();ma<-list();corr<-list();summary<-list()
for(ctx in c('Bulk_rHGP_minus_dHGP','Spatial_rHGP_minus_dHGP_tumour','Spatial_rHGP_minus_dHGP_near')){
 fit<-fits[[ctx]];m<-copy(fit$meta);m<-m[match(colnames(fit$d),id)]
 raw<-if(startsWith(ctx,'Bulk'))bulkmat[,m$id,drop=FALSE]else spmat[,m$id,drop=FALSE]
 stopifnot(all(abs(colSums(raw)-m$library_size)<1e-6))
 nf<-fit$d$samples$norm.factors
 m[,`:=`(dataset=ctx,norm_factor=nf,log2_norm_factor=log2(nf),
  top1_count_fraction=apply(raw,2,max)/library_size,
  top10_count_fraction=apply(raw,2,function(v)sum(sort(v,decreasing=TRUE)[1:10]))/library_size,
  top_gene=rownames(raw)[max.col(t(raw),ties.method='first')])]
 lc<-cpm(raw,lib.size=m$library_size,log=TRUE,prior.count=.5)
 for(k in names(markers)){
  g<-intersect(markers[[k]],rownames(raw));m[,(paste0(k,'_mean_logCPM')):=colMeans(lc[g,,drop=FALSE])]
 }
 for(v in c('top10_count_fraction','Epithelial_mean_logCPM','Liver_mean_logCPM')){
  corr[[length(corr)+1]]<-data.table(dataset=ctx,variable=v,n_profiles=nrow(m),
   spearman_with_log2_norm_factor=cor(m[[v]],m$log2_norm_factor,method='spearman'))
 }
 summary[[ctx]]<-data.table(dataset=ctx,n_profiles=nrow(m),n_genes_eligible=nrow(fit$d),
   min_factor=min(nf),max_factor=max(nf),r_minus_d_log2_factor=mean(log2(nf)[m$r==1])-mean(log2(nf)[m$r==0]),
   min_top10_fraction=min(m$top10_count_fraction),max_top10_fraction=max(m$top10_count_fraction))
 for(norm in c('TMM','library_only')){
  cp<-cpm(raw,lib.size=m$library_size*if(norm=='TMM')nf else 1)
  rr<-rowMeans(cp[,m$r==1,drop=FALSE]);dd<-rowMeans(cp[,m$r==0,drop=FALSE]);ok<-rr>0 & dd>0
  ma[[paste(ctx,norm)]]<-data.table(dataset=ctx,normalization=norm,gene=rownames(raw)[ok],
    A=(log2(rr[ok])+log2(dd[ok]))/2,M=log2(rr[ok]/dd[ok]),count_eligible=rownames(raw)[ok]%in%rownames(fit$d))
 }
 patient[[ctx]]<-m
}
fwrite(rbindlist(patient,fill=TRUE),file.path(out,'patient_normalization_diagnostics.tsv'),sep='\t')
fwrite(rbindlist(ma),file.path(out,'MA_gene_summaries.tsv.gz'),sep='\t')
fwrite(rbindlist(corr),file.path(out,'descriptive_proxy_correlations.tsv'),sep='\t')
fwrite(rbindlist(summary),file.path(out,'normalization_summary.tsv'),sep='\t')
print(rbindlist(summary));print(rbindlist(corr))
