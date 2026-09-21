suppressPackageStartupMessages({library(data.table); library(survival); library(edgeR); library(readxl)})
ROOT <- normalizePath(getwd(), mustWork=TRUE)
VERSION <- "rec_public_upgrade_2026-09-07"
OUT <- file.path(ROOT,"analysis_results",VERSION)
META <- file.path(ROOT,"metadata",VERSION)
DATA <- file.path(ROOT,"data_sources",VERSION)
OLD <- file.path(ROOT,"analysis_results/deep_biology_upgrade_2026-08-31/phase2_mechanistic_specificity")
set.seed(42)
save_tsv <- function(x,path) fwrite(as.data.table(x),path,sep="\t",na="NA",quote="auto")
gene_z <- function(x,reference=seq_len(ncol(x))) {
  mu <- rowMeans(x[,reference,drop=FALSE],na.rm=TRUE)
  sdv <- apply(x[,reference,drop=FALSE],1,sd,na.rm=TRUE)
  keep <- is.finite(sdv) & sdv>0
  sweep(sweep(x[keep,,drop=FALSE],1,mu[keep],"-"),1,sdv[keep],"/")
}
score_set <- function(z,genes,min_genes=3,min_fraction=0.8) {
  genes <- intersect(unique(genes),rownames(z))
  if(length(genes)<min_genes) return(rep(NA_real_,ncol(z)))
  a <- z[genes,,drop=FALSE]; s <- colMeans(a,na.rm=TRUE)
  s[colMeans(is.finite(a))<min_fraction | !is.finite(s)] <- NA_real_
  as.numeric(s)
}
read_sets <- function() {
  p <- file.path(META,"analysis_gene_sets.tsv")
  if(!file.exists(p))p<-file.path(META,"frozen_gene_sets_long.tsv")
  s<-fread(p);split(s$gene,s$set_id)
}
read_matrix <- function(path) {
  d<-if(grepl("\\.gz$",path))fread(cmd=paste("gzip -dc",shQuote(path)))else fread(path)
  x<-as.matrix(d[,-1,with=FALSE]);storage.mode(x)<-"double";rownames(x)<-d[[1]];x
}
safe_cox <- function(formula,d) {
  warnings<-character()
  fit<-tryCatch(withCallingHandlers(coxph(formula,data=d,x=TRUE,model=TRUE,ties="efron",singular.ok=FALSE),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
  if(inherits(fit,"error"))return(list(fit=NULL,status=conditionMessage(fit),warnings=warnings))
  valid<-(!length(coef(fit)) || (all(is.finite(coef(fit))) && all(is.finite(diag(vcov(fit)))))) && !any(grepl("infinite|did not converge|iterations",warnings,ignore.case=TRUE))
  list(fit=fit,status=if(valid)"estimable" else "unstable",warnings=warnings)
}
boot_cox <- function(formula,d,B=1000) {
  vals<-rep(NA_real_,B)
  for(i in seq_len(B)) {
    b<-safe_cox(formula,d[sample.int(nrow(d),replace=TRUE),,drop=FALSE])
    if(b$status=="estimable") vals[i]<-coef(b$fit)["REC"]
  }
  ok<-is.finite(vals)
  ci<-if(sum(ok)>=B*0.95) exp(quantile(vals[ok],c(.025,.975),names=FALSE)) else c(NA_real_,NA_real_)
  c(bootstrap_n_valid=sum(ok),bootstrap_n=B,bootstrap_low=ci[1],bootstrap_high=ci[2])
}
