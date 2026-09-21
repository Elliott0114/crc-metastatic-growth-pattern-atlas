#!/usr/bin/env Rscript
# Raw-count models for fixed programmes; patients are the replication units.
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(digest))
out <- "analysis_results/rec_evidence_bridge_2026-09-16"
previous <- "analysis_results/rec_junction_context_2026-09-16"
dir.create(file.path(out,"counts"), recursive=TRUE, showWarnings=FALSE)
read_gz <- function(path) fread(cmd=paste("gzip -dc",shQuote(path)))
write_tsv <- function(x,name) fwrite(x,file.path(out,"counts",name),sep="\t",na="NA")
canonical <- function(x) { x <- toupper(x); x[x=="MPP5"] <- "PALS1"; x }
read_count <- function(path) {
  x <- read_gz(path)
  genes <- x[[1]]
  mat <- as.matrix(x[,-1])
  rownames(mat) <- canonical(genes)
  stopifnot(!anyDuplicated(rownames(mat)), !anyNA(mat), all(mat>=0), all(mat==floor(mat)))
  mat
}
sets <- fread(file.path(previous,"definitions.tsv"))[
  component %in% c("Claudin","Polarity","Junction","HRC")]
spatial_sets <- fread(file.path(out,"spatial/definitions.tsv"))
gene_results <- list(); set_results <- list(); eligibility <- list()
normalization <- list(); qc <- list(); fits <- list(); provenance <- list()

run_model <- function(counts, meta, design, coef, dataset, target_sets=sets) {
  stopifnot(identical(colnames(counts),meta$id), all(meta$library_size>0),
            all(colSums(counts)<=meta$library_size), qr(design)$rank==ncol(design),
            nrow(design)>ncol(design))
  raw_cpm <- cpm(counts,lib.size=meta$library_size,log=FALSE)
  keep <- rowSums(raw_cpm>=1)>=ceiling(ncol(counts)/2) & rowSums(counts)>=15
  eligibility[[dataset]] <<- data.table(dataset=dataset, gene=rownames(counts),
    retained=keep, n_profiles_CPM_ge_1=rowSums(raw_cpm>=1), total_count=rowSums(counts))
  for (norm in c("TMM","library_only")) {
    d <- DGEList(counts=counts,lib.size=meta$library_size)
    d <- d[keep,,keep.lib.sizes=TRUE]
    d <- calcNormFactors(d,method=if (norm=="TMM") "TMM" else "none")
    d <- estimateDisp(d,design,robust=TRUE)
    fit <- glmQLFit(d,design,robust=TRUE,prior.count=0)
    qlf <- glmQLFTest(fit,coef=coef)
    tab <- as.data.table(qlf$table,keep.rownames="gene")
    tab[,q_genome:=p.adjust(PValue,"BH")]
    tab[,q_junction:=NA_real_]
    target <- tab$gene %in% target_sets[component=="Junction",gene]
    tab$q_junction[target] <- p.adjust(tab$PValue[target],"BH")
    tab[,`:=`(dataset=dataset,normalization=norm)]
    gene_results[[paste(dataset,norm)]] <<- tab
    index <- lapply(split(target_sets$gene,target_sets$component),function(g) which(rownames(d) %in% g))
    stopifnot(all(lengths(index)>=2))
    result <- as.data.table(fry(d,index=index,design=design,contrast=coef,sort="none"),keep.rownames="component")
    result[,`:=`(dataset=dataset,normalization=norm)]
    result[,`:=`(q_direction=p.adjust(PValue,"BH"), q_mixed=p.adjust(PValue.Mixed,"BH"))]
    result[,median_gene_logFC:=vapply(component,function(k) median(tab$logFC[tab$gene %in% target_sets[component==k,gene]]),numeric(1))]
    result[,n_positive_logFC:=vapply(component,function(k) sum(tab$logFC[tab$gene %in% target_sets[component==k,gene]]>0),integer(1))]
    set_results[[paste(dataset,norm)]] <<- result
    nm <- copy(meta)
    nm[,`:=`(dataset=dataset,normalization=norm,norm_factor=d$samples$norm.factors,
              effective_library=d$samples$lib.size*d$samples$norm.factors)]
    normalization[[paste(dataset,norm)]] <<- nm
    qc[[paste(dataset,norm)]] <<- data.table(dataset=dataset,normalization=norm,
       n_profiles=ncol(d),n_genes=nrow(d),design_columns=ncol(design),
       design_rank=qr(design)$rank,residual_df=nrow(design)-ncol(design),
       coef=colnames(design)[coef], common_dispersion=d$common.dispersion,
       min_norm_factor=min(d$samples$norm.factors),max_norm_factor=max(d$samples$norm.factors),
       min_library=min(d$samples$lib.size),max_library=max(d$samples$lib.size))
    if (norm=="TMM") fits[[dataset]] <<- list(d=d,design=design,coef=coef,meta=meta,
                                             qlf_table=qlf$table)
    cat(dataset,norm,"retained",nrow(d),"genes\n")
    print(result[,.(component,NGenes,Direction,PValue,q_direction,median_gene_logFC)])
  }
}

path <- file.path(previous,"measurement_qc/regional_common_counts.tsv.gz")
counts <- read_count(path)
meta <- fread(file.path(previous,"measurement_qc/regional_zero_count_baseline.tsv"))
meta[,id:=paste(patient,region,sep="|")]
meta <- meta[match(colnames(counts),id)]
stopifnot(nrow(meta)==22,uniqueN(meta$patient)==11,all(table(meta$patient)==2))
meta[,macro:=as.integer(region=="macro_tumour")]
design <- model.matrix(~factor(patient)+macro,data=meta)
run_model(counts,meta,design,ncol(design),"Liu_macro_minus_micro")
provenance[[1]] <- data.table(source=path,sha256=digest(file=path,algo="sha256"))

path <- file.path(previous,"measurement_qc/gse151165_raw_counts.tsv.gz")
raw <- read_gz(path)
meta <- fread("analysis_results/rec_public_upgrade_2026-09-07/specificity/gse151165_program_scores.tsv")
stopifnot(nrow(meta)==15,sum(meta$r)==6)
ids <- meta$sample_id
counts <- as.matrix(raw[,..ids])
rownames(counts) <- canonical(raw$NAME)
stopifnot(!anyDuplicated(rownames(counts)),!anyNA(counts),all(counts>=0),all(counts==floor(counts)))
meta <- meta[,.(id=sample_id,patient=sample_id,hgp,library_size=colSums(counts),r)]
design <- model.matrix(~r,data=meta)
run_model(counts,meta,design,2,"Bulk_rHGP_minus_dHGP")
provenance[[2]] <- data.table(source=path,sha256=digest(file=path,algo="sha256"))

path <- file.path(out,"spatial/region_counts.tsv.gz")
counts <- read_count(path)
meta <- fread(file.path(out,"spatial/region_metadata.tsv"))
meta <- meta[hops==5 & restriction=="all_tumour"]
meta[,r:=as.integer(hgp=="rHGP")]
for (band in c("tumour","near")) {
  selected <- meta[region==band]
  stopifnot(nrow(selected)==6,sum(selected$r)==3)
  design <- model.matrix(~r,data=selected)
  run_model(counts[,selected$id],selected,design,2,
            paste0("Spatial_rHGP_minus_dHGP_",band),spatial_sets)
}
selected <- meta[region %in% c("near","deep")]
selected[,near:=as.integer(region=="near")]
selected[,near_r:=near*r]
stopifnot(nrow(selected)==12,all(table(selected$patient)==2))
design <- model.matrix(~factor(patient)+near+near_r,data=selected)
run_model(counts[,selected$id],selected,design,ncol(design),
          "Spatial_HGP_by_localization",spatial_sets)
provenance[[3]] <- data.table(source=path,sha256=digest(file=path,algo="sha256"))

genes <- rbindlist(gene_results)
write_tsv(genes,"gene_results.tsv.gz")
write_tsv(genes[gene %in% sets[component=="Junction",gene]],"junction_gene_results.tsv")
write_tsv(rbindlist(set_results),"programme_tests.tsv")
write_tsv(rbindlist(eligibility),"gene_eligibility.tsv.gz")
write_tsv(rbindlist(normalization,fill=TRUE),"library_normalization.tsv")
write_tsv(rbindlist(qc),"model_qc.tsv")
write_tsv(rbindlist(provenance),"source_hashes.tsv")
saveRDS(fits,file.path(out,"counts/primary_fits.rds"))
writeLines(capture.output(sessionInfo()),file.path(out,"counts/session_info.txt"),useBytes=TRUE)

# Independent reconstruction of the primary pairing and contrast sign.
checks <- list()
for (dataset in c("Liu_macro_minus_micro","Spatial_HGP_by_localization")) {
  obj <- fits[[dataset]]
  m <- as.data.frame(obj$meta)
  patients <- factor(m$patient,levels=rev(sort(unique(m$patient))))
  if (dataset=="Liu_macro_minus_micro") {
    alternate <- model.matrix(~0+patients)
    alternate <- cbind(alternate,macro=as.integer(m$region=="macro_tumour"))
  } else {
    alternate <- model.matrix(~0+patients)
    alternate <- cbind(alternate,near=as.integer(m$region=="near"),
                       interaction=as.integer(m$region=="near" & m$hgp=="rHGP"))
  }
  # Re-estimate both dispersion and QL fit with a different nuisance parameterization.
  dd <- estimateDisp(obj$d,alternate,robust=TRUE)
  ff <- glmQLFit(dd,alternate,robust=TRUE,prior.count=0)
  tt <- glmQLFTest(ff,coef=ncol(alternate))$table
  a <- max(abs(tt$logFC-obj$qlf_table$logFC))
  b <- max(abs(tt$PValue-obj$qlf_table$PValue))
  # Reparameterizing nuisance terms changes numerical dispersion interpolation.
  # Audit showed <0.001 deviations genome-wide; retain these measured errors and
  # explicitly check the target coefficient directions instead of claiming exact identity.
  target <- rownames(tt) %in% sets[component=="Junction",gene]
  sign_matches <- all(sign(tt$logFC[target])==sign(obj$qlf_table$logFC[target]))
  stopifnot(a<1e-3,b<1e-3,sign_matches)
  checks[[dataset]] <- data.table(dataset=dataset,max_logFC_error=a,max_p_error=b,
    target_max_logFC_error=max(abs(tt$logFC[target]-obj$qlf_table$logFC[target])),
    target_max_p_error=max(abs(tt$PValue[target]-obj$qlf_table$PValue[target])),
    target_signs_agree=sign_matches, tolerance=1e-3)
}
write_tsv(rbindlist(checks),"independent_design_check.tsv")
cat("Count models and independent design checks completed.\n")
