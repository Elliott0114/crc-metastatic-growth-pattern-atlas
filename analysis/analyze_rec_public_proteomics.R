source("scripts/rec_public_upgrade_common.R")
dest<-file.path(OUT,"proteomics");dir.create(dest,recursive=TRUE,showWarnings=FALSE)
p<-file.path(ROOT,"data_sources/Nissen_2025_CRLM_proteomics/mmc4.xlsx")
d<-as.data.frame(read_excel(p,sheet="1 Sample Overview",na=c("","NA","NaN")))
a<-as.data.table(read_excel(p,sheet="5 Total Cohort Abundances",skip=2,na=c("","NA","NaN")))
samples<-d$Sample;stopifnot(setequal(samples,names(a)[-(1:2)]))
a<-a[!is.na(`Gene name`)&`Gene name`!=""]
a[,gene:=toupper(`Gene name`)]
a<-a[,lapply(.SD,function(x)if(all(is.na(x)))NA_real_ else mean(as.numeric(x),na.rm=TRUE)),by=gene,.SDcols=samples]
x<-as.matrix(a[,..samples]);rownames(x)<-a$gene;storage.mode(x)<-"double"
d$patient_id<-as.character(d$Patient_Nr);d$center<-as.character(d$Batch)
stopifnot(!anyDuplicated(d$Sample))
pat<-unique(d[,c("patient_id","center")]);stopifnot(!anyDuplicated(pat$patient_id))
xp<-sapply(pat$patient_id,function(id)apply(x[,d$patient_id==id,drop=FALSE],1,function(v)if(all(is.na(v)))NA_real_ else median(v,na.rm=TRUE)))
rownames(xp)<-rownames(x);colnames(xp)<-pat$patient_id
z<-gene_z(xp);sets<-read_sets()
rec<-sets$FROZEN_REC_TOP50
protein_genes<-intersect(rec,rownames(z));stopifnot(length(protein_genes)==18)
save_tsv(data.frame(set_id="REC_MEASURABLE_PROTEIN_18",gene=protein_genes),file.path(META,"rec_protein18_genes.tsv"))
coverage<-data.frame(gene=rec,measured_in_integrated_matrix=rec%in%rownames(z),n_patients_detected=vapply(rec,function(g)if(g%in%rownames(xp))sum(is.finite(xp[g,]))else 0L,integer(1)),n_samples_detected=vapply(rec,function(g)if(g%in%rownames(x))sum(is.finite(x[g,]))else 0L,integer(1)))
save_tsv(coverage,file.path(dest,"rec50_protein_coverage.tsv"))
for(line in readLines(file.path(OLD,"msigdb_c2.cp.reactome_2025.1.Hs.gmt"))) {
 f<-strsplit(line,"\t",fixed=TRUE)[[1]]
 if(f[1]=="REACTOME_TIGHT_JUNCTION_INTERACTIONS")sets[[f[1]]]<-f[-c(1,2)]
}
targets<-c("CANELLAS_CORE_HRC","REACTOME_TIGHT_JUNCTION_INTERACTIONS","HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT")
all_scores<-pat;all_scores$REC18<-score_set(z,protein_genes)
all_scores$n_REC_proteins<-colSums(is.finite(xp[protein_genes,,drop=FALSE]))
for(id in targets)all_scores[[id]]<-score_set(z,sets[[id]])
save_tsv(all_scores,file.path(dest,"patient_program_scores.tsv"))
gene_rows<-list();rows<-list();pair_data<-list()
partial_rank<-function(a) {
  if(nrow(a)<6||sd(a$x)==0||sd(a$y)==0)return(c(rho=NA_real_,p=NA_real_))
  a$rx<-rank(a$x);a$ry<-rank(a$y)
  k<-length(unique(a$center))-1L
  if(k>0) {rx<-resid(lm(rx~factor(center),a));ry<-resid(lm(ry~factor(center),a))}else{rx<-a$rx;ry<-a$ry}
  r<-cor(rx,ry);df<-nrow(a)-k-2
  c(rho=r,p=2*pt(-abs(r*sqrt(df/(1-r*r))),df=df))
}
for(id in targets) {
  partner<-intersect(sets[[id]],rownames(z));overlap<-intersect(protein_genes,partner)
  left<-setdiff(protein_genes,overlap);right<-setdiff(partner,overlap)
  gene_rows[[length(gene_rows)+1]]<-data.frame(partner=id,side="REC18",gene=left)
  gene_rows[[length(gene_rows)+1]]<-data.frame(partner=id,side="partner",gene=right)
  q<-data.frame(pat,x=score_set(z,left),y=score_set(z,right),complete_proteins=colSums(!is.finite(z[c(left,right),,drop=FALSE]))==0)
  q<-q[complete.cases(q[,c("x","y")]),,drop=FALSE]
  pair_data[[length(pair_data)+1]]<-data.frame(partner=id,q)
  specs<-c("pooled_center_adjusted",paste0("center_",sort(unique(pat$center))),paste0("omit_",sort(unique(pat$center))),"complete_proteins_only")
  for(spec in specs) {
    b<-q
    if(startsWith(spec,"center_"))b<-b[b$center==sub("center_","",spec),,drop=FALSE]
    if(startsWith(spec,"omit_"))b<-b[b$center!=sub("omit_","",spec),,drop=FALSE]
    if(spec=="complete_proteins_only")b<-b[b$complete_proteins,,drop=FALSE]
    result<-partial_rank(b);ci<-c(NA_real_,NA_real_);B<-2000L;nvalid<-0
    if(is.finite(result["rho"])) {
      groups<-split(seq_len(nrow(b)),b$center)
      draws<-replicate(B,{idx<-unlist(lapply(groups,function(ix)sample(ix,length(ix),replace=TRUE)),use.names=FALSE);partial_rank(b[idx,,drop=FALSE])["rho"]})
      nvalid<-sum(is.finite(draws));if(nvalid>=.95*B)ci<-quantile(draws,c(.025,.975),na.rm=TRUE,names=FALSE)
    }
    rows[[length(rows)+1]]<-data.frame(partner=id,specification=spec,n_patients=nrow(b),n_centers=length(unique(b$center)),n_REC_disjoint=length(left),n_partner_disjoint=length(right),n_shared_removed=length(overlap),rho=result["rho"],p=result["p"],CI_low=ci[1],CI_high=ci[2],bootstrap_valid=nvalid,bootstrap_B=B)
  }
}
models<-rbindlist(rows);models[,p_BH:=p.adjust(p,"BH"),by=specification]
save_tsv(models,file.path(dest,"disjoint_protein_correlations.tsv"));save_tsv(rbindlist(gene_rows),file.path(dest,"disjoint_protein_gene_sets.tsv"));save_tsv(rbindlist(pair_data),file.path(dest,"disjoint_pair_patient_scores.tsv"))
save_tsv(d,file.path(dest,"source_sample_manifest.tsv"))

# The deposited table supplies diagnosis/referral-based survival but no specimen
# date or diagnostic-to-resection delay. No time-aligned survival model is fit.
surv<-d[!is.na(d$`survival status`)&d$`survival status`%in%c("alive","dead")&is.finite(as.numeric(d$`survival time (Years)`))&as.numeric(d$`survival time (Years)`)>0,,drop=FALSE]
surv$event<-as.integer(surv$`survival status`=="dead")
surv$survival_model_status<-"not_estimable_specimen_time_origin_alignment_unavailable"
save_tsv(surv,file.path(dest,"survival_candidate_samples.tsv"))
u<-surv[!duplicated(surv$patient_id),,drop=FALSE]
save_tsv(data.frame(status="not_estimable",n_patients_with_outcome=nrow(u),n_deaths=sum(u$event),reason="Survival starts at diagnosis/referral of individual LM; specimen dates/delays are absent, recurrence and lesion ordinal do not recover elapsed time; repeated lesions must not duplicate patient outcomes."),file.path(dest,"survival_feasibility.tsv"))
writeLines(c("Source: Nissen et al., 2025, doi:10.1016/j.mcpro.2025.101026; Supplemental Table S5.","Author preprocessing: log2, column-median normalization, HarmonizR/ComBat across centres, proteins detected in at least 70% of samples. Missing values remain; no new imputation.","Non-survival estimand: patient median across available lesions; patient-level gene z scores; at least 80% of each measurable module required per patient, plus complete-protein sensitivity.","Correlations remove all shared proteins, adjust ranks for centre, and bootstrap patients within centre (2000 replicates). Approximate partial-rank P values; BH within each specification across four module pairs.","Protein module scores do not validate HGP, assembled junctions, or a full 50-gene signature."),file.path(dest,"methods.txt"),useBytes=TRUE)
capture.output(sessionInfo(),file=file.path(dest,"sessionInfo.txt"))
print(models[specification=="pooled_center_adjusted"]);print(data.frame(n_patients=nrow(pat),REC_proteins=length(protein_genes),survival_candidates=nrow(u),deaths=sum(u$event)))
