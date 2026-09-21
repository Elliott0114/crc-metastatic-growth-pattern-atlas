# Analysis: fixed HGP component contrasts and composition sensitivities
# Date: 2026-09-07; seed: 42; R/package versions: emitted sessionInfo.txt
source("scripts/rec_public_upgrade_common.R")
sets <- read_sets()
previous <- OUT
dest <- file.path(ROOT, "analysis_results/rec_hgp_focused_revision_2026-09-07/hgp")
dir.create(dest, recursive=TRUE, showWarnings=FALSE)
set.seed(42)
B <- 10000L
s <- as.data.frame(fread(file.path(previous, "specificity/gse151165_program_scores.tsv")))
x <- read_matrix(file.path(previous, "specificity/gse151165_logcpm.tsv"))
x <- x[, s$sample_id, drop=FALSE]
z <- gene_z(x)
programs <- c("FROZEN_REC_TOP50", "CANELLAS_CORE_HRC", "REACTOME_TIGHT_JUNCTION_INTERACTIONS",
              "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT")
stopifnot(nrow(s)==15L, !anyDuplicated(s$sample_id), sum(s$r)==6L)
stopifnot(all(vapply(programs, function(id) max(abs(score_set(z, sets[[id]])-s[[id]]))<1e-12, logical(1))))
groups <- split(seq_len(nrow(s)), s$r)
alloc <- combn(nrow(s), sum(s$r))
exact <- function(y) {
  est <- mean(y[s$r==1])-mean(y[s$r==0])
  null <- colMeans(matrix(y[alloc], nrow=nrow(alloc)))
  null <- null-(sum(y)-nrow(alloc)*null)/(length(y)-nrow(alloc))
  mean(abs(null)>=abs(est)-1e-12)
}
boot_indices <- function() unlist(lapply(groups, function(g) sample(g, length(g), replace=TRUE)))
rows <- list(); coverage <- list(); loo <- list()
for (j in seq_along(programs)) {
  id <- programs[j]; y <- s[[id]]; set.seed(42+j)
  draws <- replicate(B, {ii<-boot_indices(); mean(y[ii][s$r[ii]==1])-mean(y[ii][s$r[ii]==0])})
  ci <- quantile(draws, c(.025,.975), names=FALSE)
  rows[[j]] <- data.frame(program=id, role=if(j==1) "existing_primary_REC" else "extension_secondary_component",
    n_patients=nrow(s), n_rHGP=sum(s$r), n_dHGP=sum(s$r==0), effect=mean(y[s$r==1])-mean(y[s$r==0]),
    CI_low=ci[1], CI_high=ci[2], bootstrap_B=B, bootstrap_seed=42+j, p_allocation=exact(y),
    allocations=ncol(alloc), missing_scores=sum(!is.finite(y)))
  measured <- intersect(sets[[id]], rownames(z))
  coverage[[j]] <- data.frame(program=id, n_defined=length(sets[[id]]), n_measured=length(measured),
    measured_genes=paste(measured, collapse=";"), missing_or_constant_genes=paste(setdiff(sets[[id]], measured),collapse=";"))
  for (i in seq_len(nrow(s))) {
    keep <- seq_len(nrow(s))!=i
    for (method in c("fixed_scores", "gene_restandardization")) {
      yy <- if(method=="fixed_scores") y[keep] else score_set(gene_z(x[,keep,drop=FALSE]), sets[[id]])
      rr <- s$r[keep]
      loo[[length(loo)+1L]] <- data.frame(program=id, omitted_patient=s$sample_id[i], method=method,
        effect=mean(yy[rr==1])-mean(yy[rr==0]), n_patients=sum(keep))
    }
  }
}
effects <- rbindlist(rows); effects[, p_BH_four_components:=NA_real_]
effects[role=="extension_secondary_component", p_BH_four_components:=p.adjust(p_allocation, "BH")]
save_tsv(effects, file.path(dest,"component_contrasts.tsv"))
save_tsv(rbindlist(coverage), file.path(dest,"component_coverage.tsv"))
save_tsv(rbindlist(loo), file.path(dest,"component_leave_one_patient.tsv"))
save_tsv(s[,c("sample_id","hgp",programs)], file.path(dest,"patient_component_scores.tsv"))

# Covariate models retain effect/interval reporting; no adjusted label-allocation P.
model_rows <- list(); residual_rows <- list(); counter <- 0L
pdf(file.path(dest,"composition_model_diagnostics.pdf"), width=8, height=7)
for (id in programs[2:3]) for(covariate in c("PAN_EPITHELIAL_7","HEPATOCYTE_CONTEXT_6")) {
  overlap <- intersect(sets[[id]],sets[[covariate]])
  for(deoverlap in if(length(overlap)) c(FALSE,TRUE) else FALSE) {
    counter <- counter+1L
    removed <- if(deoverlap) overlap else character()
    d <- data.frame(y=score_set(z,setdiff(sets[[id]],removed)),r=s$r,
                    context=score_set(z,setdiff(sets[[covariate]],removed)))
    stopifnot(all(is.finite(as.matrix(d))))
    fit <- lm(y~r+context,d); est<-coef(fit)["r"]
    set.seed(4200+counter)
    draws <- replicate(B, {ii<-boot_indices(); f<-lm(y~r+context,d[ii,]); if(f$rank==3L) unname(coef(f)["r"]) else NA_real_})
    valid<-is.finite(draws); ci<-if(sum(valid)>=.95*B) quantile(draws[valid],c(.025,.975),names=FALSE) else c(NA_real_,NA_real_)
    model_rows[[counter]]<-data.frame(program=id,covariate=covariate,deoverlap=deoverlap,n_patients=nrow(d),
      n_rHGP=sum(d$r),n_dHGP=sum(d$r==0),n_program_genes=sum(setdiff(sets[[id]],removed)%in%rownames(z)),
      n_covariate_genes=sum(setdiff(sets[[covariate]],removed)%in%rownames(z)),shared_removed=paste(removed,collapse=";"),
      unadjusted_effect=mean(d$y[d$r==1])-mean(d$y[d$r==0]),effect=est,CI_low=ci[1],CI_high=ci[2],
      model_rank=fit$rank,VIF=1/(1-cor(d$r,d$context)^2),bootstrap_valid=sum(valid),bootstrap_B=B,
      bootstrap_seed=4200+counter,shapiro_residual_p=shapiro.test(resid(fit))$p.value,
      max_cooks_distance=max(cooks.distance(fit)),inference="Within-HGP patient bootstrap; no adjusted permutation P")
    residual_rows[[counter]]<-data.frame(program=id,covariate=covariate,deoverlap=deoverlap,patient=s$sample_id,
      fitted=fitted(fit),residual=resid(fit),leverage=hatvalues(fit),cooks_distance=cooks.distance(fit))
    par(mfrow=c(2,2),oma=c(0,0,2,0));plot(fit,which=c(1,2,3,5),ask=FALSE)
    mtext(paste(id,covariate,"deoverlap",deoverlap),outer=TRUE,cex=.7)
  }
}
dev.off()
save_tsv(rbindlist(model_rows),file.path(dest,"component_composition_models.tsv"))
save_tsv(rbindlist(residual_rows),file.path(dest,"composition_model_patient_diagnostics.tsv"))
capture.output(sessionInfo(),file=file.path(dest,"sessionInfo.txt"))
print(effects[,.(program,effect,CI_low,CI_high,p_allocation,p_BH_four_components)])
print(rbindlist(model_rows)[,.(program,covariate,deoverlap,effect,CI_low,CI_high,bootstrap_valid)])
