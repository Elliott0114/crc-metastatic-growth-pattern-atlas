source("scripts/rec_public_upgrade_common.R")
dest<-file.path(OUT,"specificity")
old<-fread(file.path(OLD,"rec_external_rival_models.tsv"))
removed<-grep("^(adjusted_exact_p|deoverlap_adjusted_exact_p|adjusted_fdr|deoverlap_adjusted_fdr)",names(old),value=TRUE)
old[,(removed):=NULL]
old[,adjusted_inference:="Patient bootstrap confidence intervals; adjusted label-exchange P values withdrawn"]
save_tsv(old,file.path(dest,"existing_rival_models_corrected.tsv"))
save_tsv(data.frame(removed_column=removed,reason="Exchangeability of raw group labels conditional on biological covariates is not established; raw endpoint permutation tests are retained."),file.path(dest,"statistical_correction_log.tsv"))

