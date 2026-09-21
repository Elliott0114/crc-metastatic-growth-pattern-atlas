source("scripts/rec_public_upgrade_common.R")
dest<-file.path(OUT,"clinical");dir.create(dest,recursive=TRUE,showWarnings=FALSE)
d<-as.data.frame(fread(file.path(dest,"geo_sample_metadata.tsv"),na.strings=c("","NA","NaN")))
names(d)<-make.names(names(d))
names(d)[names(d)=="patient.number"]<-"patient_id"
d$patient_id<-as.character(d$patient_id)
d$os_time<-as.numeric(d$X60.months.overall.survival..time)
d$os_event<-as.numeric(d$X60.months.overall.survival..status)
d$css_time<-as.numeric(d$X60.months.cancer.specific.survival..time)
d$css_event<-as.numeric(d$X60.months.cancer.specific.survival..status)
d$source_representative<-!is.na(d$os_time)&!is.na(d$os_event)
d$known_repeat<-!is.na(d$repeated.hepatic.resection)&as.character(d$repeated.hepatic.resection)%in%c("2","3")
d$source_excluded<-(!is.na(d$lms.from.initial.subtype.discovery)&d$lms.from.initial.subtype.discovery=="NotIncluded")|(!is.na(d$lms.from.random.forest.prediction.model)&d$lms.from.random.forest.prediction.model=="NotIncluded")
stopifnot(!anyDuplicated(d$patient_id[d$source_representative]))
stopifnot(all(na.omit(d$os_event)%in%c(0,1)),all(na.omit(d$css_event)%in%c(0,1)))
stopifnot(all(na.omit(d$os_time)>=0),all(na.omit(d$os_time)<=60),all(na.omit(d$css_time)<=60))
d$sex<-factor(d$gender,levels=c("Female","Male"))
d$synchronous<-factor(d$synchronous.liver.metastases,levels=c("No","Yes"))
d$margin<-factor(d$r.status.liver,levels=c("R0","R1","R2"))
d$extrahepatic<-factor(d$extra.heaptic.disease,levels=c("No","Yes"))
d$neoadjuvant<-factor(ifelse(is.na(d$chemotherapy.prior.to.tumor.sampling),NA,ifelse(grepl("NeoAdjuvant",d$chemotherapy.prior.to.tumor.sampling),"Yes","No")),levels=c("No","Yes"))
mut<-d$kras.mutation=="Mut" | d$nras.mutation=="Mut"
wild<-d$kras.mutation=="WT" & d$nras.mutation=="WT"
d$RAS<-factor(ifelse(!is.na(mut)&mut,"Mut",ifelse(!is.na(wild)&wild,"WT",NA)),levels=c("WT","Mut"))
d$LMS<-factor(d$lms.from.random.forest.prediction.model,levels=paste0("LMS",1:5))
d$primary_eligible<-d$source_representative & !d$known_repeat & !d$source_excluded
stopifnot(all(d$os_time[d$primary_eligible]>0),all(d$css_time[d$primary_eligible]>0))
d$exclusion_reason<-ifelse(!d$source_representative,"nonrepresentative_sample",ifelse(d$known_repeat,"known_repeat_resection",ifelse(d$source_excluded,"source_quality_exclusion","included")))
x<-read_matrix(file.path(dest,"gene_expression_log2.tsv.gz"));stopifnot(identical(colnames(x),d$sample_id))
sets<-read_sets();rec<-sets$FROZEN_REC_TOP50;stopifnot(length(rec)==50)
z<-gene_z(x,which(d$primary_eligible));d$REC_raw<-score_set(z,rec)
mu<-mean(d$REC_raw[d$primary_eligible]);sdv<-sd(d$REC_raw[d$primary_eligible]);d$REC<-(d$REC_raw-mu)/sdv
save_tsv(data.frame(gene=rec,measured=rec%in%rownames(z)),file.path(dest,"rec50_coverage.tsv"))
save_tsv(data.frame(parameter=c("REC_mean","REC_sd","n_gene_scaling_reference_patients"),value=c(mu,sdv,sum(d$primary_eligible))),file.path(dest,"score_scaling.tsv"))
save_tsv(d,file.path(dest,"sample_patient_manifest.tsv"))
primary<-d[d$primary_eligible,,drop=FALSE]
flow<-data.frame(step=c("deposited_samples","unique_patients","outcome_representatives","after_excluding_known_repeated_resection","after_source_quality_exclusion"),n=c(nrow(d),length(unique(d$patient_id)),sum(d$source_representative),sum(d$source_representative&!d$known_repeat),nrow(primary)))
save_tsv(flow,file.path(dest,"patient_flow.tsv"))
save_tsv(data.frame(field=names(primary),n_missing=vapply(primary,function(x)sum(is.na(x)),integer(1)),n_patients=nrow(primary)),file.path(dest,"clinical_missingness.tsv"))
covars<-c("sex","synchronous","margin","extrahepatic","neoadjuvant","RAS")
tab<-rbindlist(lapply(c(covars,"LMS","os_event","css_event"),function(v){t<-table(primary[[v]],useNA="ifany");data.table(variable=v,level=names(t),n=as.integer(t))}))
save_tsv(tab,file.path(dest,"clinical_characteristics.tsv"))

# Alternative specimen choices are restricted to patients in the index cohort.
initial<-d[!d$known_repeat&!d$source_excluded&d$patient_id%in%primary$patient_id,,drop=FALSE]
initial<-initial[order(initial$sample_id),,drop=FALSE]
alt<-initial[!duplicated(initial$patient_id),c("patient_id","sample_id","REC_raw")]
median_scores<-aggregate(REC_raw~patient_id,initial,median)
specimens<-data.frame(patient_id=primary$patient_id,index_sample=primary$sample_id,index_score=primary$REC_raw,
  alternative_sample=alt$sample_id[match(primary$patient_id,alt$patient_id)],alternative_score=alt$REC_raw[match(primary$patient_id,alt$patient_id)],median_initial_score=median_scores$REC_raw[match(primary$patient_id,median_scores$patient_id)])
save_tsv(specimens,file.path(dest,"sampling_sensitivity_mapping.tsv"))
all_models<-list();ph_rows<-list();coef_rows<-list();curve_rows<-list();model_patients<-list()
for(endpoint in c("OS","CSS")) {
  timevar<-if(endpoint=="OS")"os_time" else "css_time"
  eventvar<-if(endpoint=="OS")"os_event" else "css_event"
  for(sensitivity in c("primary","R0_R1","alternative_initial_sample","median_initial_samples")) {
    a<-primary;a$time<-a[[timevar]];a$event<-a[[eventvar]]
    if(sensitivity=="R0_R1")a<-a[a$margin%in%c("R0","R1"),,drop=FALSE]
    if(sensitivity=="alternative_initial_sample")a$REC<-(specimens$alternative_score-mu)/sdv
    if(sensitivity=="median_initial_samples")a$REC<-(specimens$median_initial_score-mu)/sdv
    a<-droplevels(a[complete.cases(a[,c("time","event","REC",covars,"LMS")]),,drop=FALSE])
    stopifnot(!anyDuplicated(a$patient_id))
    model_patients[[length(model_patients)+1]]<-data.frame(endpoint,sensitivity,patient_id=a$patient_id)
    formulas<-list(univariable=as.formula("Surv(time,event)~REC"),clinical=as.formula(paste("Surv(time,event)~REC+",paste(covars,collapse="+"))),clinical_LMS=as.formula(paste("Surv(time,event)~REC+",paste(c(covars,"LMS"),collapse="+"))))
    for(model in names(formulas)) {
      f<-formulas[[model]];fit<-safe_cox(f,a)
      row<-data.frame(endpoint,sensitivity,model,n_patients=nrow(a),n_events=sum(a$event),status=fit$status,warnings=paste(unique(fit$warnings),collapse=" | "))
      if(!is.null(fit$fit)) {
        b<-coef(fit$fit)["REC"];se<-sqrt(vcov(fit$fit)["REC","REC"])
        row$beta<-b;row$se<-se;row$HR<-exp(b);row$CI_low<-exp(b-1.96*se);row$CI_high<-exp(b+1.96*se);row$p<-2*pnorm(-abs(b/se));row$n_parameters<-length(coef(fit$fit))
        row$events_per_parameter<-sum(a$event)/length(coef(fit$fit))
        reduced<-safe_cox(update(f,.~.-REC),a)
        row$increment_LR_p<-if(!is.null(reduced$fit))pchisq(2*(tail(fit$fit$loglik,1)-tail(reduced$fit$loglik,1)),1,lower.tail=FALSE) else NA_real_
        if(sensitivity=="primary") {
          boot<-boot_cox(f,a);for(n in names(boot))row[[n]]<-boot[[n]]
          co<-summary(fit$fit)$coefficients
          coef_rows[[length(coef_rows)+1]]<-data.frame(endpoint,model,term=rownames(co),co,check.names=FALSE)
          ph<-tryCatch(cox.zph(fit$fit),error=function(e)NULL)
          if(!is.null(ph)) {
            ph_rows[[length(ph_rows)+1]]<-data.frame(endpoint,model,term=rownames(ph$table),ph$table,check.names=FALSE)
            if("REC"%in%rownames(ph$table)&&ph$table["REC","p"]<.05) {
              ft<-tryCatch(coxph(update(f,.~.+tt(REC)),data=a,tt=function(x,t,...)x*log(t/12)),error=function(e)NULL)
              if(!is.null(ft))for(t in c(12,36)) {
                w<-c(1,log(t/12));ix<-c("REC","tt(REC)");be<-sum(coef(ft)[ix]*w);ve<-as.numeric(t(w)%*%vcov(ft)[ix,ix]%*%w)
                curve_rows[[length(curve_rows)+1]]<-data.frame(endpoint,model,type="time_varying_PH_sensitivity",time_months=t,REC_SD=1,HR=exp(be),CI_low=exp(be-1.96*sqrt(ve)),CI_high=exp(be+1.96*sqrt(ve)))
              }
            }
          }
          for(v in seq(quantile(a$REC,.05),quantile(a$REC,.95),length.out=61))curve_rows[[length(curve_rows)+1]]<-data.frame(endpoint,model,type="continuous_HR_relative_to_REC_mean",time_months=NA_real_,REC_SD=v,HR=exp(b*v),CI_low=exp(b*v-1.96*abs(v)*se),CI_high=exp(b*v+1.96*abs(v)*se))
        }
      }
      all_models[[length(all_models)+1]]<-row
    }
  }
}
models<-rbindlist(all_models,fill=TRUE)
models[,p_BH:=p.adjust(p,"BH"),by=.(endpoint,sensitivity)]
save_tsv(models,file.path(dest,"cox_models.tsv"));save_tsv(rbindlist(ph_rows,fill=TRUE),file.path(dest,"proportional_hazards.tsv"));save_tsv(rbindlist(coef_rows,fill=TRUE),file.path(dest,"cox_all_coefficients.tsv"));save_tsv(rbindlist(curve_rows,fill=TRUE),file.path(dest,"continuous_and_timevarying_effects.tsv"));save_tsv(rbindlist(model_patients),file.path(dest,"model_patient_ids.tsv"))
capture.output(sessionInfo(),file=file.path(dest,"sessionInfo.txt"))
print(flow);print(models[sensitivity=="primary",.(endpoint,model,n_patients,n_events,HR,CI_low,CI_high,p,bootstrap_n_valid)])
