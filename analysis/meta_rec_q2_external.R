#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(metafor)})
out<-'analysis_results/rec_q2_integration_2026-09-19/external'
components<-fread(file.path(out,'study_component_effects.tsv'))
genes<-fread(file.path(out,'study_member_effects.tsv'))
results<-list();omissions<-list();influences<-list();numerical_checks<-list()
checked_fit<-function(d,key) {
 warnings<-character()
 f<-withCallingHandlers(rma.uni(yi=d$estimate,vi=d$variance,method='REML',test='knha'),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 # Independently minimize the intercept-only restricted likelihood, including tau2=0.
 objective<-function(tau) {
  v<-d$variance+tau;w<-1/v;mu<-sum(w*d$estimate)/sum(w)
  .5*(sum(log(v))+log(sum(w))+sum(w*(d$estimate-mu)^2))
 }
 upper<-max(1,max(d$estimate^2)*10,max(d$variance)*10)
 grid<-c(0,exp(seq(log(1e-12),log(upper),length.out=500)))
 values<-vapply(grid,objective,numeric(1));i<-which.min(values)
 candidate<-optimize(objective,c(grid[max(1,i-1)],grid[min(length(grid),i+1)]),tol=1e-12)$minimum
 tau<-if(objective(0)<=objective(candidate))0 else candidate
 gain<-objective(f$tau2)-objective(tau)
 if(gain>1e-7)f<-rma.uni(yi=d$estimate,vi=d$variance,tau2=tau,method='REML',test='knha')
 stopifnot(abs(objective(f$tau2)-objective(tau))<1e-6)
 numerical_checks[[key]]<<-data.table(fit=key,metafor_tau2=f$tau2,independent_tau2=tau,objective_difference=objective(f$tau2)-objective(tau),numerical_refit=gain>1e-7,warning=paste(warnings,collapse='; '))
 f
}
fit_meta<-function(d,scoring,model,endpoint) {
 d<-d[study!='Ogden_2025' & n_patients>=2 & is.finite(estimate) & is.finite(variance) & variance>0]
 if(nrow(d)<2)return(data.table(scoring=scoring,model=model,endpoint=endpoint,k=nrow(d),estimate=NA_real_,low=NA_real_,high=NA_real_,p=NA_real_))
 key<-paste(scoring,model,endpoint)
 f<-checked_fit(d,key)
 for(s in d$study) {
  dd<-d[study!=s]
  if(nrow(dd)<2)next
  ff<-checked_fit(dd,paste(key,'omit',s))
  omissions[[paste(key,s)]]<<-data.table(scoring=scoring,model=model,endpoint=endpoint,excluded_study=s,k=nrow(dd),estimate=as.numeric(ff$b),low=ff$ci.lb,high=ff$ci.ub,p=ff$pval)
 }
 influences[[key]]<<-data.table(scoring=scoring,model=model,endpoint=endpoint,study=d$study,weight_percent=as.numeric(weights(f)),study_estimate=d$estimate,study_variance=d$variance)
 data.table(scoring=scoring,model=model,endpoint=endpoint,k=nrow(d),n_patients=sum(d$n_patients),estimate=as.numeric(f$b),se=f$se,low=f$ci.lb,high=f$ci.ub,p=f$pval,tau2=f$tau2,I2=f$I2,Q=f$QE,Q_p=f$QEp,df=f$k-f$p)
}
for(sc in unique(components$scoring))for(en in unique(components$endpoint)) {
 d<-components[scoring==sc & model=='common' & endpoint==en]
 results[[paste(sc,en)]]<-fit_meta(d,sc,'common',en)
}
for(sc in unique(genes$scoring))for(g in unique(genes$gene)) {
 d<-genes[scoring==sc & model=='common' & gene==g]
 results[[paste(sc,g)]]<-fit_meta(d,sc,'common',g)
}
res<-rbindlist(results,fill=TRUE);res[,q_member_family:=NA_real_]
for(sc in unique(res$scoring)) {
 ix<-which(res$scoring==sc & res$endpoint %in% unique(genes$gene))
 res$q_member_family[ix]<-p.adjust(res$p[ix],method='BH',n=30L)
}
fwrite(res,file.path(out,'random_effects_meta.tsv'),sep='\t',na='NA')
fwrite(rbindlist(omissions,fill=TRUE),file.path(out,'leave_one_study.tsv'),sep='\t',na='NA')
fwrite(rbindlist(influences,fill=TRUE),file.path(out,'meta_weights.tsv'),sep='\t')
fwrite(rbindlist(numerical_checks,fill=TRUE),file.path(out,'meta_numerical_checks.tsv'),sep='\t',na='NA')
capture.output(sessionInfo(),file=file.path(out,'meta_session_info.txt'))
print(res[endpoint %in% c('Claudin_minus_Polarity','Claudin','Polarity')])
