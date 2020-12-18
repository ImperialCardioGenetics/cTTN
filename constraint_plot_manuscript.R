rm(list=ls())
load("./rough_TTN_exp.RData")
load("./rough_TTN_obs.RData")

get_constraint<-function(v_list){
  list=data.frame(POS=v_list)
  list=merge(list,TTN_exp,by="POS",all.x=TRUE)
  list<-subset(list,is.na(predicted_prop)==FALSE)
  list<-merge(list,rough_TTN_obs,by=c("CHROM","POS","REF","ALT"),all.x = TRUE)
  list$high_coverage<-ifelse(list$median.x>=40,"high","low")
  list$obs<-ifelse(list$AF<0.001,1,0)
  
  load("./coverage_correction_factor.RData")
  
  list<-merge(list,correction_factor,by.x = "median.x",by.y="coverage",all.x = TRUE)
  list$exp_prop<-ifelse(list$high_coverage=="low",list$predicted_prop*list$cor_factor_exp,list$predicted_prop)
  p=list$exp_prop
  
  exp_sample<-apply(sapply(p,function(x){rbinom(10000,size=1,prob=x)}),1,sum)
  sd=sd(exp_sample)
  
  obs<-sum(list$obs,na.rm = TRUE)
  p_value<-mean(exp_sample<=obs)
  exp<-sum(list$exp_prop,na.rm = TRUE)

  oe<-obs/exp
  X<-sort(exp_sample)
  e_cdf<-1:length(X)/length(X)
  exp_ci<-c(X[which(e_cdf>0.95)[1]],X[which(e_cdf>=0.05)[1]])
  oe_ci<-obs/exp_ci
  result<-c(obs,exp,oe,oe_ci[1],oe_ci[2],p_value)
  result
}



result<-data.frame(Obs_sub=NA,Exp_sub=NA,oe=NA,oe_ci_lower=NA,oe_ci_upper=NA,p_value=NA)

load("./list_var/list_1.RData")
result<-rbind(result,get_constraint(list_1))

load("./list_var/ls_motif.RData")
result<-rbind(result,get_constraint(ls_motif))

load("./list_var/ls_norm.RData")
result<-rbind(result,get_constraint(ls_norm))

load("./list_var/ls_all.RData")
result<-rbind(result,get_constraint(ls_all))

load("./list3.RData")
result<-rbind(result,get_constraint(list_3))

result<-na.omit(result)

result$list<-c("Motif_BSJ","LinearSplicing_circRNA","LinearSplicing_non-circRNA","LinearSplicing_all","Motif_Exon")

TTN_exp<-subset(TTN_exp,SYMBOL=="TTN")

get_constraint_consq<-function(TTN_var){
  list<-subset(TTN_var,is.na(predicted_prop)==FALSE)
  list<-merge(list,rough_TTN_obs,by=c("CHROM","POS","REF","ALT"),all.x = TRUE)
  list$high_coverage<-ifelse(list$median.x>=40,"high","low")
  list$obs<-ifelse(list$AF<0.001,1,0)
  
  list<-merge(list,correction_factor,by.x = "median.x",by.y="coverage",all.x = TRUE)
  list$exp_prop<-ifelse(list$high_coverage=="low",list$predicted_prop*list$cor_factor_exp,list$predicted_prop)
  p=list$exp_prop
  exp_sample<-apply(sapply(p,function(x){rbinom(10000,size=1,prob=x)}),1,sum)
  sd=sd(exp_sample)
  #var<-sum(list$obs*list$exp_prop*(1-list$exp_prop),na.rm = TRUE)
  obs<-sum(list$obs,na.rm = TRUE)
  p_value<-mean(exp_sample<=obs)
  exp<-sum(list$exp_prop,na.rm = TRUE)
  oe<-obs/exp
  X<-sort(exp_sample)
  e_cdf<-1:length(X)/length(X)
  #plot(X,e_cdf,type="s")
  exp_ci<-c(X[which(e_cdf>0.95)[1]],X[which(e_cdf>=0.05)[1]])
  oe_ci<-obs/exp_ci
  result<-c(obs,exp,oe,oe_ci[1],oe_ci[2],p_value)
  result
}


tv<-c("frameshift_variant","splice_acceptor_variant",
      "splice_donor_variant","stop_gained")
TTN_tv<-subset(TTN_exp,is.element(Consequence,tv))

tv_result<-c(get_constraint_consq(TTN_tv),"TTNtv")

TTN_mis<-subset(TTN_exp,is.element(Consequence,"missense_variant"))
mis_result<-c(get_constraint_consq(TTN_mis),"TTN-mis")

TTN_syn<-subset(TTN_exp,is.element(Consequence,"synonymous_variant"))
syn_result<-c(get_constraint_consq(TTN_syn),"TTN-syn")

result<-rbind(result,mis_result,syn_result,tv_result)

library(ggplot2)
result$oe<-round(as.numeric(result$oe),2)
result$oe_ci_lower<-round(as.numeric(result$oe_ci_lower),2)
result$oe_ci_upper<-round(as.numeric(result$oe_ci_upper),2)
ggplot(result,aes(x=list,y=oe))+geom_point(size=2,shape=23)+geom_errorbar(aes(ymin=oe_ci_lower,ymax=oe_ci_upper,width=0.2))+geom_hline(yintercept=1.0, linetype="dashed", color = "red") +labs(y= "Observed/Expected") + theme_classic()

