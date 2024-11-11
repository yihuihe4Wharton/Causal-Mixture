library(mvtnorm)
library(dplyr)
library(splines)
library(BayesTree)
library(xtable)
library(ggplot2)
library(mgcv)
library("dgof")
library(hrbrthemes)
library(ggsci)
library(ggridges)
library(GA)
library(sampling)
library(Ecume)

setwd("/Users/heyihui/Desktop/2023.02.02 codes")
modelst=c("nearest neighbour","threshold")
for (md in 3:3){
if (md==1){
	mode=1
}else{
	mode=2
	if (md==2){
		thr=1
	}else{
		thr=0.05
	}
}
t=1
r=1
rlst=c("xgboost","linCDE")
for (ed in 2:7){
edlst=c(0.25,0.5,1,2,3,4,5)
edif=edlst[ed]
losslist=NULL
esslist=NULL
varlist=NULL
anslist=NULL
keeplist=NULL
covlist=NULL
for (q in 4:4){
# losslist=NULL
# esslist=NULL
# varlist=NULL
# anslist=NULL
# keeplist=NULL
# covlist=NULL
for (p in 1:2){
set.seed(p)
cat('=====================================\n')
print(paste("Exp",p))
cat('=====================================\n')
nboot<-50
calp<-0.3
source(paste0('create_dataset',q,'.R'))
source('functions/exmatch_trt_v2.R')
source(paste0("functions/estimate_gps_v",r,".R"))
source("functions/estimate_sep_gps_v1.R")
source(paste0('functions/match_main_',modelst[mode],'.R'))
## load the data for analysis ##
load('multi_analysis_data.RData')
P<-ncol(X)
S_inv<-solve(cov(X))
d=ncol(Eobs)
N=nrow(X)
Ecf=Eobs-matrix(rep(unit,N),nrow=N,byrow=TRUE)
Ecf_rev=Eobs+matrix(rep(unit,N),nrow=N,byrow=TRUE)
subpop=1:N

  ##############
  ## MATCHING ##
  ##############
  
  ##Estimating GPS
  print("Estimating GPS")
  #X:confounder
  estimate_gps(Eobs,X,Ecf,Ecf_rev)
  #GPS:two column matrix obs_score-cf_score
  estimate_sep_gps(Eobs,X,Ecf,Ecf_rev)
  #GPS_sep: 2d column matrix
  load('multi_analysis_GPS.RData')
  load('multi_analysis_GPS_sep.RData')
  
  ##Normalizing
  #g_min=min(min(gps_obs),min(gps_cf),min(gps_rcf))
  g_min=0
  g_max=max(max(gps_obs),max(gps_cf),max(gps_rcf))
  gps_obs=(gps_obs-rep(g_min,N))/(g_max-g_min)
  gps_cf=(gps_cf-rep(g_min,N))/(g_max-g_min)
  gps_rcf=(gps_rcf-rep(g_min,N))/(g_max-g_min)
  gsobs=NULL
  gscf=NULL
  gsrcf=NULL
  for (i in 1:d){
    gs_min=0
    gs_max=max(max(gps_obs_sep[,i]),max(gps_cf_sep[,i]),max(gps_rcf_sep[,i]))
    gsobs=cbind(gsobs,(gps_obs_sep[,i]-rep(gs_min,N))/(gs_max-gs_min))
    gscf=cbind(gscf,(gps_cf_sep[,i]-rep(gs_min,N))/(gs_max-gs_min))
    gsrcf=cbind(gsrcf,(gps_rcf_sep[,i]-rep(gs_min,N))/(gs_max-gs_min))
  }
  gps_obs_sep=gsobs
  gps_cf_sep=gscf
  gps_rcf_sep=gsrcf
  ## omega values, tolerances for exact matching on E ##
  #exposure's interval
  ecut<-apply(Ecf,2,sd)*N**(-1/(d+1))*edif
  xx<-matrix(ecut,nrow=nrow(Ecf),ncol=ncol(Ecf),byrow=T)
  
  # gps_plot
  print("gps_plot")
  data=data.frame("x1"=c(log(gps_cf/gps_obs),log(gps_obs/gps_rcf)),"x2"=factor(c(rep("c/o",5000),rep("o/r",5000))))
  ggplot(data=data, aes(x=x1, y=x2,fill=factor(stat(quantile)))) +
    stat_density_ridges(
      geom="density_ridges_gradient",
      calc_ecdf=TRUE,
      quantiles=c(0.1,0.9)
    )+
    scale_fill_manual(
      name="Probability",values=c("#E2EAF6","#436FB0","#E2eAF6")
    )+
    theme_bw()+
    theme(legend.position="none")+
    labs(x="x",y="Density")
  ggsave(paste0("dataset",q,"/matching=",md,"/edif=",edif,"/gps_estimate_=",rlst[r],"/gps_ratio_plot/p=",p,".png"))
  
  
  
  
  
  ## Carry Out Matching
 print("Carry Out Matching")
 match_main(trtobs=Eobs,trtcf=Ecf,confounders=X,GPS_obs=gps_obs,GPS_cf=gps_cf,GPS_rcf=gps_rcf,GPS_obs_sep=gps_obs_sep,GPS_cf_sep=gps_cf_sep,GPS_rcf_sep=gps_rcf_sep,trtdiff=xx,S_inv=S_inv)
  load("multi_analysis_match_data.RData")
  orig.matchind.0<-as.data.frame(matchdat0)
  keep.0<-unique(orig.matchind.0$trtind)
  ploss.0=mean(real_dif[keep.0])-man_dif
  
  orig.matchind.1<-as.data.frame(matchdat1)
  keep.1<-unique(orig.matchind.1$trtind)
  ploss.1=mean(real_dif[keep.1])-man_dif
  
  orig.matchind.2<-as.data.frame(matchdat2)
  keep.2<-unique(orig.matchind.2$trtind)
  ploss.2=mean(real_dif[keep.2])-man_dif
  
  orig.matchind.3<-as.data.frame(matchdat3)
  keep.3<-unique(orig.matchind.3$trtind)
  ploss.3=mean(real_dif[keep.3])-man_dif
  
  orig.matchind.4<-as.data.frame(matchdat4)
  keep.4<-unique(orig.matchind.4$trtind)
  ploss.4=mean(real_dif[keep.4])-man_dif
 
  orig.matchind<-as.data.frame(matchdat)
  keep<-unique(orig.matchind$trtind)
  ploss=mean(real_dif[keep])-man_dif
  
  ##ESS
  print("ESS")
  num_ma0=rep(0,N)
  for (i in 1:nrow(orig.matchind.0)){
    num_ma0[orig.matchind.0[i,2]]=num_ma0[orig.matchind.0[i,2]]+1/cat0[orig.matchind.0[i,1]]
  }
  ESS0=sum(num_ma0)**2/sum(num_ma0**2)

  num_ma1=rep(0,N)
  for (i in 1:nrow(orig.matchind.1)){
    num_ma1[orig.matchind.1[i,2]]=num_ma1[orig.matchind.1[i,2]]+1/cat1[orig.matchind.1[i,1]]
  }
  ESS1=sum(num_ma1)**2/sum(num_ma1**2)

  num_ma2=rep(0,N)
  for (i in 1:nrow(orig.matchind.2)){
    num_ma2[orig.matchind.2[i,2]]=num_ma2[orig.matchind.2[i,2]]+1/cat2[orig.matchind.2[i,1]]
  }
  ESS2=sum(num_ma2)**2/sum(num_ma2**2)

  num_ma3=rep(0,N)
  for (i in 1:nrow(orig.matchind.3)){
    num_ma3[orig.matchind.3[i,2]]=num_ma3[orig.matchind.3[i,2]]+1/cat3[orig.matchind.3[i,1]]
  }
  ESS3=sum(num_ma3)**2/sum(num_ma3**2)
  
  num_ma4=rep(0,N)
  for (i in 1:nrow(orig.matchind.4)){
    num_ma4[orig.matchind.4[i,2]]=num_ma4[orig.matchind.4[i,2]]+1/cat4[orig.matchind.4[i,1]]
  }
  ESS4=sum(num_ma4)**2/sum(num_ma4**2)

  num_ma=rep(0,N)
  for (i in 1:nrow(orig.matchind)){
    num_ma[orig.matchind[i,2]]=num_ma[orig.matchind[i,2]]+1/cat[orig.matchind[i,1]]
  }
  ESS=sum(num_ma)**2/sum(num_ma**2)

  ESSw=(sum(gps_rcf/gps_obs[keep.4]))**2/sum((gps_rcf/gps_obs[keep.4])**2)

 

  
  #KS
  print("KS")
  ks_stat = c(sapply(1:d, function(i) ks_test(as.numeric(Ecf[keep,i]), as.numeric(Eobs[,i]), w_x = rep(1, length(keep)), w_y = num_ma)),
              sapply(1:P, function(i) ks_test(as.numeric(X[keep,i]), as.numeric(X[,i]), w_x = rep(1, length(keep)), w_y = num_ma)),thresh = .001)
  ks=ks_stat[(1:(d+P))*5-4]
  KS=NULL
  for (i in 1:(d+P)){
    KS=c(KS,ks[[i]])
  }

  ks_stat = c(sapply(1:d, function(i) ks_test(as.numeric(Ecf[keep.0,i]), as.numeric(Eobs[,i]), w_x = rep(1, length(keep.0)), w_y = num_ma0)),
              sapply(1:P, function(i) ks_test(as.numeric(X[keep.0,i]), as.numeric(X[,i]), w_x = rep(1, length(keep.0)), w_y = num_ma0)),thresh = .001)
  ks=ks_stat[(1:(d+P))*5-4]
  KS0=NULL
  for (i in 1:(d+P)){
    KS0=c(KS0,ks[[i]])
  }

  ks_stat = c(sapply(1:d, function(i) ks_test(as.numeric(Ecf[keep.1,i]), as.numeric(Eobs[,i]), w_x = rep(1, length(keep.1)), w_y = num_ma1)),
              sapply(1:P, function(i) ks_test(as.numeric(X[keep.1,i]), as.numeric(X[,i]), w_x = rep(1, length(keep.1)), w_y = num_ma1)),thresh = .001)
  ks=ks_stat[(1:(d+P))*5-4]
  KS1=NULL
  for (i in 1:(d+P)){
    KS1=c(KS1,ks[[i]])
  }

  ks_stat = c(sapply(1:d, function(i) ks_test(as.numeric(Ecf[keep.2,i]), as.numeric(Eobs[,i]), w_x = rep(1, length(keep.2)), w_y = num_ma2)),
              sapply(1:P, function(i) ks_test(as.numeric(X[keep.2,i]), as.numeric(X[,i]), w_x = rep(1, length(keep.2)), w_y = num_ma2)),thresh = .001)
  ks=ks_stat[(1:(d+P))*5-4]
  KS2=NULL
  for (i in 1:(d+P)){
    KS2=c(KS2,ks[[i]])
  }

  ks_stat = c(sapply(1:d, function(i) ks_test(as.numeric(Ecf[keep.3,i]), as.numeric(Eobs[,i]), w_x = rep(1, length(keep.3)), w_y = num_ma3)),
              sapply(1:P, function(i) ks_test(as.numeric(X[keep.3,i]), as.numeric(X[,i]), w_x = rep(1, length(keep.3)), w_y = num_ma3)),thresh = .001)
  ks=ks_stat[(1:(d+P))*5-4]
  KS3=NULL
  for (i in 1:(d+P)){
    KS3=c(KS3,ks[[i]])
  }

  ks_stat = c(sapply(1:d, function(i) ks_test(as.numeric(Ecf[keep.4,i]), as.numeric(Eobs[,i]), w_x = rep(1, length(keep.4)), w_y = num_ma4)),
              sapply(1:P, function(i) ks_test(as.numeric(X[keep.4,i]), as.numeric(X[,i]), w_x = rep(1, length(keep.4)), w_y = num_ma4)),thresh = .001)
  ks=ks_stat[(1:(d+P))*5-4]
  KS4=NULL
  for (i in 1:(d+P)){
    KS4=c(KS4,ks[[i]])
  }

  #png(paste0("dataset",q,"/matching=",md,"/edif=",edif,"/gps_estimate_=",rlst[r],"/ks_plot/p=",p,".png"))
  #plot(KS,main=paste0("mode=",mode,";dataset=",q,";r=",r,";thr=",thr";edif=",edif,";ks_plot"),type = "o", col = 1, xlab = "Covariate component", ylab = "ks-statistics",xlim=c(1,d+P),ylim=c(0,0.2))
  #lines(KS0,type = "o", col = 2)
  #lines(KS1,type = "o", col = 3)
  #lines(KS2,type = "o", col = 4)
  #lines(KS3,type = "o", col = 5)
  #lines(KS4,type = "o", col = 6)
  #legend("topright",pch=c(3,3),legend=c("non-match","cov","gps","sep_gps","ratio","trimmed-rat"),col=c(1,2,3,4,5,6))
  #dev.off()
  
  ## Point Estimate ##
  print("Point Estimate")
  if (t==1){
    EYcf<-tapply(Yobs[orig.matchind$ctlind],orig.matchind$trtind,mean)
    est_ans<-mean(Yobs[as.numeric(names(EYcf))]-EYcf)
    EYcf<-tapply(Yobs[orig.matchind.0$ctlind],orig.matchind.0$trtind,mean)
    est_ans0<-mean(Yobs[as.numeric(names(EYcf))]-EYcf)
    EYcf<-tapply(Yobs[orig.matchind.1$ctlind],orig.matchind.1$trtind,mean)
    est_ans1<-mean(Yobs[as.numeric(names(EYcf))]-EYcf)
    EYcf<-tapply(Yobs[orig.matchind.2$ctlind],orig.matchind.2$trtind,mean)
    est_ans2<-mean(Yobs[as.numeric(names(EYcf))]-EYcf)
    EYcf<-tapply(Yobs[orig.matchind.3$ctlind],orig.matchind.3$trtind,mean)
    est_ans3<-mean(Yobs[as.numeric(names(EYcf))]-EYcf)
    EYcf<-tapply(Yobs[orig.matchind.4$ctlind],orig.matchind.4$trtind,mean)
    est_ans4<-mean(Yobs[as.numeric(names(EYcf))]-EYcf)
  }
  else if (t==2){
    mod="gaussian"
    fmla<-as.formula(paste0('y~',paste0('s(exp',1:d,')',collapse='+'),'+',paste0('s(conf',1:P,')',collapse='+')))
    fmla_lm<-as.formula(paste0('y~',paste0('exp',1:d,collapse='+'),'+',paste0('conf',1:P,collapse='+')))
    Yobs<-as.matrix(Yobs)
    regdat<-data.frame(Yobs,Eobs,X)
    names(regdat)<-c('y',paste0('exp',1:d),paste0('conf',1:P))
    fit.1<-bam(fmla, family=mod, data=regdat,discrete=TRUE)
  
    predXi<-data.frame(Ecf[orig.matchind$trtind,],X[orig.matchind$trtind,])
    predXj<-data.frame(Ecf[orig.matchind$trtind,],X[orig.matchind$ctlind,])
    names(predXi)<-c(paste0('exp',1:d),paste0('conf',1:P))
    names(predXj)<-c(paste0('exp',1:d),paste0('conf',1:P))
    ## take difference to be used for bias correction ##
    diffmeans<-predict(fit.1, predXi, type="response")-predict(fit.1, predXj, type="response")
    ## impute counterfactuals ##
    EYcf<-tapply(Yobs[orig.matchind$ctlind]+diffmeans,orig.matchind$trtind,mean)
    ## difference in imputed counterfactuals and observed number of deaths in each area ##
    est_ans<-mean(Yobs[as.numeric(names(EYcf))]-EYcf)
  
    predXi<-data.frame(Ecf[orig.matchind.1$trtind,],X[orig.matchind.1$trtind,])
    predXj<-data.frame(Ecf[orig.matchind.1$trtind,],X[orig.matchind.1$ctlind,])
    names(predXi)<-c(paste0('exp',1:d),paste0('conf',1:P))
    names(predXj)<-c(paste0('exp',1:d),paste0('conf',1:P))
    ## take difference to be used for bias correction ##
    diffmeans<-predict(fit.1, predXi, type="response")-predict(fit.1, predXj, type="response")
    ## impute counterfactuals ##
    EYcf<-tapply(Yobs[orig.matchind.1$ctlind]+diffmeans,orig.matchind.1$trtind,mean)
    ## difference in imputed counterfactuals and observed number of deaths in each area ##
    est_ans1<-mean(Yobs[as.numeric(names(EYcf))]-EYcf)
  
    predXi<-data.frame(Ecf[orig.matchind.0$trtind,],X[orig.matchind.0$trtind,])
    predXj<-data.frame(Ecf[orig.matchind.0$trtind,],X[orig.matchind.0$ctlind,])
    names(predXi)<-c(paste0('exp',1:d),paste0('conf',1:P))
    names(predXj)<-c(paste0('exp',1:d),paste0('conf',1:P))
    ## take difference to be used for bias correction ##
    diffmeans<-predict(fit.1, predXi, type="response")-predict(fit.1, predXj, type="response")
    ## impute counterfactuals ##
    EYcf<-tapply(Yobs[orig.matchind.0$ctlind]+diffmeans,orig.matchind.0$trtind,mean)
    ## difference in imputed counterfactuals and observed number of deaths in each area ##
    est_ans0<-mean(Yobs[as.numeric(names(EYcf))]-EYcf)
  
    predXi<-data.frame(Ecf[orig.matchind.2$trtind,],X[orig.matchind.2$trtind,])
    predXj<-data.frame(Ecf[orig.matchind.2$trtind,],X[orig.matchind.2$ctlind,])
    names(predXi)<-c(paste0('exp',1:d),paste0('conf',1:P))
    names(predXj)<-c(paste0('exp',1:d),paste0('conf',1:P))
    ## take difference to be used for bias correction ##
    diffmeans<-predict(fit.1, predXi, type="response")-predict(fit.1, predXj, type="response")
    ## impute counterfactuals ##
    EYcf<-tapply(Yobs[orig.matchind.2$ctlind]+diffmeans,orig.matchind.2$trtind,mean)
    ## difference in imputed counterfactuals and observed number of deaths in each area ##
    est_ans2<-mean(Yobs[as.numeric(names(EYcf))]-EYcf)
  
    predXi<-data.frame(Ecf[orig.matchind.3$trtind,],X[orig.matchind.3$trtind,])
    predXj<-data.frame(Ecf[orig.matchind.3$trtind,],X[orig.matchind.3$ctlind,])
    names(predXi)<-c(paste0('exp',1:d),paste0('conf',1:P))
    names(predXj)<-c(paste0('exp',1:d),paste0('conf',1:P))
    ## take difference to be used for bias correction ##
    diffmeans<-predict(fit.1, predXi, type="response")-predict(fit.1, predXj, type="response")
    ## impute counterfactuals ##
    EYcf<-tapply(Yobs[orig.matchind.3$ctlind]+diffmeans,orig.matchind.3$trtind,mean)
    est_ans3=mean(Yobs[as.numeric(names(EYcf))]-EYcf)
    
    predXi<-data.frame(Ecf[orig.matchind.4$trtind,],X[orig.matchind.4$trtind,])
    predXj<-data.frame(Ecf[orig.matchind.4$trtind,],X[orig.matchind.4$ctlind,])
    names(predXi)<-c(paste0('exp',1:d),paste0('conf',1:P))
    names(predXj)<-c(paste0('exp',1:d),paste0('conf',1:P))
    ## take difference to be used for bias correction ##
    diffmeans<-predict(fit.1, predXi, type="response")-predict(fit.1, predXj, type="response")
    ## impute counterfactuals ##
    EYcf<-tapply(Yobs[orig.matchind.4$ctlind]+diffmeans,orig.matchind.4$trtind,mean)
    # lower=quantile(log(gps_cf/gps_obs)[as.numeric(names(EYcf))],0.05)
    # upper=quantile(log(gps_cf/gps_obs)[as.numeric(names(EYcf))],0.95)
    # est_ans3<-weighted.mean(Yobs[as.numeric(names(EYcf))]-EYcf,weights=pnorm(log(gps_cf/gps_obs)[as.numeric(names(EYcf))]-lower,sd=0.1)*pnorm(upper-log(gps_cf/gps_obs)[as.numeric(names(EYcf))],sd=0.1))
    est_ans4=mean(Yobs[as.numeric(names(EYcf))]-EYcf)
  }
  est_answ=mean(Yobs-Yobs*gps_rcf/gps_obs)
    
  real_ans_origin=mean(real_dif)
  real_ans=mean(real_dif[keep])
  real_ans0=mean(real_dif[keep.0])
  real_ans1=mean(real_dif[keep.1])
  real_ans2=mean(real_dif[keep.2])
  real_ans3=mean(real_dif[keep.3])
  real_ans4=mean(real_dif[keep.4])
  
 print("Save Loss")
covlist=rbind(covlist,c(ploss,ploss.0,ploss.1,ploss.2,ploss.3,ploss.4)/man_dif)
  
  losslist=rbind(losslist,c((est_ans-real_ans)/real_ans,(est_ans0-real_ans0)/real_ans0,(est_ans1-real_ans1)/real_ans1,(est_ans2-real_ans2)/real_ans2,(est_ans3-real_ans3)/real_ans3,(est_ans4-real_ans4)/real_ans4,(est_answ-man_dif)/man_dif))
  anslist=rbind(anslist,c(real_ans_origin,est_ans,real_ans,est_ans0,real_ans0,est_ans1,real_ans1,est_ans2,real_ans2,est_ans3,real_ans3,est_ans4,real_ans4,est_answ))
  esslist=rbind(esslist,c(ESS,ESS0,ESS1,ESS2,ESS3,ESS4,ESSw))
  keeplist=rbind(keeplist,c(length(keep),length(keep.0),length(keep.1),length(keep.2),length(keep.3),length(keep.4)))
  save(anslist,losslist,esslist,keeplist,covlist,file=paste0("new_mat=",md,";est=",r,";edif=",edif,";simulation result.RData"))
}
}
finlist=(1+cbind(covlist,rep(0,nrow(covlist))))*losslist+cbind(covlist,rep(0,nrow(covlist)))
save(finlist,file=paste0("mat=",md,";est=",r,";edif=",edif,";final_result.RData"))
}
}
#colnames(anslist)<-c('real_ans_origin','est_ans','real_ans','est_ans0','real_ans0','est_ans1','real_ans1','est_ans2','real_ans2','est_ans3','real_ans3','est_ans4','real_ans4','est_answ') 
#save(losslist,esslist,keeplist,anslist,covlist,file = paste0('mode=',md,';r=',r,'_summary.RData'))
# unkeep=NULL
# for (i in 1:5000){
#   if (i %in% keep.3){}
#   else{
#     unkeep=c(unkeep,i)
#   }
# }

