## write a little function to (1) select treatment matches and (2) among treatment matches, select confounder matches ##
match_main<-function(trtobs,trtcf,confounders,GPS_obs,GPS_cf,GPS_rcf,GPS_obs_sep,GPS_cf_sep,GPS_rcf_sep,trtdiff,S_inv){
  
  ## get inverse covariance matrix of counfounders to compute mahalanobis distance ##
  #sds<-apply(trtobs,2,sd) 
  subpop=1:5000
  all_trt_match_ind<-mapply(exmatch_trt_v2,split(trtcf,1:nrow(trtcf)),split(trtdiff,1:nrow(trtdiff)),MoreArgs=list(trtobs=trtobs))
  save(all_trt_match_ind,file='middle_match_data.RData')
  
  source(paste0('functions/match_conf_',modelst[mode],'.R'))
  source(paste0('functions/match_gps_',modelst[mode],'.R'))
  source(paste0('functions/match_sep_gps_',modelst[mode],'.R'))
  source(paste0('functions/match_rat_gps_',modelst[mode],'.R'))
  lambda=1
  m1.ind<-match_gps(trtobs,trtcf,trt_match=all_trt_match_ind,GPS_obs=GPS_obs,GPS_cf=GPS_cf,GPS_rcf=GPS_rcf,lambda=lambda)
  #save(m1.ind,file='exp_data.RData')
  m2.ind<-match_sep_gps(trtobs,trtcf,trt_match=all_trt_match_ind,GPS_obs_sep=GPS_obs_sep,GPS_cf_sep=GPS_cf_sep,GPS_rcf_sep=GPS_rcf_sep,lambda=lambda)
  #save(m2.ind,file='/Users/heyihui/Desktop/Wu\ documents/11.26meeting/data_test/final_match_condition_data.RData')
  m3.ind<-match_rat_gps(trtobs,trtcf,all_trt_match_ind,GPS_obs,GPS_cf,GPS_rcf,lambda)
  m0.ind<-mapply(match_conf,1:N,all_trt_match_ind,split(confounders,1:nrow(confounders)),MoreArgs=list(confounders=confounders,S_inv=S_inv))
  #save(m0.ind,file='/Users/heyihui/Desktop/Wu\ documents/11.26meeting/data_test/final_match_confounder_data.RData')
  
  ## construct a dataset with treatment indices in one column and their matched control indices in the other ##
  matchdat=NULL
  matchdat0<-NULL
  matchdat1<-NULL
  matchdat2<-NULL
  matchdat3=NULL
  for (i in 1:length(all_trt_match_ind)){
    if (sum(is.na(all_trt_match_ind[[i]]))==length(all_trt_match_ind[[i]])){}
    else{
      if (i %in% subpop){
        matchdat<-rbind(matchdat,cbind(rep(i,length(all_trt_match_ind[[i]])),all_trt_match_ind[[i]]))
      }else{}
    }
  }
  for (i in 1:length(m0.ind)){
    if (sum(is.na(m0.ind[[i]]))==length(m0.ind[[i]])){}
    else{
      matchdat0<-rbind(matchdat0,cbind(rep(i,length(m0.ind[[i]])),m0.ind[[i]]))
    }
  }
  for (i in 1:length(m1.ind)){
    if (sum(is.na(m1.ind[[i]]))==length(m1.ind[[i]])){}
    else{
      matchdat1<-rbind(matchdat1,cbind(rep(i,length(m1.ind[[i]])),m1.ind[[i]]))
    }
  }
  for (i in 1:length(m2.ind)){
    if (sum(is.na(m2.ind[[i]]))==length(m2.ind[[i]])){}
    else{
      matchdat2<-rbind(matchdat2,cbind(rep(i,length(m2.ind[[i]])),m2.ind[[i]]))
    }
  }
  for (i in 1:length(m3.ind)){
    if (sum(is.na(m3.ind[[i]]))==length(m3.ind[[i]])){}
    else{
      matchdat3<-rbind(matchdat3,cbind(rep(i,length(m3.ind[[i]])),m3.ind[[i]]))
    }
  }
  err=lambda**2*(log(GPS_cf/GPS_obs)[matchdat3[,2]]-log(GPS_obs/GPS_rcf)[matchdat3[,1]])**2+(1-lambda)**2*rowSums((trtobs[matchdat3[,2],]-trtcf[matchdat3[,1],])**2)
  up=quantile(err,0.95)
  lo=quantile(err,0.05)
  matchdat4=matchdat3[which((lambda**2*(log(GPS_cf/GPS_obs)[matchdat3[,2]]-log(GPS_obs/GPS_rcf)[matchdat3[,1]])**2+(1-lambda)**2*rowSums((trtobs[matchdat3[,2],]-trtcf[matchdat3[,1],])**2)>lo)
                      &(lambda**2*(log(GPS_cf/GPS_obs)[matchdat3[,2]]-log(GPS_obs/GPS_rcf)[matchdat3[,1]])**2+(1-lambda)**2*rowSums((trtobs[matchdat3[,2],]-trtcf[matchdat3[,1],])**2)<up)),]
  colnames(matchdat)<-c('trtind','ctlind') 
  colnames(matchdat0)<-c('trtind','ctlind')
  colnames(matchdat1)<-c('trtind','ctlind')
  colnames(matchdat2)<-c('trtind','ctlind')
  colnames(matchdat3)<-c('trtind','ctlind')
  colnames(matchdat4)<-c('trtind','ctlind')
  cat=sapply(all_trt_match_ind,function(x) length(x))
  cat0=sapply(m0.ind,function(x) length(x))
  cat1=sapply(m1.ind,function(x) length(x))
  cat2=sapply(m2.ind,function(x) length(x))
  cat3=sapply(m3.ind,function(x) length(x))
  cat4=cat3
  save(matchdat,matchdat0,matchdat1,matchdat2,matchdat3,matchdat4,cat,cat0,cat1,cat2,cat3,cat4,file='multi_analysis_match_data.RData')
  
}
