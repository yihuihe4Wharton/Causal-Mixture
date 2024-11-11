estimate_sep_gps<-function(trtobs,confounders,trtcf,trtcf_rev){
  library("xgboost")
  sink()
  n=N%/%2
  for (j in 1:ncol(trtobs)){
    if (j==1){
      GPS_mod1 <-xgboost(data = confounders[1:n,],
                         label = trtobs[1:n,1],nrounds=50,verbose=0)
      mod_mean=predict(GPS_mod1,confounders[(n+1):N,])
      mod_sd<- sd(trtobs[(n+1):N,1] - predict(GPS_mod1,confounders[(n+1):N,]))
      gps_obs1=dnorm(trtobs[(n+1):N,1],mean=mod_mean,sd=mod_sd)
      gps_cf1=dnorm(trtcf[(n+1):N,1],mean=mod_mean,sd=mod_sd)
      gps_rcf1=dnorm(trtcf_rev[(n+1):N,1],mean=mod_mean,sd=mod_sd)
      GPS_mod2 <-xgboost(data = confounders[(n+1):N,],
                         label = trtobs[(n+1):N,1],nrounds=50,verbose=0)
      mod_mean=predict(GPS_mod2,confounders[1:n,])
      mod_sd<- sd(trtobs[1:n,1] - predict(GPS_mod2,confounders[1:n,]))
      gps_obs2=dnorm(trtobs[1:n,1],mean=mod_mean,sd=mod_sd)
      gps_cf2=dnorm(trtcf[1:n,1],mean=mod_mean,sd=mod_sd)
      gps_rcf2=dnorm(trtcf_rev[1:n,1],mean=mod_mean,sd=mod_sd)
    } else {
      GPS_mod1 <-xgboost(data = cbind(confounders[1:n,],trtobs[1:n,1:(j-1)]),
                         label = trtobs[1:n,j],nrounds=50,verbose=0)
      mod_mean=predict(GPS_mod1,cbind(confounders[(n+1):N,],trtobs[(n+1):N,1:(j-1)]))
      mod_sd<- sd(trtobs[(n+1):N,j] - predict(GPS_mod1,cbind(confounders[(n+1):N,],trtobs[(n+1):N,1:(j-1)])))
      gps_obs1=cbind(gps_obs1,dnorm(trtobs[(n+1):N,j],mean=mod_mean,sd=mod_sd))
      gps_cf1=cbind(gps_cf1,dnorm(trtcf[(n+1):N,j],mean=mod_mean,sd=mod_sd))
      gps_rcf1=cbind(gps_rcf1,dnorm(trtcf_rev[(n+1):N,j],mean=mod_mean,sd=mod_sd))
      GPS_mod2 <-xgboost(data = cbind(confounders[(n+1):N,],trtobs[(n+1):N,1:(j-1)]),
                         label = trtobs[(n+1):N,j],nrounds=50,verbose=0)
      mod_mean=predict(GPS_mod2,cbind(confounders[1:n,],trtobs[1:n,1:(j-1)]))
      mod_sd<- sd(trtobs[1:n,j] - predict(GPS_mod2,cbind(confounders[1:n,],trtobs[1:n,1:(j-1)])))
      gps_obs2=cbind(gps_obs2,dnorm(trtobs[1:n,j],mean=mod_mean,sd=mod_sd))
      gps_cf2=cbind(gps_cf2,dnorm(trtcf[1:n,j],mean=mod_mean,sd=mod_sd))
      gps_rcf2=cbind(gps_rcf2,dnorm(trtcf_rev[1:n,j],mean=mod_mean,sd=mod_sd))
    }
  }
  gps_obs_sep=rbind(gps_obs2,gps_obs1)
  gps_cf_sep=rbind(gps_cf2,gps_cf1)
  gps_rcf_sep=rbind(gps_rcf2,gps_rcf1)
  save(gps_obs_sep,gps_cf_sep,gps_rcf_sep,file='multi_analysis_GPS_sep.RData')
}