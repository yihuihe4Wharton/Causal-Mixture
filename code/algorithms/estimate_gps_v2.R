estimate_gps<-function(trtobs,confounders,trtcf,trtcf_rev){
  library("xgboost")
  library(LinCDE)
  sink()
  n=N%/%2
  for (j in 1:ncol(trtobs)){
    if (j==1){
      model1= LinCDE.boost(X = confounders[1:n,], y = trtobs[1:n,1],verbose=F)
      gps_obs1=predict.LinCDE(model1,X=unname(confounders[(n+1):N,]),y=trtobs[(n+1):N,1])
      gps_cf1=predict.LinCDE(model1,X=unname(confounders[(n+1):N,]),y=trtcf[(n+1):N,1])
      gps_rcf1=predict.LinCDE(model1,X=unname(confounders[(n+1):N,]),y=trtcf_rev[(n+1):N,1])
      model2= LinCDE.boost(X = confounders[(n+1):N,], y = trtobs[(n+1):N,1], verbose = F)
      gps_obs2=predict.LinCDE(model2,X=unname(confounders[1:n,]),y=trtobs[1:n,1])
      gps_cf2=predict.LinCDE(model2,X=unname(confounders[1:n,]),y=trtcf[1:n,1])
      gps_rcf2=predict.LinCDE(model2,X=unname(confounders[1:n,]),y=trtcf_rev[1:n,1])
    } else {
      model1= LinCDE.boost(X = cbind(confounders[1:n,],trtobs[1:n,1:(j-1)]), y = trtobs[1:n,j], verbose = F)
      gps_obs1=gps_obs1*predict.LinCDE(model1,X=unname(cbind(confounders[(n+1):N,],trtobs[(n+1):N,1:(j-1)])),y=trtobs[(n+1):N,j])
      gps_cf1=gps_cf1*predict.LinCDE(model1,X=unname(cbind(confounders[(n+1):N,],trtcf[(n+1):N,1:(j-1)])),y=trtcf[(n+1):N,j])
      gps_rcf1=gps_rcf1*predict.LinCDE(model1,X=unname(cbind(confounders[(n+1):N,],trtcf_rev[(n+1):N,1:(j-1)])),y=trtcf_rev[(n+1):N,j])
      model2= LinCDE.boost(X = cbind(confounders[(n+1):N,],trtobs[(n+1):N,1:(j-1)]), y = trtobs[(n+1):N,j], verbose = F)
      gps_obs2=gps_obs2*predict.LinCDE(model2,X=unname(cbind(confounders[1:n,],trtobs[1:n,1:(j-1)])),y=trtobs[1:n,j])
      gps_cf2=gps_cf2*predict.LinCDE(model2,X=unname(cbind(confounders[1:n,],trtcf[1:n,1:(j-1)])),y=trtcf[1:n,j])
      gps_rcf2=gps_rcf2*predict.LinCDE(model2,X=unname(cbind(confounders[1:n,],trtcf_rev[1:n,1:(j-1)])),y=trtcf_rev[1:n,j])
    }
  }
  gps_obs=c(gps_obs2,gps_obs1)
  gps_cf=c(gps_cf2,gps_cf1)
  gps_rcf=c(gps_rcf2,gps_rcf1)
  save(gps_obs,gps_cf,gps_rcf,file='multi_analysis_GPS.RData')
}