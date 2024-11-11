match_rat_gps<-function(trtobs,trtcf,trt_match,GPS_obs,GPS_cf,GPS_rcf,lambda=1){
  subpop=1:5000
  tot_dif=NULL
  for (j in 1:length(trt_match)){
    if (j %in% subpop){
    foo2=unlist(trt_match[j])
    if (length(foo2)==1){
      if (is.na(foo2)==1){}
      else{
        ratio_dif<-abs(log(rep(GPS_cf[j]/GPS_obs[j],time=length(foo2)))-log(GPS_obs[foo2]/GPS_rcf[foo2]))
        tot_dif=c(tot_dif,ratio_dif)
      }
      }
    else{
      ratio_dif<-abs(log(rep(GPS_cf[j]/GPS_obs[j],time=length(foo2)))-log(GPS_obs[foo2]/GPS_rcf[foo2]))
      tot_dif=c(tot_dif,ratio_dif)
    }
    }else{}
  }
  mdqtl=quantile(tot_dif,thr)  
  #print(mdqtl)
  ans=NULL
  for (j in 1:length(trt_match)){
    if (j %in% subpop){
    foo2=unlist(trt_match[j])
    if (length(foo2)==1){
      if (is.na(foo2)==1){
        m.ind=NA}
      else{
        ratio_dif<-abs(log(rep(GPS_cf[j]/GPS_obs[j],time=length(foo2)))-log(GPS_obs[foo2]/GPS_rcf[foo2]))
        m.ind=unname(foo2[which(ratio_dif<=mdqtl)])
      }
    }
    else{
      ratio_dif<-abs(log(rep(GPS_cf[j]/GPS_obs[j],time=length(foo2)))-log(GPS_obs[foo2]/GPS_rcf[foo2]))
      m.ind=unname(foo2[which(ratio_dif<=mdqtl)])
    }
    }else{
      m.ind=NA
    }
  ans=c(ans,list(m.ind))
}
return(ans)
}