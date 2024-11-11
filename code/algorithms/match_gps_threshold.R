match_gps<-function(trtobs,trtcf,trt_match,GPS_obs,GPS_cf,GPS_rcf,lambda=1){
  tot_dif=NULL
  for (j in 1:length(trt_match)){
    if (j %in% subpop){
    foo2=unlist(trt_match[j])
    if (length(foo2)==1){
      if (is.na(foo2)==1){}
      else{
        e1_dif<-abs(rep(GPS_cf[j],time=length(foo2))-GPS_obs[foo2])
        e2_dif<-abs(rep(GPS_obs[j],time=length(foo2))-GPS_rcf[foo2])
        tot_dif=c(tot_dif,e1_dif**2+e2_dif**2)
      }
      }
    else{
    e1_dif<-abs(rep(GPS_cf[j],time=length(foo2))-GPS_obs[foo2])
    e2_dif<-abs(rep(GPS_obs[j],time=length(foo2))-GPS_rcf[foo2])
    tot_dif=c(tot_dif,e1_dif**2+e2_dif**2)
    }
    }else{}
  }
  mdqtl=quantile(tot_dif,thr)  
  ans=NULL
  for (j in 1:N){
    if (j %in% subpop){
    foo2=unlist(trt_match[j])
    if (length(foo2)==1){
      if (is.na(foo2)==1){
        m.ind=NA}
      else{
        e1_dif<-abs(rep(GPS_cf[j],time=length(foo2))-GPS_obs[foo2])
        e2_dif<-abs(rep(GPS_obs[j],time=length(foo2))-GPS_rcf[foo2])
        #w_dif=abs(trtobs[foo2,]-matrix(rep(trtcf[j,],length(foo2)),nrow=length(foo2),byrow=TRUE))
        m.ind=unname(foo2[which(e1_dif**2+e2_dif**2<=mdqtl)])
      }
    }
    else{
    foo2=unlist(trt_match[j])
    e1_dif<-abs(rep(GPS_cf[j],time=length(foo2))-GPS_obs[foo2])
    e2_dif<-abs(rep(GPS_obs[j],time=length(foo2))-GPS_rcf[foo2])
    m.ind=unname(foo2[which(e1_dif**2+e2_dif**2<=mdqtl)])
    }
    }else{
      m.ind=NA
    }
  ans=c(ans,list(m.ind))
}
return(ans)
}