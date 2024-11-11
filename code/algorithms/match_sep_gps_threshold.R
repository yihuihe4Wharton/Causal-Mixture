match_sep_gps<-function(trtobs,trtcf,trt_match,GPS_obs_sep,GPS_cf_sep,GPS_rcf_sep,lambda=1){
  #foo2 comes from trt_match, foo3 is the GPS of to-be-matched unit
  quan=NULL
  for (j in 1:length(trt_match)){
    if (j %in% subpop){
    foo2=unlist(trt_match[j])
    if (length(foo2)==1){
      if (is.na(foo2)==1){}
      else{
        e1_dif=0
        e2_dif=0
        for (i in 1:d){
          e1_dif=e1_dif+abs(rep(GPS_cf_sep[j,i],time=length(foo2))-GPS_obs_sep[foo2,i])**2
          e2_dif=e2_dif+abs(rep(GPS_obs_sep[j,i],time=length(foo2))-GPS_rcf_sep[foo2,i])**2
        }
        quan=c(quan,e1_dif+e2_dif)
      }
    }else{
      e1_dif=0
      e2_dif=0
      for (i in 1:d){
        e1_dif=e1_dif+abs(rep(GPS_cf_sep[j,i],time=length(foo2))-GPS_obs_sep[foo2,i])**2
        e2_dif=e2_dif+abs(rep(GPS_obs_sep[j,i],time=length(foo2))-GPS_rcf_sep[foo2,i])**2
      }
      quan=c(quan,e1_dif+e2_dif)
    }
    }else{}
  }
  mdqtl=quantile(quan,thr)
  ans=NULL
  for (j in 1:length(trt_match)){
    if (j %in% subpop){
    foo2=unlist(trt_match[j])
    if (length(foo2)==1){
      if (is.na(foo2)==1){
        m.ind=NA
      }
      else{
        e1_dif=0
        e2_dif=0
        for (i in 1:d){
          e1_dif=e1_dif+abs(rep(GPS_cf_sep[j,i],time=length(foo2))-GPS_obs_sep[foo2,i])**2
          e2_dif=e2_dif+abs(rep(GPS_obs_sep[j,i],time=length(foo2))-GPS_rcf_sep[foo2,i])**2
        }
        m.ind=unname(foo2[which(e1_dif+e2_dif<mdqtl)])
      }
    }else{
      e1_dif=0
      e2_dif=0
      for (i in 1:d){
        e1_dif=e1_dif+abs(rep(GPS_cf_sep[j,i],time=length(foo2))-GPS_obs_sep[foo2,i])**2
        e2_dif=e2_dif+abs(rep(GPS_obs_sep[j,i],time=length(foo2))-GPS_rcf_sep[foo2,i])**2
      }
      #MD_cut<-quantile(MD_i,probs=mdqtl)
      m.ind=unname(foo2[which(e1_dif+e2_dif<mdqtl)])
    }
    }else{
      m.ind=NA
    }
    ans=c(ans,list(m.ind))
    names(ans)[j] <- j
  }
  return(ans)
}