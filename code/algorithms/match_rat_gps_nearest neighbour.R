match_rat_gps<-function(trtobs,trtcf,trt_match,GPS_obs,GPS_cf,GPS_rcf,lambda=0.01){
  ans=NULL
  subpop=1:5000
  for (j in 1:length(trt_match)){
    if (j %in% subpop){
    foo2=unlist(trt_match[j])
    if (length(foo2)==1){
      if (is.na(foo2)==1){
        m.ind=NA}
      else{
        w_dif=abs(trtobs[foo2,]-matrix(rep(trtcf[j,],length(foo2)),nrow=length(foo2),byrow=TRUE))
        ratio_dif<-lambda**2*(log(rep(GPS_cf[j]/GPS_obs[j],time=length(foo2)))-log(GPS_obs[foo2]/GPS_rcf[foo2]))**2+(1-lambda)**2*rowSums(w_dif**2)
        m.ind=unname(foo2[order(ratio_dif)[1]])
        if (m.ind==j){
          m.ind=unname(foo2[order(ratio_dif)[2]])
        }
      }
    }
    else{
      w_dif=abs(trtobs[foo2,]-matrix(rep(trtcf[j,],length(foo2)),nrow=length(foo2),byrow=TRUE))
      ratio_dif<-lambda**2*(log(rep(GPS_cf[j]/GPS_obs[j],time=length(foo2)))-log(GPS_obs[foo2]/GPS_rcf[foo2]))**2+(1-lambda)**2*rowSums(w_dif**2)
      m.ind=unname(foo2[order(ratio_dif)[1]])
      if (m.ind==j){
        m.ind=unname(foo2[order(ratio_dif)[2]])
      }
    }
    }else{
      m.ind=NA
    }
  ans=c(ans,list(m.ind))
  }
  
return(ans)
}