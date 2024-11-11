match_gps<-function(trtobs,trtcf,trt_match,GPS_obs,GPS_cf,GPS_rcf,lambda){
  ans=NULL
  for (j in 1:length(trt_match)){
    foo2=unlist(trt_match[j])
  if (sum(is.na(foo2))==length(foo2)){
    m.ind<-NA
  } else{
    e1_dif<-abs(rep(GPS_cf[j],time=length(foo2))-GPS_obs[foo2])
    e2_dif<-abs(rep(GPS_obs[j],time=length(foo2))-GPS_rcf[foo2])
    w_dif=abs(trtobs[foo2,]-matrix(rep(trtcf[j,],length(foo2)),nrow=length(foo2),byrow=TRUE))
    tot_dif=NULL
    for (i in 1:length(foo2)){
      tot_dif=c(tot_dif,((1-lambda)**2*sum(w_dif[i,]**2)+lambda**2*(e1_dif[i]**2+e1_dif[i]**2)))
    }
    #MD_cut<-quantile(MD_i,probs=mdqtl)
    m.ind=unname(foo2[order(tot_dif)[1]])
    if (m.ind==j){
      m.ind=unname(foo2[order(tot_dif)[2]])}
  }
  ans=c(ans,m.ind)
  names(ans)[j]=j
}
return(ans)
}