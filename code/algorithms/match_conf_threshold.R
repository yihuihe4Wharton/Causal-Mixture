match_conf<-function(trt_match,confounders,S_inv){
  X=confounders
  quan0=NULL
  for (j in 1:N){
    if (j %in% subpop){
    foo2=unlist(trt_match[j])
    if (length(foo2)==1){
      if (is.na(foo2)==1){}
      else{
          temp=matrix(rep(X[j,],length(foo2)),nrow=length(foo2),byrow=T)
          quan0<-c(quan0,sqrt(rowSums((temp-X[foo2,])%*%S_inv*(temp-X[foo2,]))))
      }
    }else{
        temp=matrix(rep(X[j,],length(foo2)),nrow=length(foo2),byrow=T)
        quan0<-c(quan0,sqrt(rowSums((temp-X[foo2,])%*%S_inv*(temp-X[foo2,]))))
    }
    }else{}
    }
  mdqtl=quantile(quan0,thr)
  ans=NULL
  for (j in 1:N){
  if (j %in% subpop){
  foo2=unlist(trt_match[j])
  if (length(foo2)==1){
    if (is.na(foo2)==1){
      m.ind=NA
    }
    else{
        temp=matrix(rep(X[j,],length(foo2)),nrow=length(foo2),byrow=T)
        dif=sqrt(rowSums((temp-X[foo2,])%*%S_inv*(temp-X[foo2,])))
        m.ind=unname(foo2[which(dif<=mdqtl)])
      }
    }
  else{
    temp=matrix(rep(X[j,],length(foo2)),nrow=length(foo2),byrow=T)
    dif=sqrt(rowSums((temp-X[foo2,])%*%S_inv*(temp-X[foo2,])))
    m.ind=unname(foo2[which(dif<=mdqtl)])
  }
  }else{
    m.ind=NA
  }
  ans=c(ans,list(m.ind))
  }
  return(ans)
  }
