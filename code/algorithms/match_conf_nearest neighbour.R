match_conf<-function(j,foo2,foo3,confounders,S_inv){
  foo2=unlist(foo2)
  if (sum(is.na(foo2))==length(foo2)){
    m.ind<-NA
  } else{
    temp<-matrix(foo3,nrow=nrow(confounders),ncol=ncol(confounders),byrow=T)
    
    MD_i<-sqrt(rowSums(as.matrix(temp-confounders)%*%S_inv*(as.matrix(temp-confounders))))
    #MD_cut<-quantile(MD_i,probs=mdqtl)
    m.ind=unname(foo2[order(MD_i[foo2])[1]])
    if (m.ind==j){
      m.ind=unname(foo2[order(MD_i[foo2])[2]])}
  }
  return(m.ind)
}