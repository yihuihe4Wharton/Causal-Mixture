outcome_estimate<-function(T1,T2,X,Y,pop_sec,pop_test){
  mod="gaussian"
  N=nrow(T1)
  pop_train=(1:N)[!((1:N) %in% pop_sec)]
  d1=ncol(T1[pop_train,])
  d2=ncol(X)
  fmla<-as.formula(paste0('y~',paste0('s(exp',1:d1,')',collapse='+'),'+',paste0('s(conf',1:d2,')',collapse='+')))
  regdat<-data.frame(Y[pop_train,],T1[pop_train,],X[pop_train,])
  names(regdat)<-c('y',paste0('exp',1:d1),paste0('conf',1:d2))
  fit.1<-bam(fmla, family=mod, data=regdat)
  case1<-data.frame(T1[pop_test,],X[pop_test,])
  names(case1)<-c(paste0('exp',1:d1),paste0('conf',1:d2))
  case2<-data.frame(T2[pop_test,],X[pop_test,])
  names(case2)<-c(paste0('exp',1:d1),paste0('conf',1:d2))
  return(cbind(predict(fit.1, case1),predict(fit.1, case2)))
}