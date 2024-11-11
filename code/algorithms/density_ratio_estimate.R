density_ratio_estimate<-function(T1,T2,X,pop_sec,pop_test,alpha=1){
  N=nrow(X)
  d=ncol(X)+ncol(T1)
  train_split=split(sample((1:N)[!((1:N) %in% pop_sec)]),1:2)
  pop_n0=train_split[[1]]
  pop_n1=train_split[[2]]
  data_n0=cbind(T1,X)[pop_n0,]
  data_n1=cbind(T2,X)[pop_n1,]
  #M=5
  M=round(alpha*(nrow(data_n0))**(2/(2+d)))
  rmh=NULL
  id=nrow(data_n0)+1
  for (i in 1:length(pop_test)){
    match_result=nn2(rbind(data_n0,c(T1[pop_test[i],],X[pop_test[i],])),data_n1,k=M)
    #print(nrow(rbind(data_n0,c(T1[pop_test[i],],X[pop_test[i],]))))
    Kmh=sum(apply(match_result[[1]],1,function(x) id %in% x))
    rmh=c(rmh,length(pop_n0)/length(pop_n1)*Kmh/M)
  }
  return(rmh)
}