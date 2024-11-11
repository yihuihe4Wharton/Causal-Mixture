### HIGH DIMENSIONAL EXPOSURE, GAUSSION LINEAR IN ASSIGNMENT MODEL, NON-LINEAR IN OUTCOME MODEL
library(mvtnorm)
library(dplyr)
library(GA)
library(truncnorm)

N<-5000

X=rmvnorm(n=N,mean=rep(0.8,20),sigma=diag(20)*0.8+matrix(rep(0.2,400),nrow=20))

random_term=matrix(rtruncnorm(n=5*N,a=0,b=20,mean=10,sd=15),ncol=5)
beta=matrix(rep(c(rep(1,3),rep(0.4,17)),5),nrow=20,byrow=F)
Eobs=X%*%beta+random_term

unit=apply(Eobs,2,sd)/4
subp=which(apply(random_term-matrix(rep(unit,N),nrow=N,byrow=T),1,min)>=0)
Ecf=Eobs-matrix(rep(unit,N),nrow=N,byrow=TRUE)
Ercf=Eobs+matrix(rep(unit,N),nrow=N,byrow=TRUE)

alpha=matrix(c(5,10,7.5,5,5,rep(1,15),1,1,0.6,1,0.6),nrow=25,byrow=F)
Yobs=sqrt(abs(cbind(X,Eobs)%*%alpha+5))
Ycf=sqrt(abs(cbind(X,Ecf)%*%alpha+5))
real_dif=Yobs-Ycf
#man_dif=mean(real_dif[subp])
man_dif=mean(Ycf[subp])
#print(length(subp))
save(N,Eobs,Ecf,Ercf,unit,X,Yobs,real_dif,man_dif,subp,file="working_data/dataset4.Rdata")