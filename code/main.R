library(mvtnorm)
library(dplyr)
library(splines)
library(BayesTree)
library(xtable)
library(ggplot2)
library(mgcv)
library("dgof")
library(hrbrthemes)
library(ggsci)
library(ggridges)
library(GA)
library(sampling)
library(Ecume)
library(RANN)
setwd("/Users/heyihui/Desktop/conti-inference")
source(paste0('functions/algorithms/density_ratio_estimate.R'))
source(paste0('functions/algorithms/outcome_estimate.R'))
main=function(q,K,r=2,m=1,alpha){

### LOAD DATA ###
load(paste0("working_data/dataset",q,".Rdata"))

### SAMPLE CUTTING ###
pop_in=subp
N_in=length(pop_in)
  
### SAMPLE SPLITING ###
pop_split=split(sample(1:N),1:K)
pop_actual=sapply(1:K,function(k) which(pop_split[[k]] %in% pop_in))
rh=rep(0,N)
qh=matrix(0,nrow=N,ncol=2)
for (k in 1:K){

### ESTIMATE RATIO ###
rh_k=density_ratio_estimate(Eobs,Ercf,X,pop_split[[k]],pop_actual[[k]],alpha)
rh[pop_actual[[k]]]=rh_k

### ESTIMATE OUTCOME ###
qh_k=outcome_estimate(Eobs,Ecf,X,Yobs,pop_split[[k]],pop_actual[[k]])
qh[pop_actual[[k]],]=qh_k
}

### FINAL ESTIMATE ###
#point estimate
rlh=as.vector(rh*(Yobs-qh[,1])+qh[,2])
re=mean(rlh[pop_in])
#variance estimate
vh=mean(sapply(1:K,function(k) mean((rlh[pop_actual[[k]]]-rep(re,length(pop_actual[[k]])))**2)))
return(cbind(re,vh))
}

