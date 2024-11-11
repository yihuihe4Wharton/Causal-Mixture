source("functions/main.R")
### CREATE DATASET ###
for (q in 1:1){
for (i in 1:1){
set.seed(i)
source(paste0('functions/data_setup/create_dataset',q,'.R'))
load(paste0("working_data/dataset",q,".Rdata"))
true=man_dif
### DO ESTIMATE ###
result=NULL
for (a in 5:5){
  result=rbind(result,main(q=q,r=1,K=5,alpha=0.125*2**a))
}
### PRESENT ERROR ###
#print(result)
#err=mean(apply(result,1,function(x) (x[1]-man_dif)**2))
#cov=mean(apply(result,1,function(x) abs(x[1]-man_dif)<1.96*x[2]**0.5))
### PRESENT COVERAGE RATE ###
}
}