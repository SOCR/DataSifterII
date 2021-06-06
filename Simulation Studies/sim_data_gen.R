#Data generation
args <- commandArgs(TRUE)
source("rpacksifter-2.R")
source("repSifter_reem.R")
source("sim_funcs.R")
library(reshape2)

sim.data <- function(sample,mispct,seed){
  #Linear model noiselevel = s
  datals <- data.gen(genmodel = "l",noisel = "s",sample=sample,seed=seed)
  dsls <- repSifter(data=datals,mispct=mispct,misw="perij",lnames=c("Y","K"),timevar="visit",ID="ID",reem = FALSE,maxit = 3)
  dsreemls <- repSifter(data=datals,mispct=mispct,misw="perij",lnames=c("Y","K"),timevar="visit",ID="ID",reem = TRUE,maxit = 10)
  mils <- mi2l.pan(data=datals,mispct=mispct)
  
  
  #Linear model noiselevel = l
  datall <- data.gen(genmodel = "l",noisel = "l",sample=sample,seed=seed)
  dsll <- repSifter(data=datall,mispct=mispct,misw="perij",lnames=c("Y","K"),timevar="visit",ID="ID",reem = FALSE,maxit = 3)
  dsreemll <- repSifter(data=datall,mispct=mispct,misw="perij",lnames=c("Y","K"),timevar="visit",ID="ID",reem = TRUE,maxit = 10)
  mill <- mi2l.pan(data=datall,mispct=mispct)
  
  
  #NonLinear model noiselevel = s
  datans <- data.gen(genmodel = "n",noisel = "s",sample=sample,seed=seed)
  dsns <- repSifter(data=datans,mispct=mispct,misw="perij",lnames=c("Y","K"),timevar="visit",ID="ID",reem = FALSE,maxit = 3)
  dsreemns <- repSifter(data=datans,mispct=mispct,misw="perij",lnames=c("Y","K"),timevar="visit",ID="ID",reem = TRUE,maxit = 10)
  mins <- mi2l.pan(data=datans,mispct=mispct)
  
  
  #NonLinear model noiselevel = l
  datanl <- data.gen(genmodel = "n",noisel = "l",sample=sample,seed=seed)
  dsnl <- repSifter(data=datanl,mispct=mispct,misw="perij",lnames=c("Y","K"),timevar="visit",ID="ID",reem = FALSE,maxit = 3)
  dsreemnl <- repSifter(data=datanl,mispct=mispct,misw="perij",lnames=c("Y","K"),timevar="visit",ID="ID",reem = TRUE,maxit = 10)
  minl <- mi2l.pan(data=datanl,mispct=mispct)
  
  
  return(list(ls=list(raw=datals,dsreem=dsreemls[[1]],dsglmm=dsls[[1]],mils = mils),
              ll=list(raw=datall,dsreem=dsreemll[[1]],dsglmm=dsll[[1]],mils = mill),
              ns=list(raw=datans,dsreem=dsreemns[[1]],dsglmm=dsns[[1]],mils = mins),
              nl=list(raw=datanl,dsreem=dsreemnl[[1]],dsglmm=dsnl[[1]],mils = minl)
  ))
}

library(doParallel)
library(foreach)
library(parallel)
library(doParallel)
cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)
clusterExport(cl, list(as.vector(lsf.str())),envir=environment())
clusterEvalQ(cl, c(library("dplyr"),library("randomForest"),library("plyr"),library("reshape2"),
                  library("glmmLasso"),
                 library("glmnet"),
                  library("lme4"),
                  library("missForest"),
                  library("REEMtree"),
                  library("MASS"),
                  library("nlme")))




result <- foreach(j = c(1:100),.combine='rbind') %dopar%{
sim5002 <- sim.data(sample=500,mispct=0.2,seed=j)
sim10002 <- sim.data(sample=1000,mispct=0.2,seed=j)
sim <- rbind(sim5002,sim10002)
}

save(sim5002,sim10002,file = paste("simdata_mis",args,".Rdata",sep = ""))

