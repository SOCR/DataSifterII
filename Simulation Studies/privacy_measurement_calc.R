source("rpacksifter-2.R")
source("repSifter_reem.R")
source("funcs_sim.R")
datals <- data.gen(genmodel = "l",noisel = "s",sample=500,seed=123)
mils <- mi2l.pan(data=datals,mispct=0.2)


#Function to recognize the static location in multiple imputed datasets.
static_loc <- function(mi_object){
  allimp <- complete(mi_object,action ="long")
  s_locY <- allimp %>% group_by(.id)%>% summarise(static = (length(unique(Y))==1))
  s_locK <- allimp %>% group_by(.id)%>% summarise(static = (length(unique(K))==1))
  return(list(Y=s_locY,K=s_locK))
}


#Function to calculate the abs(hat.K - K) on average. 
#The hat.K is based on whether the location has static values in MI data 
#hat.K only depends on x4,x5 and visit for nonstatic locations and depends on Y,x4,x5 on static locations.
#for DataSifter all cells are nonstatic.

dis_risk <-function(rawdata,miobj,siftdata,siftdata1,simtype=c("l","n"),testrn=50){
#Datasets ready
  n=nrow(rawdata)
  modelform <- ifelse(simtype=="l","K~Y+x4+x5+visit+(1|ID)","K~sin(Y)+exp(cos(x4))+Y*abs(x5)+visit+(1|ID)")
  modelformY <- ifelse(simtype=="l","Y~x1+x2+x3+visit+(1|ID)","Y~sin(x1)+(x2)^2+x1*abs(x3)+visit+(1|ID)")
  mi_sloc<-static_loc(miobj)
  mi_slocY <- mi_sloc$Y[1:testrn,]
  mi_slocK <- mi_sloc$K[1:testrn,]
  
  midata <- complete(miobj,action ="long")
  midata <- midata %>% mutate(rowid=rep(1:n,max(midata$.imp)))
  midata <- midata %>% filter(rowid %in% c(1:testrn))
  siftdata <- siftdata %>% mutate(rowid=1:n)
  siftdata$ID <- as.numeric(siftdata$ID)
  siftdata <- siftdata[1:testrn,]
  siftdata1 <- siftdata1 %>% mutate(rowid=1:n)
  siftdata1$ID <- as.numeric(siftdata1$ID)
  siftdata1 <- siftdata1[1:testrn,]
#Construct models based on D_{-i}  
  models <- sapply(1:testrn, function(x) lmer(modelform,data=rawdata[-x,]))
  modelsY <- sapply(1:testrn, function(x) lmer(modelformY,data=rawdata[-x,]))
    #{sid <- sample(c(1:n)[-x],floor(0.5*n))
  rawdata <- rawdata[1:testrn,]
#predictions based on models
  cidrow <- which(names(midata)=="rowid")
  pred_mi <- data.frame(rowid=midata$rowid,k_i=apply(midata,1,function(x) {
    x<-data.frame(matrix(x,nrow = 1))
    colnames(x)<- colnames(midata)
    if(mi_slocK$static[x[,cidrow]]==TRUE){
      rawdata$K[x[,cidrow]]
    }else{
      predict(models[[x[,cidrow]]],newdata=x,allow.new.levels=TRUE)}}),
    y_i=apply(midata,1,function(x) {
      x<-data.frame(matrix(x,nrow = 1))
      colnames(x)<- colnames(midata)
      if(mi_slocY$static[x[,cidrow]]==TRUE){
        rawdata$Y[x[,cidrow]]
      }else{
        predict(modelsY[[x[,cidrow]]],newdata=x,allow.new.levels=TRUE)}}))
  pred_mi <- pred_mi %>% group_by(rowid) %>% summarise(mik_i=mean(k_i),miy_i=mean(y_i))
  cidrowsif <- which(names(siftdata)=="rowid")
  pred_ds <- data.frame(rowid=siftdata$rowid,K=rawdata$K,Y=rawdata$Y,dsk_i=apply(siftdata,1,function(x) {
    x<-data.frame(matrix(x,nrow = 1))
    names(x)<- names(siftdata)
    predict(models[[unlist(x[,cidrowsif])]],newdata=x,allow.new.levels=TRUE)}),
    dsy_i=apply(siftdata,1,function(x) {
      x<-data.frame(matrix(x,nrow = 1))
      names(x)<- names(siftdata)
      predict(modelsY[[unlist(x[,cidrowsif])]],newdata=x,allow.new.levels=TRUE)}))
  pred_ds1 <- data.frame(rowid=siftdata1$rowid,ds1k_i=apply(siftdata1,1,function(x) {
    x<-data.frame(matrix(x,nrow = 1))
    names(x)<- names(siftdata1)
    predict(models[[unlist(x[,cidrowsif])]],newdata=x,allow.new.levels=TRUE)}),
    ds1y_i=apply(siftdata1,1,function(x) {
      x<-data.frame(matrix(x,nrow = 1))
      names(x)<- names(siftdata)
      predict(modelsY[[unlist(x[,cidrowsif])]],newdata=x,allow.new.levels=TRUE)}))
  predfull <- pred_ds %>% left_join(pred_ds1,by="rowid") %>% left_join(pred_mi,by="rowid") %>% 
              mutate(midr_K=abs(mik_i-K),dsdr_K=abs(dsk_i-K),ds1dr_K=abs(ds1k_i-K),
                     midr_Y=abs(miy_i-Y),dsdr_Y=abs(dsy_i-Y),ds1dr_Y=abs(ds1y_i-Y))
  return(predfull)
} 

####TEST######
mispct=0.2
sample=500
seed=123
datals <- data.gen(genmodel = "l",noisel = "s",sample=sample,seed=seed)
dsreemls <- repSifter(data=datals,mispct=mispct,misw="perij",lnames=c("Y","K"),timevar="visit",ID="ID",reem = TRUE,maxit = 3)
mils <- mi2l.pan(data=datals,mispct=mispct)
dr <- dis_risk(rawdata = datals,miobj=mils,siftdata = dsreemls[[1]],siftdata1=dsreemls[[1]],simtype = "l")
mean(dr$midr)
mean(dr$dsdr)
mean(dr$ds1dr)

####Actual numbers#####
setwd("~/Box Sync/Longitudinal Sifter/MI_data_obfuscation/result/")

lsdr_full_500 <- data.frame()
lsdr_full_1000 <- data.frame()
lldr_full_500 <- data.frame()
lldr_full_1000 <- data.frame()
nsdr_full_500 <- data.frame()
nsdr_full_1000 <- data.frame()
nldr_full_500 <- data.frame()
nldr_full_1000 <- data.frame()

for(i in c(1:100)){
  load(paste("simdata",i,".Rdata",sep = "")) #load sim500 sim1000 sim2000
  #load n=500
  lsraw500<-sim500$ls$raw
  llraw500<-sim500$ll$raw
  nsraw500<-sim500$ns$raw
  nlraw500<-sim500$nl$raw
  
  lsreem500<-sim500$ls$dsreem
  llreem500<-sim500$ll$dsreem
  nsreem500<-sim500$ns$dsreem
  nlreem500<-sim500$nl$dsreem
  
  lsmi500<-sim500$ls$mi
  llmi500<-sim500$ll$mi
  nsmi500<-sim500$ns$mi
  nlmi500<-sim500$nl$mi
  
  #load n=1000
  lsraw1000<-sim1000$ls$raw
  llraw1000<-sim1000$ll$raw
  nsraw1000<-sim1000$ns$raw
  nlraw1000<-sim1000$nl$raw
  
  lsreem1000<-sim1000$ls$dsreem
  llreem1000<-sim1000$ll$dsreem
  nsreem1000<-sim1000$ns$dsreem
  nlreem1000<-sim1000$nl$dsreem
  
  lsmi1000<-sim1000$ls$mi
  llmi1000<-sim1000$ll$mi
  nsmi1000<-sim1000$ns$mi
  nlmi1000<-sim1000$nl$mi
  load(paste("lmm/simdata_lmm",i,".Rdata",sep = ""))
  lslmm500<-sim500$ls$dslmm
  lllmm500<-sim500$ll$dslmm
  nslmm500<-sim500$ns$dslmm
  nllmm500<-sim500$nl$dslmm
  
  lslmm1000<-sim1000$ls$dslmm
  lllmm1000<-sim1000$ll$dslmm
  nslmm1000<-sim1000$ns$dslmm
  nllmm1000<-sim1000$nl$dslmm
  #Calculate DR for n=500
  lsdr500 <- dis_risk(rawdata = lsraw500,miobj=lsmi500,siftdata = lsreem500,siftdata1=lslmm500,simtype = "l")
  lldr500 <- dis_risk(rawdata = llraw500,miobj=llmi500,siftdata = llreem500,siftdata1=lllmm500,simtype = "l")
  nsdr500 <- dis_risk(rawdata = nsraw500,miobj=nsmi500,siftdata = nsreem500,siftdata1=nslmm500,simtype = "n")
  nldr500 <- dis_risk(rawdata = nlraw500,miobj=nlmi500,siftdata = nlreem500,siftdata1=nllmm500,simtype = "n")
  lsdr_full_500<- rbind(lsdr_full_500,cbind(lsdr500,rep=i))
  lldr_full_500<- rbind(lldr_full_500,cbind(lldr500,rep=i))
  nsdr_full_500<- rbind(nsdr_full_500,cbind(nsdr500,rep=i))
  nldr_full_500<- rbind(nldr_full_500,cbind(nldr500,rep=i))
  #Calculate DR for n=100
  lsdr1000 <- dis_risk(rawdata = lsraw1000,miobj=lsmi1000,siftdata = lsreem1000,siftdata1=lslmm1000,simtype = "l")
  lldr1000 <- dis_risk(rawdata = llraw1000,miobj=llmi1000,siftdata = llreem1000,siftdata1=lllmm1000,simtype = "l")
  nsdr1000 <- dis_risk(rawdata = nsraw1000,miobj=nsmi1000,siftdata = nsreem1000,siftdata1=nslmm1000,simtype = "n")
  nldr1000 <- dis_risk(rawdata = nlraw1000,miobj=nlmi1000,siftdata = nlreem1000,siftdata1=nllmm1000,simtype = "n")
  lsdr_full_1000<- rbind(lsdr_full_1000,cbind(lsdr1000,rep=i))
  lldr_full_1000<- rbind(lldr_full_1000,cbind(lldr1000,rep=i))
  nsdr_full_1000<- rbind(nsdr_full_1000,cbind(nsdr1000,rep=i))
  nldr_full_1000<- rbind(nldr_full_1000,cbind(nldr1000,rep=i))
  print(i)
}

#save(lsdr_full_500,
#     lsdr_full_1000,
#     lldr_full_500,
#     lldr_full_1000,
#     nsdr_full_500,
#     nsdr_full_1000,
#     nldr_full_500,
#     nldr_full_1000,file = "../disclosurerisk.Rdata")#1-100

namechar <- list(c("lsdr","500"),
           c("lsdr","1000"),
           c("lldr","500"),
           c("lldr","1000"),
           c("nsdr","500"),
           c("nsdr","1000"),
           c("nldr","500"),
           c("nldr","1000"))
j=1
library(dplyr)
fulldata <- data.frame()
temp <- lapply(list(lsdr_full_500,lsdr_full_1000,lldr_full_500,lldr_full_1000,nsdr_full_500,nsdr_full_1000,nldr_full_500,nldr_full_1000),function(i){
  i%>% group_by(rep) %>% summarize(midr_K=mean(midr_K),midr_Y=mean(midr_Y),dsreem_K=mean(dsdr_K),dsreem_Y=mean(dsdr_Y),dsglmm_K=mean(ds1dr_K),dsglmm_Y=mean(ds1dr_Y)) 
  #temp$setting=namechar[[j]][1]
  #temp$n=namechar[[j]][2]
  #fulldata <- rbind(fulldata,temp)
  #j=j+1
})

for(i in c(1:length(namechar))){
  temp[[i]]$setting=namechar[[i]][1]
  temp[[i]]$n=namechar[[i]][2]
}

full_data <- do.call(rbind,temp)

#save(full_data,file = "../full_privacy_data.Rdata")
##############Plotting#################
library(ggplot2)
library(reshape2)
library(dplyr)
load("full_privacy_data.Rdata")
full_data$setting<-ifelse(full_data$setting=="lsdr","Linear Small Noise",
                          ifelse(full_data$setting=="lldr","Linear Large Noise",
                                 ifelse(full_data$setting=="nsdr","Nonlinear Small Noise",
                                        "Nonlinear Large Noise")))

full_data_long <- melt(full_data,id.vars=c("rep","setting","n"))
full_data_long <- full_data_long %>% mutate(lvar = sub("^[a-z]*_","",variable),
                          model = sub("(dr)?_[A-Z]","",variable))

full_data_long$setting <- factor(full_data_long$setting,levels=unique(full_data_long$setting))
full_data_long$n <- factor(full_data_long$n,levels=c(500,1000),labels = c("n=500","n=1000"))

ggplot(data = full_data_long,aes(x=lvar,y=value,fill=model))+geom_boxplot(lwd=0.2)+
  facet_grid(n~setting,scales="free_y")+
  scale_fill_manual(name="Method:",values = c("green","lightblue","red"))+
  xlab("Longitudinal Variable Type")+
  ylab("Mean Privacy Measurement")+
  theme_bw()

View(full_data_long %>% group_by(setting,n,variable) %>% summarise(mean(value)))
