#Summarize the findings in simulation section
#Function to gather results
library(reshape2)
library(dplyr)
source("rpacksifter-2.R")
source("repSifter_reem.R")
source("funcs_sim.R")

sim_result_lmm <- function(zdslmm,dataori,datatest,
                           modelform=c("Y~x1+x2+x3+visit+(1|ID)","K~Y+x4+x5+visit+(1|ID)"),
                           modelbeta=list(c(-1,-0.5,-0.3,0.8),c(0.2,-1,0.2,2)),
                           gentype =c("l","n")){
  result <- data.frame()
  if(gentype=="l"){
    #models for original and ds data.
    ori_y <- lmer(modelform[1],data=dataori)
    ori_k <- lmer(modelform[2],data=dataori)
    ori_ciy <- confint(ori_y)[-c(1:3),]
    ori_cik <- confint(ori_k)[-c(1:3),]
    
    ds_lmm_y <- lmer(modelform[1],data=zdslmm)
    ds_lmm_k <- lmer(modelform[2],data=zdslmm)
    ds_lmm_ciy <- confint(ds_lmm_y)[-c(1:3),]
    ds_lmm_cik <- confint(ds_lmm_k)[-c(1:3),]
    
    #CIs
    ori_ciy_bool <- ori_ciy[,1]<=modelbeta[[1]]&ori_ciy[,2]>=modelbeta[[1]]
    ori_cik_bool <- ori_cik[,1]<=modelbeta[[2]]&ori_cik[,2]>=modelbeta[[2]]
    ds_lmm_ciy_bool <- ds_lmm_ciy[,1]<=modelbeta[[1]]&ds_lmm_ciy[,2]>=modelbeta[[1]]
    ds_lmm_cik_bool <- ds_lmm_cik[,1]<=modelbeta[[2]]&ds_lmm_cik[,2]>=modelbeta[[2]]
    result <- data.frame(cover=c(ori_ciy_bool,ori_cik_bool,ds_lmm_ciy_bool,
                                 ds_lmm_cik_bool),
                         datatype=rep(c("original","DataSifterLMM"),each=2*length(ori_cik_bool)),
                         variable=rep(rep(c("Y","K"),each=length(ori_cik_bool)),2)
    )
  }else{
    ori_y <- lmer(modelform[1],data=dataori)
    ori_k <- lmer(modelform[2],data=dataori)
    ds_lmm_y <- lmer(modelform[1],data=zdslmm)
    ds_lmm_k <- lmer(modelform[2],data=zdslmm)
  }
  #Prediction accuracy
  ori_predy <- predict(ori_y,newdata=datatest,allow.new.levels=TRUE)
  ori_predk <- predict(ori_k,newdata=datatest,allow.new.levels=TRUE)
  ds_lmm_predy <- predict(ds_lmm_y,newdata=datatest,allow.new.levels=TRUE)
  ds_lmm_predk <- predict(ds_lmm_k,newdata=datatest,allow.new.levels=TRUE)
  pred.acc <- data.frame(absdiff=c(abs(ori_predy-datatest$Y),
                                   abs(ori_predk-datatest$K),
                                   abs(ds_lmm_predy-datatest$Y),
                                   abs(ds_lmm_predk-datatest$K)),
                         datatype=rep(c("original","DataSifterLMM"),each=2*nrow(datatest)),
                         variable=rep(rep(c("Y","K"),each=nrow(datatest)),2))
  
  if(dim(result)[1]!=0){
    result <- list(result,pred.acc)
  }else{result=pred.acc}
  return(result)
}


gatherresult_lmm <- function(datalist,seed,testsample){
  testdatals <- data.gen(genmodel = "l",noisel = "s",sample=testsample,seed=seed) #Same test set seed
  testdatall <- data.gen(genmodel = "l",noisel = "l",sample=testsample,seed=seed) #Same test set seed
  testdatans <- data.gen(genmodel = "n",noisel = "s",sample=testsample,seed=seed) #Same test set seed
  testdatanl <- data.gen(genmodel = "n",noisel = "l",sample=testsample,seed=seed) #Same test set seed
  resultl_s <- sim_result_lmm(zdslmm=datalist[[1]][[2]],dataori=datalist[[1]][[1]],datatest=testdatals,
                          modelform=c("Y~x1+x2+x3+visit+(1|ID)","K~Y+x4+x5+visit+(1|ID)"),
                          modelbeta=list(c(-1,-0.5,-0.3,0.8),c(0.2,-1,0.2,2)),
                          gentype ="l")
  resultl_l <- sim_result_lmm(zdslmm=datalist[[2]][[2]],dataori=datalist[[2]][[1]],datatest=testdatall,
                          modelform=c("Y~x1+x2+x3+visit+(1|ID)","K~Y+x4+x5+visit+(1|ID)"),
                          modelbeta=list(c(-1,-0.5,-0.3,0.8),c(0.2,-1,0.2,2)),
                          gentype ="l")
  resultn_s <- sim_result_lmm(zdslmm=datalist[[3]][[2]],dataori=datalist[[3]][[1]],datatest=testdatans,
                          modelform=c("Y~sin(x1)+(x2)^2+x1*abs(x3)+visit+(1|ID)","K~sin(Y)+exp(cos(x4))+Y*abs(x5)+visit+(1|ID)"),
                          modelbeta=list(c(-1,-0.5,-0.3,0.8),c(0.2,-1,0.2,2)),
                          gentype ="n")
  resultn_l <- sim_result_lmm(zdslmm=datalist[[4]][[2]],dataori=datalist[[4]][[1]],datatest=testdatanl,
                          modelform=c("Y~sin(x1)+(x2)^2+x1*abs(x3)+visit+(1|ID)","K~sin(Y)+exp(cos(x4))+Y*abs(x5)+visit+(1|ID)"),
                          modelbeta=list(c(-1,-0.5,-0.3,0.8),c(0.2,-1,0.2,2)),
                          gentype ="n")
  temp <- list(resultl_s,resultl_l,resultn_s,resultn_l)
  covert1 <- cbind.data.frame(temp[[1]][[1]],covname = rep(c("x1","x2","x3","visit","Y","x4","x5","visit"),2),seed=seed,simtype="ls")
  covert1 <- dcast(data = covert1,formula=datatype+variable+seed+simtype ~ covname,value.var="cover")
  covert2 <- cbind.data.frame(temp[[2]][[1]],covname = rep(c("x1","x2","x3","visit","Y","x4","x5","visit"),2),seed=seed,simtype="ll")
  covert2 <- dcast(data = covert2,formula=datatype+variable+seed+simtype ~ covname,value.var="cover")
  covert <- rbind(covert1,covert2)
  
  absdiff1 <- cbind.data.frame(temp[[1]][[2]] %>% group_by(datatype,variable) %>% summarise(mean(absdiff)),seed=seed,simtype="ls")
  absdiff2 <- cbind.data.frame(temp[[2]][[2]] %>% group_by(datatype,variable) %>% summarise(mean(absdiff)),seed=seed,simtype="ll")
  absdiff3 <- cbind.data.frame(temp[[3]] %>% group_by(datatype,variable) %>% summarise(mean(absdiff)),seed=seed,simtype="ns")
  absdiff4 <- cbind.data.frame(temp[[4]] %>% group_by(datatype,variable) %>% summarise(mean(absdiff)),seed=seed,simtype="nl")
  
  absdiff <- rbind(absdiff1,absdiff2,absdiff3,absdiff4)
  names(absdiff)[3] = "absdiff"
  
  res <- absdiff %>% left_join(covert,by=c("datatype","variable","seed","simtype"))
  res$sample = testsample
  return(res)
}


#Load results from DSLMM
setwd("result/lmm")
filelist <- list.files()
result500 <- data.frame()
result1000 <- data.frame()
result2000 <- data.frame()
for(i in 1:length(filelist)){
  load(filelist[i])
  result500 <- rbind(result500,gatherresult_lmm(datalist=sim500,seed=12345,testsample=500))
  result1000 <- rbind(result1000,gatherresult_lmm(datalist=sim1000,seed=12345,testsample=1000))
  result2000 <- rbind(result2000,gatherresult_lmm(datalist=sim2000,seed=12345,testsample=2000))
}

#save(result500,result1000,result2000,file = "../full_result_0.2lmm.Rdata")
load("~/Box Sync/Longitudinal Sifter/MI_data_obfuscation/full_result_0.2lmm.Rdata")


r500_imp_lmm <- result500 %>% filter(datatype!="original")%>%group_by(datatype,variable,simtype) %>% summarise(meanabsdiff = mean(absdiff))
r500_ci_lmm <- result500 %>% filter(simtype %in% c("ll","ls")&datatype!="original")%>%group_by(datatype,variable,simtype) %>% 
  summarise(visit=mean(visit),x1=mean(x1,na.rm=TRUE),x2=mean(x2,na.rm=TRUE),x3=mean(x3,na.rm=TRUE),
            x4=mean(x4,na.rm=TRUE),x5=mean(x5,na.rm=TRUE),Y=mean(Y,na.rm=TRUE))

r1000_imp_lmm <- result1000 %>% filter(datatype!="original")%>%group_by(datatype,variable,simtype) %>% summarise(meanabsdiff = mean(absdiff))
r1000_ci_lmm <- result1000 %>% filter(simtype %in% c("ll","ls")&datatype!="original")%>%group_by(datatype,variable,simtype) %>% 
  summarise(visit=mean(visit),x1=mean(x1,na.rm=TRUE),x2=mean(x2,na.rm=TRUE),x3=mean(x3,na.rm=TRUE),
            x4=mean(x4,na.rm=TRUE),x5=mean(x5,na.rm=TRUE),Y=mean(Y,na.rm=TRUE))

r2000_imp_lmm <- result2000 %>% filter(datatype!="original")%>%group_by(datatype,variable,simtype) %>% summarise(meanabsdiff = mean(absdiff))
r2000_ci_lmm <- result2000 %>% filter(simtype %in% c("ll","ls")&datatype!="original")%>%group_by(datatype,variable,simtype) %>% 
  summarise(visit=mean(visit),x1=mean(x1,na.rm=TRUE),x2=mean(x2,na.rm=TRUE),x3=mean(x3,na.rm=TRUE),
            x4=mean(x4,na.rm=TRUE),x5=mean(x5,na.rm=TRUE),Y=mean(Y,na.rm=TRUE))


#Load results from MI, DSREEMtree
load("~/Box Sync/Longitudinal Sifter/MI_data_obfuscation/full_result_0.2.Rdata")

r500_imp <- result500 %>%filter(datatype!="DataSifterLMM") %>%group_by(datatype,variable,simtype) %>% summarise(meanabsdiff = mean(absdiff))
r500_imp <- rbind(r500_imp,r500_imp_lmm)
View(r500_imp)

r500_ci <- result500 %>% filter(simtype %in% c("ll","ls")&datatype!="DataSifterLMM")%>%group_by(datatype,variable,simtype) %>% 
  summarise(visit=mean(visit),x1=mean(x1,na.rm=TRUE),x2=mean(x2,na.rm=TRUE),x3=mean(x3,na.rm=TRUE),
            x4=mean(x4,na.rm=TRUE),x5=mean(x5,na.rm=TRUE),Y=mean(Y,na.rm=TRUE))
r500_ci <- rbind(r500_ci,r500_ci_lmm)
View(r500_ci)

r1000_imp <- result1000 %>%filter(datatype!="DataSifterLMM") %>%group_by(datatype,variable,simtype) %>% summarise(meanabsdiff = mean(absdiff))
r1000_imp <- rbind(r1000_imp,r1000_imp_lmm)
r1000_ci <- result1000 %>% filter(simtype %in% c("ll","ls")&datatype!="DataSifterLMM")%>%group_by(datatype,variable,simtype) %>% 
  summarise(visit=mean(visit),x1=mean(x1,na.rm=TRUE),x2=mean(x2,na.rm=TRUE),x3=mean(x3,na.rm=TRUE),
            x4=mean(x4,na.rm=TRUE),x5=mean(x5,na.rm=TRUE),Y=mean(Y,na.rm=TRUE))
r1000_ci <- rbind(r1000_ci,r1000_ci_lmm)
View(r1000_ci)

r2000_imp <- result2000 %>%filter(datatype!="DataSifterLMM") %>%group_by(datatype,variable,simtype) %>% summarise(meanabsdiff = mean(absdiff))
r2000_imp <- rbind(r2000_imp,r2000_imp_lmm)
r2000_ci <- result2000 %>% filter(simtype %in% c("ll","ls")&datatype!="DataSifterLMM")%>%group_by(datatype,variable,simtype) %>% 
  summarise(visit=mean(visit),x1=mean(x1,na.rm=TRUE),x2=mean(x2,na.rm=TRUE),x3=mean(x3,na.rm=TRUE),
            x4=mean(x4,na.rm=TRUE),x5=mean(x5,na.rm=TRUE),Y=mean(Y,na.rm=TRUE))
r2000_ci <- rbind(r2000_ci,r2000_ci_lmm)
View(r1000_imp)

########Extra testing case######
#REEM,MI, n=100,p=0.2 vs n=500,p=0.5
gatherresult <- function(datalist,seed,testsample){
  testdatals <- data.gen(genmodel = "l",noisel = "s",sample=testsample,seed=seed) #Same test set seed
  testdatall <- data.gen(genmodel = "l",noisel = "l",sample=testsample,seed=seed) #Same test set seed
  testdatans <- data.gen(genmodel = "n",noisel = "s",sample=testsample,seed=seed) #Same test set seed
  testdatanl <- data.gen(genmodel = "n",noisel = "l",sample=testsample,seed=seed) #Same test set seed
  resultl_s <- sim_result(zmi=datalist[[1]][[3]],zdsreem = datalist[[1]][[2]],zdslmm=datalist[[1]][[2]],dataori=datalist[[1]][[1]],datatest=testdatals,
                          modelform=c("Y~x1+x2+x3+visit+(1|ID)","K~Y+x4+x5+visit+(1|ID)"),
                          modelbeta=list(c(-1,-0.5,-0.3,0.8),c(0.2,-1,0.2,2)),
                          gentype ="l")
  resultl_l <- sim_result(zmi=datalist[[2]][[3]],zdsreem = datalist[[2]][[2]],zdslmm=datalist[[2]][[2]],dataori=datalist[[2]][[1]],datatest=testdatall,
                          modelform=c("Y~x1+x2+x3+visit+(1|ID)","K~Y+x4+x5+visit+(1|ID)"),
                          modelbeta=list(c(-1,-0.5,-0.3,0.8),c(0.2,-1,0.2,2)),
                          gentype ="l")
  resultn_s <- sim_result(zmi=datalist[[3]][[3]],zdsreem = datalist[[3]][[2]],zdslmm=datalist[[3]][[2]],dataori=datalist[[3]][[1]],datatest=testdatans,
                          modelform=c("Y~sin(x1)+(x2)^2+x1*abs(x3)+visit+(1|ID)","K~sin(Y)+exp(cos(x4))+Y*abs(x5)+visit+(1|ID)"),
                          modelbeta=list(c(-1,-0.5,-0.3,0.8),c(0.2,-1,0.2,2)),
                          gentype ="n")
  resultn_l <- sim_result(zmi=datalist[[4]][[3]],zdsreem = datalist[[4]][[2]],zdslmm=datalist[[4]][[2]],dataori=datalist[[4]][[1]],datatest=testdatanl,
                          modelform=c("Y~sin(x1)+(x2)^2+x1*abs(x3)+visit+(1|ID)","K~sin(Y)+exp(cos(x4))+Y*abs(x5)+visit+(1|ID)"),
                          modelbeta=list(c(-1,-0.5,-0.3,0.8),c(0.2,-1,0.2,2)),
                          gentype ="n")
  temp <- list(resultl_s,resultl_l,resultn_s,resultn_l)
  covert1 <- cbind.data.frame(temp[[1]][[1]],covname = rep(c("x1","x2","x3","visit","Y","x4","x5","visit"),4),seed=seed,simtype="ls")
  covert1 <- dcast(data = covert1,formula=datatype+variable+seed+simtype ~ covname,value.var="cover")
  covert2 <- cbind.data.frame(temp[[2]][[1]],covname = rep(c("x1","x2","x3","visit","Y","x4","x5","visit"),4),seed=seed,simtype="ll")
  covert2 <- dcast(data = covert2,formula=datatype+variable+seed+simtype ~ covname,value.var="cover")
  covert <- rbind(covert1,covert2)
  
  absdiff1 <- cbind.data.frame(temp[[1]][[2]] %>% group_by(datatype,variable) %>% summarise(mean(absdiff)),seed=seed,simtype="ls")
  absdiff2 <- cbind.data.frame(temp[[2]][[2]] %>% group_by(datatype,variable) %>% summarise(mean(absdiff)),seed=seed,simtype="ll")
  absdiff3 <- cbind.data.frame(temp[[3]] %>% group_by(datatype,variable) %>% summarise(mean(absdiff)),seed=seed,simtype="ns")
  absdiff4 <- cbind.data.frame(temp[[4]] %>% group_by(datatype,variable) %>% summarise(mean(absdiff)),seed=seed,simtype="nl")
  
  absdiff <- rbind(absdiff1,absdiff2,absdiff3,absdiff4)
  names(absdiff)[3] = "absdiff"
  
  res <- absdiff %>% left_join(covert,by=c("datatype","variable","seed","simtype"))
  res$sample = testsample
  return(res)
}


setwd("extra_sim")
filelist <- list.files()
result5005 <- data.frame()
result1002 <- data.frame()
for(i in 1:length(filelist)){
  load(filelist[i])
  result5005 <- rbind(result5005,gatherresult(datalist=sim5005,seed=12345,testsample=500))
  result1002 <- rbind(result1002,gatherresult(datalist=sim1002,seed=12345,testsample=100))
 print(i)
  }
save(result5005,result1002,file = "../sim_extra_reem.Rdata")

r500_imp <- result5005 %>%filter(datatype!="DataSifterLMM") %>%group_by(datatype,variable,simtype) %>% summarise(meanabsdiff = mean(absdiff))

View(r500_imp)

r500_ci <- result5005 %>% filter(simtype %in% c("ll","ls")&datatype!="DataSifterLMM")%>%group_by(datatype,variable,simtype) %>% 
  summarise(visit=mean(visit),x1=mean(x1,na.rm=TRUE),x2=mean(x2,na.rm=TRUE),x3=mean(x3,na.rm=TRUE),
            x4=mean(x4,na.rm=TRUE),x5=mean(x5,na.rm=TRUE),Y=mean(Y,na.rm=TRUE))
View(r500_ci)

r100_imp <- result1002 %>%filter(datatype!="DataSifterLMM") %>%group_by(datatype,variable,simtype) %>% summarise(meanabsdiff = mean(absdiff))

View(r100_imp)

r100_ci <- result1002 %>% filter(simtype %in% c("ll","ls")&datatype!="DataSifterLMM")%>%group_by(datatype,variable,simtype) %>% 
  summarise(visit=mean(visit),x1=mean(x1,na.rm=TRUE),x2=mean(x2,na.rm=TRUE),x3=mean(x3,na.rm=TRUE),
            x4=mean(x4,na.rm=TRUE),x5=mean(x5,na.rm=TRUE),Y=mean(Y,na.rm=TRUE))
View(r100_ci)

#lmm
setwd("lmm")
filelist <- list.files()
result5005_lmm <- data.frame()
result1002_lmm <- data.frame()
for(i in 1:length(filelist)){
  load(filelist[i])
  result5005_lmm <- rbind(result5005_lmm,gatherresult_lmm(datalist=sim5005,seed=12345,testsample=500))
  result1002_lmm <- rbind(result1002_lmm,gatherresult_lmm(datalist=sim1002,seed=12345,testsample=100))
  print(i)
}

r500_imp_lmm <- result5005_lmm %>% filter(datatype!="original")%>%group_by(datatype,variable,simtype) %>% summarise(meanabsdiff = mean(absdiff))
r500_ci_lmm <- result5005_lmm %>% filter(simtype %in% c("ll","ls")&datatype!="original")%>%group_by(datatype,variable,simtype) %>% 
  summarise(visit=mean(visit),x1=mean(x1,na.rm=TRUE),x2=mean(x2,na.rm=TRUE),x3=mean(x3,na.rm=TRUE),
            x4=mean(x4,na.rm=TRUE),x5=mean(x5,na.rm=TRUE),Y=mean(Y,na.rm=TRUE))

r100_imp_lmm <- result1002_lmm %>% filter(datatype!="original")%>%group_by(datatype,variable,simtype) %>% summarise(meanabsdiff = mean(absdiff))
r100_ci_lmm <- result1002_lmm %>% filter(simtype %in% c("ll","ls")&datatype!="original")%>%group_by(datatype,variable,simtype) %>% 
  summarise(visit=mean(visit),x1=mean(x1,na.rm=TRUE),x2=mean(x2,na.rm=TRUE),x3=mean(x3,na.rm=TRUE),
            x4=mean(x4,na.rm=TRUE),x5=mean(x5,na.rm=TRUE),Y=mean(Y,na.rm=TRUE))

save(r500_imp_lmm,r500_ci_lmm,r100_imp_lmm,r100_ci_lmm,file = "extra_results_lmm.Rdata")
