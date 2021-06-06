#Functions for simulations
require(mice)
require(lattice)
require(pan)
library(dplyr)
library(glmmLasso)
library(glmnet)
library(doParallel)
library(lme4)
library(missForest)
library(REEMtree)
library(MASS)
library(nlme)
#Generate simulation data based on
#1. generative model type; 2. noise level; 3. patient sample size; 4. seed.
data.gen <- function(genmodel = c("l","n"),noisel = c("s","l"),sample,seed){
  if(genmodel=="l"){
  set.seed(seed)
  x1 <- rnorm(sample,mean=1,sd=5)
  x2 <- rnorm(sample,mean=4,sd=10)
  x3 <- rnorm(sample,mean=2,sd=4)
  x4 <- rnorm(sample,mean=6,sd=20)
  x5 <- rnorm(sample,mean=10,sd=30)
  
  
  ni <- sample(c(1:10),sample,replace = TRUE)
  visit <- unlist(sapply(ni,function(x) 1:x))
  e <- rnorm(length(visit),0,2)
  e1 <- rnorm(length(visit),0,2)
  ID <- unlist(sapply(1:sample,function(x) rep(x,ni[x])))
  visitdata<- data.frame(ID,visit,e,e1)
  
  #Generate bi and eij, j=1,2,3,4
  b <- rnorm(sample,0,1)
  b1 <- rnorm(sample,0,1)
  ID <- c(1:sample)
  
  #Combining
  sim1<-data.frame(ID,x1,x2,x3,x4,x5,b,b1)
  sim1 <- sim1 %>% right_join(visitdata,by="ID")
  sim1 <- sim1 %>% mutate(Y=1-x1-0.5*x2-0.3*x3+0.8*visit+b+e)
  sim1 <- sim1 %>% mutate(K=-15+0.2*Y-x4+0.2*x5+2*visit+b1+e1)
  
  if(noisel=="s"){
  noise <- data.frame(sapply(1:5,function(x) rnorm(dim(sim1)[1],x*2.73,10)))
  names(noise) <- paste("x",c(6:(5+5)),sep = "")
  }else{
  noise <- data.frame(sapply(1:20,function(x) rnorm(dim(sim1)[1],x*2.73,10)))
  names(noise) <- paste("x",c(6:(5+20)),sep = "")
  }
  sim1 <- cbind(sim1,noise)
  sim1 <- sim1 %>% dplyr::select(-c(b1,e1,e,b))
  }else{
    set.seed(seed)
    x1 <- runif(sample,min=-5,max=5)
    x2 <- runif(sample,min=-14,max=20)
    x3 <- runif(sample,min=-4,max=8)
    x4 <- runif(sample,min=-10,max=20)
    x5 <- runif(sample,min=-20,max=30)
    #Try uniform distribution and F distribution.
    ni <- sample(c(1:10),sample,replace = TRUE)
    visit <- unlist(sapply(ni,function(x) 1:x))
    #e <- rnorm(length(visit),0,4)
    #e1 <- rnorm(length(visit),0,4)
    e <- rnorm(length(visit),0,8)
    e1 <- rnorm(length(visit),0,8)
    ID <- unlist(sapply(1:sample,function(x) rep(x,ni[x])))
    visitdata<- data.frame(ID,visit,e,e1)
    
    #Generate bi and eij, j=1,2,3,4
    b <- rnorm(sample,0,3)
    b1 <- rnorm(sample,0,4)
    ID <- c(1:sample)
    
    #Combining
    sim1<-data.frame(ID,x1,x2,x3,x4,x5,b,b1)
    sim1 <- sim1 %>% right_join(visitdata,by="ID")
    sim1 <- sim1 %>% mutate(Y=10+3*sin(x1)-0.2*(x2)^2-0.1*x1*abs(x3)+1*visit+b+e)
    sim1 <- sim1 %>% mutate(K=2+0.05*sin(Y)+0.4*exp(cos(x4))-0.02*Y*abs(x5)+2*visit+b1+e1)
    
    if(noisel=="s"){
      noise <- data.frame(sapply(1:5,function(x) rnorm(dim(sim1)[1],x*2.73,10)))
      names(noise) <- paste("x",c(6:(5+5)),sep = "")
    }else{
      noise <- data.frame(sapply(1:20,function(x) rnorm(dim(sim1)[1],x*2.73,10)))
      names(noise) <- paste("x",c(6:(5+20)),sep = "")
    }
    sim1 <- cbind(sim1,noise)
    sim1 <- sim1 %>% dplyr::select(-c(b1,e1,e,b))
  }
  
  misindY<-numeric()
  ni <- sim1 %>% group_by(ID) %>% dplyr::mutate(ni=n())
  ni <- ni %>%  filter(row_number()==1) %>% dplyr::select(ID,ni) %>% mutate(b1=rnorm(1,0,2))
  b1 <- do.call("c",sapply(1:nrow(ni), function(x) rep(ni$b1[x],ni$ni[x])))
  set.seed(seed)
  for(i in 1:dim(sim1)[1]){misindY[i]<-rbinom(1,1,1/(1+exp(2-10*sim1$x1[i]+10*sim1$visit[i]-b1[i])))}
  table(misindY)
  
  misindK<-numeric()
  for(i in 1:dim(sim1)[1]){misindK[i]<-rbinom(1,1,1/(1+exp(-3+4*sim1$x5[i]+12*sim1$visit[i]-b1[i])))}
  table(misindK)
  
  sim1$Y[(misindY==1)]<-rep(NA,length(which(misindY==1)))
  sim1$K[(misindK==1)]<-rep(NA,length(which(misindK==1)))
  return(sim1)
}



#2. generate missing and impute back with MI mutilevel model 2l.pan create 5 complete datasets.
mi2l.pan <- function(data,mispct){
  dmis <- data
  lnames <- c("Y","K")
  dmis[,lnames] <- prodNA(data.frame(dmis[,lnames]),noNA = mispct)
  ini <- mice(dmis, maxit = 0)
  pred <- ini$pred
  pred["Y", ] <- c(-2, 1, 1, 1, 1, 1, 2,0,0,rep(2,ncol(dmis)-9)) 
  pred["K", ] <- c(-2, 1, 1, 1, 1, 1, 2,2,0,rep(2,ncol(dmis)-9)) 
  meth <- ini$method
  meth[c(8,9)] <- c("2l.pan","2l.pan")
  imp5 <- mice(dmis, pred = pred, meth = meth, print = FALSE)
  return(imp5)
}

#3. comparing results
sim_result <- function(zmi,zdsreem,zdslmm,dataori,datatest,
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
  
  ds_reem_y <- lmer(modelform[1],data=zdsreem)
  ds_reem_k <- lmer(modelform[2],data=zdsreem)
  ds_lmm_y <- lmer(modelform[1],data=zdslmm)
  ds_lmm_k <- lmer(modelform[2],data=zdslmm)
  ds_reem_ciy <- confint(ds_reem_y)[-c(1:3),]
  ds_reem_cik <- confint(ds_reem_k)[-c(1:3),]
  ds_lmm_ciy <- confint(ds_lmm_y)[-c(1:3),]
  ds_lmm_cik <- confint(ds_lmm_k)[-c(1:3),]
  
  
  #pooled model for MI results 
  mi_y <- with(zmi,lmer(modelform[1]))
  mi_k <- with(zmi,lmer(modelform[2]))
  mi_ciy <- summary(pool(mi_y),conf.int=TRUE)[-1,c("2.5 %","97.5 %")]
  mi_cik <- summary(pool(mi_k),conf.int=TRUE)[-1,c("2.5 %","97.5 %")]
  
  #CIs
  ori_ciy_bool <- ori_ciy[,1]<=modelbeta[[1]]&ori_ciy[,2]>=modelbeta[[1]]
  ori_cik_bool <- ori_cik[,1]<=modelbeta[[2]]&ori_cik[,2]>=modelbeta[[2]]
  ds_reem_ciy_bool <- ds_reem_ciy[,1]<=modelbeta[[1]]&ds_reem_ciy[,2]>=modelbeta[[1]]
  ds_reem_cik_bool <- ds_reem_cik[,1]<=modelbeta[[2]]&ds_reem_cik[,2]>=modelbeta[[2]]
  ds_lmm_ciy_bool <- ds_lmm_ciy[,1]<=modelbeta[[1]]&ds_lmm_ciy[,2]>=modelbeta[[1]]
  ds_lmm_cik_bool <- ds_lmm_cik[,1]<=modelbeta[[2]]&ds_lmm_cik[,2]>=modelbeta[[2]]
  mi_ciy_bool <- mi_ciy[,1]<=modelbeta[[1]]&mi_ciy[,2]>=modelbeta[[1]]
  mi_cik_bool <- mi_cik[,1]<=modelbeta[[2]]&mi_cik[,2]>=modelbeta[[2]]
  result <- data.frame(cover=c(ori_ciy_bool,ori_cik_bool,ds_reem_ciy_bool,
                                ds_reem_cik_bool,ds_lmm_ciy_bool,
                               ds_lmm_cik_bool,mi_ciy_bool,mi_cik_bool),
                        datatype=rep(c("original","DataSifterREEM","DataSifterLMM","MI"),each=2*length(ori_cik_bool)),
                        variable=rep(rep(c("Y","K"),each=length(ori_cik_bool)),4)
                        )
  }else{
    ori_y <- lmer(modelform[1],data=dataori)
    ori_k <- lmer(modelform[2],data=dataori)
    ds_reem_y <- lmer(modelform[1],data=zdsreem)
    ds_reem_k <- lmer(modelform[2],data=zdsreem)
    ds_lmm_y <- lmer(modelform[1],data=zdslmm)
    ds_lmm_k <- lmer(modelform[2],data=zdslmm)
    mi_y <- with(zmi,lmer(modelform[1]))
    mi_k <- with(zmi,lmer(modelform[2]))
  }
  #Prediction accuracy
  ori_predy <- predict(ori_y,newdata=datatest,allow.new.levels=TRUE)
  ori_predk <- predict(ori_k,newdata=datatest,allow.new.levels=TRUE)
  ds_reem_predy <- predict(ds_reem_y,newdata=datatest,allow.new.levels=TRUE)
  ds_reem_predk <- predict(ds_reem_k,newdata=datatest,allow.new.levels=TRUE)
  ds_lmm_predy <- predict(ds_lmm_y,newdata=datatest,allow.new.levels=TRUE)
  ds_lmm_predk <- predict(ds_lmm_k,newdata=datatest,allow.new.levels=TRUE)
  mi_predy <- rowMeans(sapply(mi_y$analyses,function(x) predict(x,newdata=datatest,allow.new.levels=TRUE)))
  mi_predk <- rowMeans(sapply(mi_k$analyses,function(x) predict(x,newdata=datatest,allow.new.levels=TRUE)))
  pred.acc <- data.frame(absdiff=c(abs(ori_predy-datatest$Y),
                          abs(ori_predk-datatest$K),
                          abs(ds_reem_predy-datatest$Y),
                          abs(ds_reem_predk-datatest$K),
                          abs(ds_lmm_predy-datatest$Y),
                          abs(ds_lmm_predk-datatest$K),
                          abs(mi_predy-datatest$Y),
                          abs(mi_predk-datatest$K)),
                         datatype=rep(c("original","DataSifterREEM","DataSifterLMM","MI"),each=2*nrow(datatest)),
                         variable=rep(rep(c("Y","K"),each=nrow(datatest)),4))

  if(dim(result)[1]!=0){
    result <- list(result,pred.acc)
  }else{result=pred.acc}
  return(result)
  }
