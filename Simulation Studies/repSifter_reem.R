#repSifter on RE-EMtree
library(dplyr)
library(glmmLasso)
library(glmnet)
library(doParallel)
library(lme4)
library(missForest)
library(REEMtree)
library(MASS)
library(nlme)

#RE-EMtree
#function to calculate RE-EMtree model
RE_EMmodel<-function(data,outcome = "mis",ID=ID,random=NULL,minsplit=20,cp=0.01,xval=10){
  if(is.null(random)) random=formula(paste("~1|",ID,sep = ""))
  form <- as.formula(paste(outcome,"~",paste(names(data)[-which(names(data) %in% c(ID,outcome))],collapse = "+")))
  reem <- REEMtree::REEMtree(formula=form, data=data, random = random,
                                     tree.control = rpart.control(minsplit = minsplit,cp=cp,xval = xval))
  return(reem)
}

#Calculate sampling weights for first imputation - output weights and full observation rows for each varaible.
#Params:
#data = original data (numerical variables only) with missing values.
#lcolnames = a vector of longitudinal variables names.
#w = type of sampling weights. "peri" is to assign weights on subject level. "perij" is to assign weights on subject and time level.
#timevar = the time variable or cluster varaible name.
calWeights <- function(data,lcolnames,ID,w=c("peri","perij"),timevar,method="param",
                       random=NULL,minsplit=20,cp=0.01,xval=10){
  #Input inspection
  if(class(lcolnames)!="character") stop("Invalid longitudinal column names.")
  if(!w %in% c("peri","perij")) stop("Weights is eighter per person(“peri”) or per person time(“perij”).")
  nlcol <- length(lcolnames)
  wlist <- list(nlcol)
  lInd <- which(names(data) %in% lcolnames)
  idInd <- which(names(data) == ID)
  #Build models for missingness per variable
  if(w=="peri"){
    #limit the data into one patient 1 line.
    w_rowidx <- lapply(lcolnames, function(x) {
      misid <- data[which(is.na(data[,x])),ID]
      #datasub <- data[!data[,ID] %in% unique(misid),]
      datasub <- data %>% group_by(ID) %>% filter(row_number()==1) %>%
        mutate(mis=!(ID %in% unique(misid))) %>% ungroup()%>%
        dplyr::select(-c(lcolnames,ID,timevar))
      if(method=="param"){
      modelmis <- cv.glmnet(x=data.matrix(datasub[,-which(names(datasub)=="mis")]),
                            y=as.numeric(datasub$mis), alpha=1,family = "binomial")
      modelmis <- glmnet(x=data.matrix(datasub[,-which(names(datasub)=="mis")]),
                         y=as.numeric(datasub$mis),lambda = modelmis$lambda.min)
      #print(modelmis$beta)
      testdata<-as.matrix(data[!data[,ID] %in% unique(misid),]%>%
                            dplyr::select(names(data)[names(data) %in% names(datasub)]))
      weights <- predict(modelmis,newx=data.matrix(testdata) ,type = "response")
      } else {
      ########THERE ARE BUGS HERE!!!###########
      mtyinit <- floor(sqrt(ncol(datasub)))
      mtryseq <- seq(2:min(mtryint+20,ncol(datasub)))
      selectmtry <- sapply(mtryseq, function(x) randomForest(factor(mis)~.,data = datasub,ntree=1000,mtry=x)$err.rate)
      mtryopt <- which.min(selectmtry)
      modelmis <- randomForest(factor(mis)~.,data = datasub,ntree=1000,mtry=mtryopt)
      testdata<-as.matrix(data[!data[,ID] %in% unique(misid),]%>%
                            dplyr::select(names(data)[names(data) %in% names(datasub)]))
      weights <- predict(modelmis,newdata=testdata ,type = "response")
      }
      rownumb <- as.numeric(rownames(weights))
      list(1/weights,rownumb)
    })
  } else {
    w_rowidx <- lapply(lcolnames, function(x) {
      misid <- data[which(is.na(data[,x])),ID]
      misfreqtb <- data.frame(ID=names(table(misid)),freq=as.numeric(table(misid)))
      datasub <- data
      mis <- ifelse(is.na(datasub[,x]),0,1)
      datasub <- cbind(datasub,mis) %>% dplyr::select(-lcolnames)
      datasub[,ID] <- as.factor(as.character(datasub[,ID]))
      #Here ID must equal to "ID"
     if(method=="param"){
       modelmis <- filter.lambda(data=datasub,outcome = "mis",ID=ID,lambdaseq = seq(250,0,-50),family = binomial(link = "logit"),startPQL = FALSE)[[1]]
       weights <- predict(modelmis,newdata = datasub[mis !=0,],type = "response")
       }else{
       modelmis <- RE_EMmodel(data=datasub,outcome = "mis",ID=ID,random=random,
                              minsplit=minsplit,cp=cp,xval=xval)
       weights <- predict(modelmis,newdata = datasub[mis !=0,],
                          type = "response",id=datasub[mis !=0,ID],EstimateRandomEffects=TRUE)
     }
      #print(weights)
      rownumb <- which(datasub$mis !=0)
      list(1/weights,rownumb)
    })
  }
  return(w_rowidx)
}




repSifterImp <- function(data,mispct,lnames,timevar,ID,reem = FALSE,maxit=10,crit=0.05,minsplit=20,cp=0.01,xval=10){
  #input checks
  if(class(lnames)!="character"){stop("Invalid longitudinal variable names.")}
  if(class(timevar)!="character"){stop("Invalid time variable names.")}
  if(!"data.frame" %in% class(data)){stop("Invalid data input.")}
  if(class(ID)!="character"|length(ID)!=1){stop("Invalid ID variable. It should be a single variable name.")}
  if(mispct<0|mispct>1){stop("Invalid artificial missing percentage.")}
  if(mispct>0.5){warning("Missingness might be too large.")}
  #Note that the original data has no missing at all.
  names(data)[which(names(data)==ID)]<-"ID"
  data$ID<-as.factor(data$ID)
  ID<-"ID"
  #1. Assign missing to the data
  ln <- length(lnames)
  dmis<-data
  dmis[,lnames] <- prodNA(data.frame(dmis[,lnames]),noNA = mispct)
  #2. Record missing cells; extract complete data for each variable
  impCVCF <- CVCF(data = dmis, lnames = lnames,ID=ID, timevar = timevar)
  mislist <- list()
  splitdata <- list()
  for(i in 1:ln){
    #Record missing index
    mislist[[i]]<-which(is.na(dmis[,lnames[i]]))
    #Cut the data into missing and completedata to build data
    splitdata[[i]]<-list()
    splitdata[[i]][[1]]<-cbind(dmis[-mislist[[i]],lnames[i]],
                               impCVCF[-mislist[[i]],-which(names(dmis) %in% lnames[i])])#complete portion of data
    names(splitdata[[i]][[1]])[1]<-lnames[i]
    splitdata[[i]][[2]]<-data[mislist[[i]],lnames[i]]#True missing values
    splitdata[[i]][[3]]<-impCVCF[mislist[[i]],-which(names(dmis) %in% lnames[i])]#Covariates for missing
    splitdata[[i]][[4]]<-impCVCF[mislist[[i]],lnames[i]]#CVCF first imputed value
  }
  #Vector of longitudinal variables that needs further imputation  
  vec <- 1:ln 
  ind <- 1:nrow(data)
  iter <- 0
  lambdaseq <- list()
  
  ####Spectial case when only one longitudinal variable is present.
  if(length(vec)==1){
    target <- splitdata[[1]][[2]]
    tempimp <- splitdata[[1]][[4]]
    if(reem==TRUE){
      model <- RE_EMmodel(data = splitdata[[1]][[1]],
                          outcome = lnames[1],ID = ID,cp=cp,minsplit = minsplit,xval = xval)
    }else{
        family <- checkLink(splitdata[[1]][[1]],lnames[1])
        temp <- filter.lambda(splitdata[[1]][[1]],lnames[1],ID,seq(250,0,-50),0,family)
        model <- temp[[1]]
        lambdaseq[[1]] <- temp[[2]]
      } 
    test.data<-data.frame(splitdata[[i]][[2]],splitdata[[i]][[3]])
    names(test.data)[1] = lnames[i]
    #Model prediction
    if(reem==TRUE){
      temp.pred<-predict(model,newdata=test.data,id=test.data$ID,EstimateRandomEffects=TRUE)
    }else{  
      temp.pred<-predict(model,newdata=test.data,type="response")
      if(family$family=="binomial"){
        temp.pred <- sapply(1:length(temp.pred),function(x)rbinom(1,1,prob=temp.pred[x]))
      }else if(family$family!="gaussian"){
        temp.pred <- round(temp.pred,0)}
    }
    ##When predict please make sure that the ID variable is converted.
    splitdata[[1]][[4]]<-temp.pred
    dmis[mislist[[1]],lnames[1]] <- splitdata[[1]][[4]]
  }
  
  ####Case when more than 1 longitudinal variables are present.
  while(length(vec)>1){
    for(i in vec){
      #initialize values
      target <- splitdata[[i]][[2]]
      tempimp <- splitdata[[i]][[4]]
      #Initiate lambdas
      if(iter>maxit|diff.imp(target,tempimp)<crit) {
        vec <- vec[-which(vec==i)]
        next
      }
      if(reem==TRUE){
      model <- RE_EMmodel(data = splitdata[[i]][[1]],outcome = lnames[i],ID = ID,
                          cp=cp,minsplit = minsplit,xval = xval)
      } else {
      if(iter==0){
        family=checkLink(splitdata[[i]][[1]],lnames[i])
        temp <- filter.lambda(splitdata[[i]][[1]],lnames[i],ID,seq(250,0,-50),50,family)
        model <- temp[[1]]
        #print(model)
        lambdaseq[[i]] <- temp[[2]]
      } else {
        if(lambdaseq[[i]][1]!=-1){
        family=checkLink(splitdata[[i]][[1]],lnames[i])
        temp <- filter.lambda(splitdata[[i]][[1]],lnames[i],ID,lambdaseq[[i]],family=family)
        model <- temp[[1]]
        }else{
        preds<-names(splitdata[[i]][[1]])[-which(names(splitdata[[i]][[1]]) %in% c(lnames[i],"ID"))]
        preds<-paste(preds,collapse =  "+")
        form <- as.formula(paste(lnames[i],"~",preds,sep = ""))  
        model <-  glmmPQL(form,random = ~1|ID,data=splitdata[[i]][[1]],family = family) 
        }
      }
      }
      test.data<-data.frame(splitdata[[i]][[2]],splitdata[[i]][[3]])
      names(test.data)[1] = lnames[i]
      
      #Model prediction
      if(reem==TRUE){
        temp.pred<-predict(model,newdata=test.data,id=test.data$ID,EstimateRandomEffects=TRUE)
        }else{  
        temp.pred<-predict(model,newdata=test.data,type="response")
                           #,drop.unused.levels = FALSE)
        if(family$family=="binomial"){
        temp.pred <- sapply(1:length(temp.pred),function(x)rbinom(1,1,prob=temp.pred[x]))
      }else if(family$family!="gaussian"){
        temp.pred <- round(temp.pred,0)}
      }
      ##When predict please make sure that the ID variable is converted.
      splitdata[[i]][[4]]<-temp.pred
      #update all other's missingness
      fullvec <- dmis[,lnames[i]]
      fullvec[mislist[[i]]] <- temp.pred
      vectemp<-vec[which(!vec %in% i)]
      if(!is.null(vectemp)){
        for(j in vectemp) {
          splitdata[[j]][[1]][,lnames[i]] <- fullvec[!ind %in% mislist[[j]]]
          splitdata[[j]][[3]][,lnames[i]] <- fullvec[ind %in% mislist[[j]]]
        }#Update all the covariates
      }
    }
    iter=iter+1
    for(k in 1:ln) {dmis[mislist[[k]],lnames[k]] <- splitdata[[k]][[4]]}#After one round of iteration for all variables update all the imputed values
  }
  return(dmis)  
} 

# calWeights <- function(data,lcolnames,ID,w=c("peri","perij"),timevar){
#   #Input inspection
#   if(class(lcolnames)!="character") stop("Invalid longitudinal column names.")
#   if(!w %in% c("peri","perij")) stop("Weights is eighter per person(peri) or per person time(perij).")
#   nlcol <- length(lcolnames)
#   wlist <- list(nlcol)
#   lInd <- which(names(data) %in% lcolnames)
#   idInd <- which(names(data) == ID)
#   #Build models for missingness per variable
#   if(w=="peri"){
#     #limit the data into one patient 1 line.
#     w_rowidx <- lapply(lcolnames, function(x) {
#       misid <- data[which(is.na(data[,x])),ID]
#       #datasub <- data[!data[,ID] %in% unique(misid),]
#       datasub <- data %>% group_by(ID) %>% filter(row_number()==1) %>%
#         mutate(mis=!(ID %in% unique(misid))) %>% ungroup()%>%
#         dplyr::select(-c(lcolnames,ID,timevar))
#       modelmis <- cv.glmnet(x=data.matrix(datasub[,-which(names(datasub)=="mis")]),
#                             y=as.numeric(datasub$mis), alpha=1,family = "binomial")
#       #modelmis <- cv.glmnet(x=as.matrix(datasub[,-which(names(datasub)=="mis")]),
#       #                      y=as.numeric(datasub$mis), alpha=1,family = "binomial")
#       #print(as.matrix(coef(modelmis)))
#       modelmis <- glmnet(x=data.matrix(datasub[,-which(names(datasub)=="mis")]),
#                          y=as.numeric(datasub$mis),lambda = modelmis$lambda.min)
#       #modelmis <- glmnet(x=as.matrix(datasub[,-which(names(datasub)=="mis")]),
#       #                   y=as.numeric(datasub$mis),lambda = modelmis$lambda.min)
# 
#       #print(modelmis$beta)
#       testdata<-as.matrix(data[!data[,ID] %in% unique(misid),]%>%
#                             dplyr::select(names(data)[names(data) %in% names(datasub)]))
#       weights <- predict(modelmis,newx=data.matrix(testdata) ,type = "response")
#       ##########wrong dimention for newx!!!!!!!!!!!!!
#       rownumb <- as.numeric(rownames(weights))
#       list(1/weights,rownumb)
#     })
#   } else {
#     w_rowidx <- lapply(lcolnames, function(x) {
#       misid <- data[which(is.na(data[,x])),ID]
#       misfreqtb <- data.frame(ID=names(table(misid)),freq=as.numeric(table(misid)))
#       datasub <- data
#       mis <- ifelse(is.na(datasub[,x]),0,1)
#       datasub <- cbind(datasub,mis) %>% dplyr::select(-lcolnames)
#       datasub[,ID] <- as.factor(as.character(datasub[,ID]))
#       #Here ID must equal to "ID"
#       modelmis <- filter.lambda(data=datasub,outcome = "mis",ID=ID,lambdaseq = seq(0,500,5),family = binomial(link = "logit"),startPQL = FALSE)[[1]]
#       weights <- predict(modelmis,newdata = datasub[mis !=0,],type = "response")
#       #print(weights)
#       rownumb <- which(datasub$mis !=0)
#       list(1/weights,rownumb)
#     })
#   }
#   return(w_rowidx)
# }


#Weighted GlmmLASSO - output glmm model for imputation.
####
#Params:
#x = data frame or matrix of candidate predictors cross-sectional and longitudinal from the full data.
#y = vector of longitudinal outcome of interest.
#weights = sampling weights
#rowid = observed row id this will match with the length of weights.
#lvars = Longitudinal variable names.
#family = family object to specify the link function for outcome.
#withlong = a logical parameter to decide if the weighted linear mixed model depends on longitudinal predictors.

wlmmLasso <- function(x,y,yname,weights,rowid,ID="ID",lvars,
                      family=gaussian(link = "identity"),withlong=FALSE){
  #Input validation testa
  if(!class(x) %in% c("data.frame","matrix")) stop("Invalid x input.")
  if(!class(y) %in% c("vector","data.frame","matrix","numeric","character","integer")) stop("Invalid y input.")
  if(length(weights) != length(rowid)) stop("Invalid weights and observed data selection.")
  
  #Variable selection with glm(can only select static variables.)
  #if(class(y) %in% c("data.frame","matrix")){y=c(y)}
  xname<-colnames(x)
  if("data.frame" %in% class(y)) y <- as.numeric(y[,1])
  
  if(withlong==FALSE){
    delnames <- names(x)[names(x) %in% lvars]
    train.data<-cbind.data.frame(x[rowid,],y[rowid])
    colnames(train.data)=c(xname,yname)
    x <- x %>% dplyr::select(-delnames)
  } else {
    train.data<-cbind.data.frame(x[rowid,],y[rowid])
    colnames(train.data)=c(xname,yname)
  }
  
  x.train <- x[rowid,]
  x.test <- x[-rowid,]
  y.train <- y[rowid]
  y.test <- y[-rowid]
  varselect <- cv.glmnet(x=as.matrix(x.train %>% dplyr::select(-ID)),y=y.train,weights = weights, alpha=1)
  varselect <- glmnet(x=as.matrix(x.train %>% dplyr::select(-ID)),
                      y=y.train,weights = weights, alpha=1,lambda = varselect$lambda.min)
  coefs = as.matrix(coef(varselect)) # convert to a matrix (p by 1)
  #print(coefs)
  ix = which(abs(coefs[,1]) > 0)
  svars<-row.names(coefs)[ix]
  if(svars[1]== "(Intercept)") svars<-svars[-1]
  #Fit GLMM with the selected variables and weights.
  glmmform <- as.formula(paste(yname,"~",paste(svars,collapse = "+"),paste("+(1|",ID,")",sep = "")))
  glmmmodel <- glmer(glmmform,data = train.data,weights = weights,family = family)
  return(glmmmodel)
}

#Closest Value Carry Forward function for first guessing.
####
# data = data for the closest or mean value imputataion.
# lnames = a vector of the longitudinal varaible names.
# timevar = name of the time variable. 
CVCF<-function(data, lnames, timevar, ID){
  #data=dataname, lnames=longitudinal variable names, timevar= time variable name
  #Check if the inputs are valid?
  if(!"data.frame"%in% class(data)){stop("Invalid data input.")}
  if(class(lnames)!="character"){stop("Invalid longitudinal variable names.")}
  if(class(timevar)!="character"){stop("Invalid time variable names.")}
  #Calculate the mean of each time point
  ln <- length(lnames)
  time_means <- data %>% group_by(get(timevar)) %>%
    summarize_at(lnames,funs(mean(.,na.rm = TRUE)))
  names(time_means)[1]<-timevar
  #Obtain the missingness for each id and each timepoint each varaible.
  #When categorical is observed the rowSums would be non NA all the time.
  misind <- data %>% dplyr::select(ID,timevar,lnames) %>% mutate(rsum=rowSums(.[-c(1,2)])) %>% 
    filter(is.na(rsum)) %>% dplyr::select(ID,lnames,timevar)
  misind$misvar <- apply(misind,1,function(x) paste(names(misind)[which(is.na(x))],collapse = ","))
  workdata <- data %>% dplyr::select(ID,lnames,timevar) %>% 
    filter(ID %in% misind[,1]) %>% arrange(ID,get(timevar))
  #Check if the misvar is which one or multiple vectors
  #Slow version
  for(i in 1:nrow(misind)){
    misvarvec <- unlist(strsplit(misind$misvar[i],","))
    tempdata <- workdata %>% filter(ID %in% misind[i,1])
    tempdatab <- tempdata %>% filter(get(timevar)<misind[i,timevar])
    tempdataa <- tempdata %>% filter(get(timevar)>misind[i,timevar])
    for(j in misvarvec){
      if(dim(tempdatab %>% filter(!is.na(get(j))))[1]>0){
        tempdatab_sub <- tempdatab %>% filter(!is.na(get(j)))
        misind[i,which(names(misind)%in% j)] <- tempdatab_sub[nrow(tempdatab_sub),j]
      }else if(dim(tempdataa %>% filter(!is.na(get(j))))[1]>0){
        tempdataa_sub <- tempdataa %>% filter(!is.na(get(j)))
        misind[i,which(names(misind)%in% j)] <- tempdataa_sub[1,j]    
      } else{
        misind[i,which(names(misind)%in% j)] <- as.numeric(time_means %>% filter(get(timevar)==misind[i,timevar]) %>% dplyr::select(j))    
      }
    }
  }
  for(i in 1:nrow(misind)){
    var <- unlist(strsplit(as.character(misind[i,"misvar"]),","))
    varind <- which(names(data) %in% var)
    id <- as.character(misind[i,ID])
    t <- misind[i,timevar]
    data[which(as.vector(data[,ID]) == c(id) & as.vector(data[,timevar]) == c(t)),var] <- misind[i,var]
  }
  return(data)
}

#Stopping criteria: Function for calculating the difference between imputed and original values:
####
#ori = original numerical vector
#imp = imputed numerical vector
diff.imp<-function(ori,imp){
  return(sum((imp-ori)^2)/sum(ori^2))
}


#The bic of the GLMM lasso model, given specific lambda. 
####
# lambdainp = lambda parameter for the GLMM lasso model
# data = data frame with complete longitudinal and cross-sectional data.
# family = family object to specify the link function for outcome.
try_glmmlasso<-function(lambdainp,data,family,form,Delta.start=NULL,q.start=NULL){
  d = nrow(data)+1
  if(is.null(Delta.start)&is.null(q.start)){
    Delta.start<-as.matrix(t(rep(0,d)))
    q.start<-0.1
}
  glm1 <- try(glmmLasso(form, rnd = list(ID=~1),  
                        family = family, data = data, lambda=lambdainp,switch.NR=FALSE,final.re=TRUE,
                        control=list(print.iter=FALSE,start=Delta.start,q.start=q.start)), silent=TRUE)  
  bic<-ifelse(class(glm1)!="try-error",glm1$bic,Inf)
  if(class(glm1)!="try-error"){
    Delta.start <- glm1$Deltamatrix[glm1$conv.step,]
    q.start<-glm1$Q_long[[glm1$conv.step+1]]
  }
  return(bic)
}

filter.lambda<-function(data,outcome,ID,lambdaseq,range=0,family=gaussian(link = "identity"),startPQL=TRUE){
  if(min(lambdaseq)<0) stop("Lambda must be greater than 0.")
  #Delete sigular values
  sname <- names(data)[which(sapply(data,function(x) length(unique(x))==1))]
  #Define formula
  preds<-names(data)[-which(names(data) %in% c(outcome,"ID",sname))]
  preds<-paste(preds,collapse =  "+")
  form <- as.formula(paste(outcome,"~",preds,sep = ""))
  #If the glmmLasso wouldnt converge--> we dont consider use lasso to penalize.
  BIC_vec<-c(Inf)
  #Starting values
  d = length(strsplit(as.character(form),"\\+")[[3]])#Length of the varaibles needed.
  intercept_model <- glmmPQL(as.formula(paste(gsub("^~.$","",form)[2],"~1")),
                               random = ~1|ID,data=data,family = family) 
  Delta.start<-c(as.numeric(intercept_model$coef$fixed),rep(0,d),as.numeric(t(intercept_model$coef$random$ID)))
  q.start<-as.numeric(VarCorr(intercept_model)[1,1])
  #Parallel process
  ncl <- detectCores()
  if(length(lambdaseq)<ncl){ncl=length(lambdaseq)}

  cl <- parallel::makeCluster(ncl,setup_strategy = "sequential")
  registerDoParallel(cl)
  clusterEvalQ(cl, { library(parallel)
                     library(doParallel)
                     library(glmmLasso)})
  clusterExport(cl, list("try_glmmlasso"),envir = .GlobalEnv)
  clusterExport(cl, list("lambdaseq","data","family"),envir = environment())
  if(startPQL==TRUE){
  BIC_vec<-foreach(k = 1:length(lambdaseq),.combine = 'rbind') %dopar%{
    try_glmmlasso(lambdainp=lambdaseq[k],data=data,family = family,form = form,
                  Delta.start = Delta.start,q.start = q.start)
  #BIC_vec<-do.call(rbind,mclapply(1:length(lambdaseq),
  #    function(k)try_glmmlasso(lambdainp=lambdaseq[k],data=data,family = family,
  #                             form = form,Delta.start = Delta.start,q.start = q.start),mi.cores=ncl))
      
  }}

  if(min(BIC_vec)==Inf){
    d = nrow(data)+1
    Delta.start<-as.matrix(t(rep(0,d)))
    q.start<-0.1
    BIC_vec<-foreach(k = 1:length(lambdaseq),.combine = 'rbind') %dopar%{
      try_glmmlasso(lambdainp=lambdaseq[k],data=data,family = family,form = form,
                    Delta.start = Delta.start,q.start = q.start)
    }
    #BIC_vec<-do.call(rbind,mclapply(1:length(lambdaseq),
    #function(k)try_glmmlasso(lambdainp=lambdaseq[k],data=data,family = family,
    #                         form = form,Delta.start = Delta.start,q.start = q.start),mi.cores=ncl))
    Delta.start = Delta.start
    q.start = q.start
  }
  parallel::stopCluster(cl)
  if(min(BIC_vec)==Inf){
    glm1 <- glmmPQL(form,random = ~1|ID,data=data,family = family)   
    seqs <- -1
  }else{
  opt<-which.min(BIC_vec)
  blambda<-lambdaseq[opt]
  blambda<-blambda[1]
  print(blambda)
  glm1 <- try(glmmLasso(form, rnd = list(ID=~1),
                        family = family, data = data, lambda=lambdaseq[opt],switch.NR=FALSE,final.re=TRUE,
                        control=list(print.iter=FALSE,start=Delta.start,q.start=q.start)), silent=TRUE)
  if(range==0){ 
    seqs=blambda
  } else if(is.na(blambda)){
    stop("No solution")
    }else {seqs=seq(max(blambda-range,0),blambda+range,1)
    }
  }
  return(list(glm1,seqs))
}


#Define link function for GLMM
####
# data = data frame object with outcome variable.
# name = the name of the outcome varaible.
checkLink<-function(data,name){#Problem with age given log link
  outcome<-data[,name]
  if(length(levels(factor(outcome)))==2){
    link <- binomial(link = "logit")
  }else if(class(outcome)=="factor" & length(levels(outcome))>2 & length(levels(outcome))<20){
    link <- binomial(link = "logit")
  #} else if(class(outcome)%in%c("numeric","integer") & sum(outcome%%1)==0&min(outcome)>=0){
  #  link <- poisson(link = "log")
  }else link <- gaussian(link = "identity")
  return(link)
}

####Main Function####

repSifter <- function(data,mispct,misw="perij",lnames,timevar,ID,reem = FALSE,maxit=10,crit=0.05,minsplit=20,cp=0.01,xval=10,cal.weights.method="param"){
  #input checks
  if(class(lnames)!="character"){stop("Invalid longitudinal variable names.")}
  if(class(timevar)!="character"){stop("Invalid time variable names.")}
  if(!"data.frame" %in% class(data)){stop("Invalid data input.")}
  if(class(ID)!="character"|length(ID)!=1){stop("Invalid ID variable. It should be a single variable name.")}
  if(mispct<0|mispct>1){stop("Invalid artificial missing percentage.")}
  if(mispct>0.5){warning("Missingness might be too large.")}
  #Change the ID name to "ID" to avoid problems.
  names(data)[which(names(data)==ID)]<-"ID"
  data$ID<-as.factor(data$ID)
  ID<-"ID"
  #Check missing
  orimis <- data[,lnames]
  orimis <- sapply(orimis, function(x) is.na(sum(as.numeric(x))))
  if(length(orimis[which(orimis==1)])>0){
    orimis <- names(orimis[which(orimis==1)])
    #First imputation for missing repeated measures
    lw <- calWeights(data = data,lcolnames = orimis,ID = ID,w = misw,timevar = timevar,method=cal.weights.method)
    ln <- length(lw)
    
    impCVCF <- CVCF(data = data, lnames = orimis,ID=ID, timevar = timevar) #mean imputation
    obslist <- list()
    mislist <- list()
    splitdata <- list()
    for(i in 1:ln){
      #Record missing index
      obslist[[i]]<-lw[[i]][[2]]
      mislist[[i]]<-which(is.na(data[,orimis[i]]))
      #Cut the data into missing and completedata to build data
      splitdata[[i]]<-list()
      splitdata[[i]][[1]]<-cbind(data[obslist[[i]],orimis[i]],
                                 impCVCF[obslist[[i]],-which(names(data) %in% orimis[i])])
      #complete portion of data
      names(splitdata[[i]][[1]])[1]<-orimis[i]
      splitdata[[i]][[2]]<-data[mislist[[i]],orimis[i]]#totally NA
      splitdata[[i]][[3]]<-impCVCF[mislist[[i]],-which(names(data) %in% orimis[i])]#Covariates for missing
      splitdata[[i]][[4]]<-impCVCF[mislist[[i]],orimis[i]]#CVCF first imputed value
      splitdata[[i]][[5]]<-rep(1,dim(splitdata[[i]][[1]])[1]) #Imputed observd value
    }
    #vec Vector of longitudinal variables that needs further imputation  
    vec <- 1:ln 
    ind <- 1:nrow(data)
    iter <- 0
    lambdaseq <- list()
    while(length(vec)>1){
      for(i in vec){
        #initialize values
        target <- splitdata[[i]][[1]][,orimis[1]] # observed values
        tempimp <- splitdata[[i]][[5]] # imputed observed values
        #Initiate lambdas
        if(iter>maxit|diff.imp(target,tempimp)<crit) {
          vec <- vec[-which(vec==i)]
          next
        }
        family=checkLink(splitdata[[i]][[1]],orimis[i])
        #print(family)
        #weighted GLMM model
        model <- wlmmLasso(x=impCVCF[,-which(names(data) %in% orimis[i])],
                           y=impCVCF[,orimis[i]],
                           weights = lw[[i]][[1]],
                           rowid = lw[[i]][[2]],
                           lvars = orimis,
                           ID = "ID",
                           family = family,
                           withlong=TRUE,
                           yname=orimis[i])
        splitdata[[i]][[5]] <- fitted(model)  
        test.data <- data.frame(splitdata[[i]][[2]],splitdata[[i]][[3]])
        names(test.data)[1] = orimis[i]
        temp.pred<-predict(model,newdata=test.data,type="response",allow.new.levels = TRUE)
        if(family$family=="binomial"){
          temp.pred <- sapply(1:length(temp.pred),function(x)rbinom(1,1,prob=temp.pred[x]))
        }else if(family$family!="gaussian"){
          temp.pred <- round(temp.pred,0)}
        ##When predict please make sure that the ID variable is converted.
        splitdata[[i]][[4]]<-temp.pred
        #update all other's missingness
        fullvec <- data[,orimis[i]]
        fullvec[mislist[[i]]] <- temp.pred
        vectemp<-vec[which(!vec %in% i)]
        if(!is.null(vectemp)){
          for(j in vectemp) {
            splitdata[[j]][[1]][,orimis[i]] <- fullvec[!ind %in% mislist[[j]]]
            splitdata[[j]][[3]][,orimis[i]] <- fullvec[ind %in% mislist[[j]]]
          }#Update all the covariates
        }
      }
      iter=iter+1
      for(k in 1:ln) {data[mislist[[k]],orimis[k]] <- splitdata[[k]][[4]]}
      #After one round of iteration for all variables update all the imputed values
    }
  }
  siftdata <- repSifterImp(data=data,mispct = mispct,lnames = lnames,
                           timevar = timevar,ID="ID",reem=reem,maxit = maxit,crit = crit,minsplit = minsplit,cp = cp,xval = xval)
  return(list(siftdata,data))
}


#####Check params###################################
#sdata = Sifted data frame.
#oridata = original data frame.
#covtrue = ture coefficients of the generative model. 
#form = true generative model formula according to "covtrue" order.
#family = family object to specify the link function for the outcome.
betarate_g<-function(sdata,oridata,covtrue=c(10,0.2,5,1,2),form=as.formula("Y~x1+x2+x3+visit+(1|ID)"),family=gaussian(link = "identity")){
  if(length(covtrue)!=length(unlist(strsplit(gsub(" ","",as.character(form)[3]),"[+]")))) warning("Different length between true covariates and formula.")
  #Beta rates
  covrate<-data.frame()
  m0 <- glmer(form,data=oridata,family=family)
  m1 <- glmer(form,data=sdata,family=family)
  modelori<-summary(m0)
  modelsif<-summary(m1)
  covrate<-abs(modelori$coefficients[,1]-modelsif$coefficients[,1])/abs(modelori$coefficients[,1])
  covrate_t<-abs(covtrue-modelsif$coefficients[,1])/abs(covtrue)
  covrate_ori_t<-abs(covtrue-modelori$coefficients[,1])/abs(covtrue)
  #Normal Confidence Interval
  c0 <- confint(m0,parm="beta_",level = 0.95)
  c0 <- c0[complete.cases(c0),]
  
  c1 <- confint(m1,parm="beta_",level = 0.95)
  c1 <- c1[complete.cases(c1),]
  #Bootstrap CI
  #c0 <- bootci(data=oridata,form = form)
  #c1 <- bootci(data=sdata,form = form)
  if(dim(c0)[1]!=dim(c1)[1]) warning("Dimension of CIs do not match!")
  overlap <- sapply(1:dim(c0)[1], function(x) 1-(c1[x,1]>c0[x,2]|c0[x,1]>c1[x,2]))
  trueinci <- sapply(1:dim(c1)[1], function(x){ (c1[x,1]<=covtrue[x]&c1[x,2]>=covtrue[x])})
  #print(trueinci)
  trueinorici <- sapply(1:dim(c0)[1], function(x) (c0[x,1]<=covtrue[x]&c0[x,2]>=covtrue[x]))
  #print(trueinorici)
  return(list(covrate,covrate_t,covrate_ori_t,overlap,trueinci,c1,trueinorici)) #returns the covariance rate compared to original data[[1]] and the trueth[[2]]
}


###########Calculate PIFV #######
#originaldata = original data frame
#sifteddata = sifted data frame
#casenum1 = row number in original data.
#casenum2 = row number is sifted data.
pctfeatures_match<-function(originaldata,sifteddata,casenum1,casenum2){
  pctmatch<-length(which(originaldata[casenum1,]==sifteddata[casenum2,]))/(ncol(originaldata)-length(which(is.na(originaldata[casenum1,]))))
  return(pctmatch)
}

#Calculate PIFV for entire dataset.
#originaldata = original data frame
#sifteddata = sifted data frame
pifv <- function(originaldata,sifteddata){
  if(nrow(originaldata)!=nrow(sifteddata)) stop("Dimensions does not match.")
  PIFV<-sapply(1:nrow(originaldata), function(x)
    pctfeatures_match(originaldata,sifteddata,x,x))
  return(PIFV)
}


##########Construct dummies for factor cross-sectional data########
create.dummy <- function(data,cs.names=NULL){
  names.list <- names(data)
  if(is.null(cs.names)) cs.names <- names(data)
  datacs <- data %>% dplyr::select(all_of(cs.names))
  datalg <- data %>% dplyr::select(-all_of(cs.names))
  fac_vars <- which(sapply(datacs,class)=="factor")
  char_vars <- which(sapply(datacs,class)=="character")
  if(length(c(fac_vars,char_vars))==0){stop("No factor or characters.")}
  if(length(char_vars>0)){
  datacs[,char_vars] <- sapply(datacs[,char_vars],function(x) as.factor(x))}
  fac_mtx<-model.matrix(as.formula(paste("~",paste(names(datacs)[c(fac_vars,char_vars)],collapse = "+"))),data=datacs)
  fac_mtx <- fac_mtx[,-1]#Delete intercept
  if(length(which(colSums(fac_mtx)==0))>0){
    fac_mtx <- fac_mtx[,-which(colSums(fac_mtx)==0)]}
  datacs <- datacs[,-fac_vars]
  datacs <- data.frame(datacs,fac_mtx)
  data <- data.frame(datalg,datacs)
  return(data)
}

###########Create synthetic data##############
#Long to wide multivariable format
myspread <- function(df, key, value) {
  # quote key
  keyq <- rlang::enquo(key)
  # break value vector into quotes
  valueq <- rlang::enquo(value)
  s <- rlang::quos(!!valueq)
  df %>% gather(variable, value, !!!s) %>%
    unite(temp, !!keyq, variable) %>%
    spread(temp, value)
}

#Variable range calculation.
neighbor.range.sift <- function(sortedlist, varname,N,sdvec){
  varlist <- sortedlist[[which(names(sortedlist)==varname)]]
  nonnarows <- which(!is.na(varlist[,2]))
  sdvar <- sdvec[varname]
  sift.var <- sapply(nonnarows, function(rownum){
    cur.rank <- varlist[rownum,1]
    cur.num <- varlist[rownum,2]
    range_neighbor <- range(varlist[which(varlist[,1] %in% c((cur.rank-N):(cur.rank+N))),2],na.rm = TRUE) 
    if(sum(is.na(range_neighbor))==1 & sdvar<(range_neighbor[2]-range_neighbor[1])){
      sifted_num <- runif(1,cur.num-0.5*sdvar,cur.num+0.5*sdvar)
    }else{
      sifted_num <- runif(1,range_neighbor[1],range_neighbor[2])
    }
    sifted_num
  })
  sifted_var <- rep(NA,dim(varlist)[1])
  sifted_var[nonnarows] <-sift.var
  return(sifted_var)
}

#Function to create synthetic dataset.
synthetic <- function(data,lnames,timevar,ID,N){#Data must be in long format
  data.wide <- data %>% myspread(timevar,lnames) #ok
  wide_longivars <- grep(paste(paste("_",lnames,sep = ""),collapse = "|"),names(data.wide),value = TRUE)
  sub.data <- data.wide %>% dplyr::select(wide_longivars)
  sortedlist <- lapply(sub.data,function(x) cbind(rank(x,na.last = FALSE,ties.method = "first"),x))
  sdvec <- sapply(sortedlist,function(x) sd(x[,2],na.rm = TRUE))
  
  range <- ifelse(N%%2==0,N/2,(N-1)/2)
  sift.data.syn <-sapply(wide_longivars, function(x)
    neighbor.range.sift(sortedlist=sortedlist,varname = x,N=range,sdvec=sdvec[x]))
  
  data.wide[,wide_longivars]<-sift.data.syn
  
  data.long <- data.wide %>% 
    gather(v, value, wide_longivars) %>% 
    separate(v, c(timevar, "col")) %>% 
    arrange(ID) %>% 
    spread(col, value)
  data.long <- data.long %>% dplyr::select(names(data))
  
  return(list(data.long,data))
}



