#' Artificial Missing Introduce and Imputation Step
#'
#' @param data A data frame contains original data to be processed. The data must be in long format. Missingness is allowed in time-varying varaibles.
#' @param mispct Percent of artificial missing that should be introduced for obfuscation. 20\%-30\% is recommended for utility preservation.
#' @param misw Type of sampling weights or missingness level. "peri" is to consider weights on subject level, which means any subjects with partial missing would be excluded from complete cases. "perij" is to consider weights on subject and time level. Only subjects with all time points missing would be excluded from complete cases.
#' @param lnames A vector of longitudinal variables names.
#' @param timevar The time variable or cluster varaible name.
#' @param ID Name of the ID variable in the dataset.
#' @param maxit Maximal iteration. The default is 10 times.
#' @param crit Critical value for the stopping criteria. The default is 0.05, which stops the algorithm when the absolute deviance of the imputed and original value is within 5\% of the original values.
#'
#' @import dplyr glmmLasso glmnet doParallel lme4 missForest REEMtree MASS nlme
#' @export


repSifterImp <- function(data,mispct,lnames,timevar,ID,maxit=10,crit=0.05){
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
      family <- checkLink(splitdata[[1]][[1]],lnames[1])
      temp <- filter.lambda(splitdata[[1]][[1]],lnames[1],ID,seq(500,0,-10),0,family)
      model <- temp[[1]]
      lambdaseq[[1]] <- temp[[2]]
    test.data<-data.frame(splitdata[[i]][[2]],splitdata[[i]][[3]])
    names(test.data)[1] = lnames[i]
    #Model prediction
      temp.pred<-predict(model,newdata=test.data,type="response")
      if(family$family=="binomial"){
        temp.pred <- sapply(1:length(temp.pred),function(x)rbinom(1,1,prob=temp.pred[x]))
      }else if(family$family!="gaussian"){
        temp.pred <- round(temp.pred,0)}
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
        if(iter==0){
          family=checkLink(splitdata[[i]][[1]],lnames[i])
          temp <- filter.lambda(splitdata[[i]][[1]],lnames[i],ID,seq(500,0,-50),50,family)
          model <- temp[[1]]
          print(model)
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

      test.data<-data.frame(splitdata[[i]][[2]],splitdata[[i]][[3]])
      names(test.data)[1] = lnames[i]

      #Model prediction
        temp.pred<-predict(model,newdata=test.data,type="response")
        #,drop.unused.levels = FALSE)
        if(family$family=="binomial"){
          temp.pred <- sapply(1:length(temp.pred),function(x)rbinom(1,1,prob=temp.pred[x]))
        }else if(family$family!="gaussian"){
          temp.pred <- round(temp.pred,0)}

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
