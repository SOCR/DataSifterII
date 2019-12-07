#' DataSifter II Algorithm for Time-varying Data
#'
#' Create a informative privacy-preserving time-varying dataset that guarantees subjects' privacy while preserving the information contained in the original dataset.
#'
#' @param data A data frame contains original data to be processed. The data must be in long format. Missingness is allowed in time-varying varaibles.
#' @param mispct Percent of artificial missing that should be introduced for obfuscation. 20\%-30\% is recommended for utility preservation.
#' @param misw Type of sampling weights or missingness level. "peri" is to consider weights on subject level, which means any subjects with partial missing would be excluded from complete cases. "perij" is to consider weights on subject and time level. Only subjects with all time points missing would be excluded from complete cases.
#' @param lnames A vector of longitudinal variables names.
#' @param timevar The time variable or cluster varaible name.
#' @param ID Name of the ID variable in the dataset.
#' @param maxit Maximal iteration. The default is 10 times.
#' @param crit Critical value for the stopping criteria. The default is 0.05, which stops the algorithm when the absolute deviance of the imputed and original value is within 5\% of the original values.
#' @param cal.weights.method Raw data missingness model for calculating IPW weights. If method = "param", the function utilize logistic regression ("peri") or GLMM ("perij") for missingness model. If method = "nonparam", the function utilize random forest ("peri") for missingness model.
#'
#' @return
#' \itemize{
#'   \item siftdata - Sifted data frame.
#'   \item data - Original data frame.
#' }
#'
#' @export
repSifter <- function(data,mispct,misw="perij",lnames,timevar,ID,maxit=10,crit=0.05,cal.weights.method="param"){
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
                           timevar = timevar,ID="ID",maxit = maxit,crit = crit)
  return(list(siftdata,data))
}
