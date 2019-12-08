#' Optimal lambda and corresponding final GLMM LASSO model
#'
#' Select the optimal lambda that minimizes BIC. The starting values of the parameters play an essential part in GLMM LASSO model fitting.
#' We first consider glmmPQL fitted values. Then, we consider 0s for sparse signal scenario.
#'
#' @param data A long format complete data frame.
#' @param outcome Outcome variable of interest.
#' @param ID Name of the ID variable in the dataset.
#' @param lambdaseq A vector of lambdas to consider.
#' @param range The next lambda sequence to consider, given the best lambda from current lambdaseq.
#' @param family Link function for the outcome.
#' @param startPQL Logical. If TRUE initiate parameter estimates are generated from glmmPQL. IF FALSE, initiate parameter estimates are 0s.
#'
#' @return
#' \itemize{
#'    \item glm1 - The optimal GLMM under LASSO regularization with smallest BIC.
#'    \item seqs - The next Lambda sequence to be considered, when further search is needed.
#' }
#' @import dplyr glmmLasso glmnet doParallel lme4 missForest REEMtree MASS nlme
#' @export
#'
filter_lambda<-function(data,outcome,ID,lambdaseq,range=0,family=gaussian(link = "identity"),startPQL=TRUE){
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
  ncl <- detectCores()-1
  if(length(lambdaseq)<ncl){ncl=length(lambdaseq)}

  cl <- parallel::makeCluster(ncl)
  registerDoParallel(cl)
  clusterEvalQ(cl, c(library(parallel), library(doParallel), library(glmmLasso))  )
  clusterExport(cl, list("try_glmmlasso"),envir = .GlobalEnv)
  clusterExport(cl, list("lambdaseq","data","family"),envir = environment())
  if(startPQL==TRUE){
    BIC_vec<-foreach(k = 1:length(lambdaseq),.combine = 'rbind') %dopar%{
      try_glmmlasso(lambdainp=lambdaseq[k],data=data,family = family,form = form,
                    Delta.start = Delta.start,q.start = q.start)
    }}

  if(min(BIC_vec)==Inf){
    d = nrow(data)+1
    Delta.start<-as.matrix(t(rep(0,d)))
    q.start<-0.1
    BIC_vec<-foreach(k = 1:length(lambdaseq),.combine = 'rbind') %dopar%{
      try_glmmlasso(lambdainp=lambdaseq[k],data=data,family = family,form = form,
                    Delta.start = Delta.start,q.start = q.start)
    }
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
