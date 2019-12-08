#' Calculating the BIC of a GLMM lasso model, given specific lambda.
#' @param lambdainp Lambda parameter for the GLMM lasso model.
#' @param data A data frame with complete longitudinal and cross-sectional data.
#' @param family A family object to specify the link function for outcome.
#' @param form Formula for GLMM lasso model.
#' @param Delta.start A scalar specifying starting value for GLMM lasso fixed effects parameters.
#' @param q.start A scalar or matrix of suitable dimension, specifying starting values for the random-effects variance-covariance matrix.
#' @import dplyr glmmLasso glmnet doParallel lme4 missForest REEMtree MASS nlme
#' @export
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
  # if(class(glm1)!="try-error"){
  #   Delta.start <- glm1$Deltamatrix[glm1$conv.step,]
  #   q.start<-glm1$Q_long[[glm1$conv.step+1]]
  # }
  return(bic)
}
