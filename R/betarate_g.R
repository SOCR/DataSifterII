#' Sifted Data Utility Check
#'
#' Check the deviance of sifted and original data parameter estimates under a parametric model. Model specification is needed for the check.
#'
#' @param sdata The Sifted data frame.
#' @param oridata The original data frame.
#' @param covtrue Ture or estimated coefficients of the target model.
#' @param form Model formula according to "covtrue" order.
#' @param family Family object to specify the link function for the outcome.
#'
#' @return
#' \itemize{
#'    \item covrate - Absolute deviance between the sifted and original coefficients over the absolute value of original coefficients.
#'    \item covrate_t - Absolute deviance between the sifted and true coefficients over the absolute value of true coefficients.
#'    \item covrate_ori_t - Absolute deviance between the original and true coefficients over the absolute value of true coefficients.
#'    \item overlap - If the sifted and original CIs overlap.
#'    \item trueinci - If the true coefficients are in the sifted CIs.
#'    \item c1 - Confidence Intervals for the Sifted model.
#'    \item trueinorici - If the true coefficients are in the original CIs.
#' }
#'
#' @import dplyr glmmLasso glmnet doParallel lme4 missForest REEMtree MASS nlme
#' @export
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
