#' Weighted GlmmLASSO
#'
#' GLMM LASSO with weights for each record. This can be used when we have a biased sample.
#'
#' @param x Data frame or matrix of candidate predictors time-invariant and longitudinal from the full data.
#' @param y Vector of longitudinal outcome of interest.
#' @param weights Sampling weights.
#' @param rowid Observed row id this will match with the length of weights.
#' @param lvars Longitudinal variable names.
#' @param family family object to specify the link function for outcome.
#' @param withlong a logical parameter to decide if the weighted linear mixed model depends on longitudinal predictors.
#' @import dplyr glmmLasso glmnet doParallel lme4 missForest REEMtree MASS nlme
#' @return The GLMM LASSO model under given weights.
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
  glmmmodel <- suppressWarnings(glmer(glmmform,data = train.data,weights = weights,family = family))
  return(glmmmodel)
}
