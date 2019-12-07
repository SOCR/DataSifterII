#' Define link function for GLMM
#'
#' @param data A data frame in long format with outcome variable.
#' @param name the name of the outcome varaible.
#'
#' @return Proper link function.
#'
#' @export
checkLink<-function(data,name){
  outcome<-data[,name]
  if(length(levels(factor(outcome)))==2){
    link <- binomial(link = "logit")
  }else if(class(outcome)=="factor" & length(levels(outcome))>2 & length(levels(outcome))<20){
    link <- binomial(link = "logit")
  }else link <- gaussian(link = "identity")
  return(link)
}
