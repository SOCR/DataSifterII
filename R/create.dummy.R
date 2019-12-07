#' Construct dummy variables for factor time-invariant vectors.
#'
#' @param data A data frame in long format with time-varying and time-invariant variables.
#' @param cs.names A vecotor of names for the time-invariant variables.
#'
#' @return A long format data with time-invariant factor/categorical variables changed to dummy variables.
#'
#' @export
create.dummy <- function(data,cs.names=NULL){
  names.list <- names(data)
  if(is.null(cs.names)) cs.names <- names(data)
  datacs <- data %>% dplyr::select(cs.names)
  datalg <- data %>% dplyr::select(-cs.names)
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
