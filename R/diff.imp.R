#' Calculating the difference between imputed and original values.
#'
#' @param ori Original numerical vector.
#' @param imp Imputed numerical vector.
#'
#' @return Sum of squared deviance between original and sifted vector over the sum of squared original vector.
#'
#' @export
diff.imp<-function(ori,imp){
  return(sum((imp-ori)^2)/sum(ori^2))
}
