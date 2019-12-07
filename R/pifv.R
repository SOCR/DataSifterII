#' Percent of Identical Feature Values
#'
#' Calculate the Percent of Identical Feature Values per record for the Sifted dataset.
#'
#' @param originaldata Original data frame.
#' @param sifteddata Sifted data frame.
#'
#' @return A vector of Percent of Identical Feature Values for each row in the Sifted data frame.
#'
#' @export
pifv <- function(originaldata,sifteddata){
  if(nrow(originaldata)!=nrow(sifteddata)) stop("Dimensions does not match.")

  pctfeatures_match<-function(originaldata,sifteddata,casenum1,casenum2){
    pctmatch<-length(which(originaldata[casenum1,]==sifteddata[casenum2,]))/(ncol(originaldata)-length(which(is.na(originaldata[casenum1,]))))
    return(pctmatch)
  }

  PIFV<-sapply(1:nrow(originaldata), function(x)
    pctfeatures_match(originaldata,sifteddata,x,x))
  return(PIFV)
}
