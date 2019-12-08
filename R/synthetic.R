#' Create fully synthetic longitudinal dataset.
#'
#' When data privacy is the major issue of time-varing data sharing, we consider this fully synthetic version of DataSifter II.
#' It obfuscate longitudinal data by generating random numbers within the neighborhood of a longitudinal variable at a given visit time.
#' Thus, the within-subject covariance structure might be changed by this process.
#'
#' @param data A data frame in long format.
#' @param lnames A vector of longitudinal variables names.
#' @param timevar The time variable or cluster varaible name.
#' @param ID Name of the ID variable in the dataset.
#' @param N The number of closest records to be considered as neighborhood.
#'
#' @return Fully synthetic Sifted data.
#' @import dplyr glmmLasso glmnet doParallel lme4 missForest REEMtree MASS nlme
#' @export
synthetic <- function(data,lnames,timevar,ID,N){#Data must be in long format

  myspread <- function(df, key, value) {
    # quote key
    keyq <- rlang::enquo(key)
    # break value vector into quotes
    valueq <- rlang::enquo(value)
    s <- rlang::quos(!!valueq)
    df %>% gather(variable, value, !!!s) %>%
      unite(temp, !!keyq, variable) %>%
      spread(temp, value)
  }

  #Variable range calculation.
  neighbor.range.sift <- function(sortedlist, varname,N,sdvec){
    varlist <- sortedlist[[which(names(sortedlist)==varname)]]
    nonnarows <- which(!is.na(varlist[,2]))
    sdvar <- sdvec[varname]
    sift.var <- sapply(nonnarows, function(rownum){
      cur.rank <- varlist[rownum,1]
      cur.num <- varlist[rownum,2]
      range_neighbor <- range(varlist[which(varlist[,1] %in% c((cur.rank-N):(cur.rank+N))),2],na.rm = TRUE)
      if(sum(is.na(range_neighbor))==1 & sdvar<(range_neighbor[2]-range_neighbor[1])){
        sifted_num <- runif(1,cur.num-0.5*sdvar,cur.num+0.5*sdvar)
      }else{
        sifted_num <- runif(1,range_neighbor[1],range_neighbor[2])
      }
      sifted_num
    })
    sifted_var <- rep(NA,dim(varlist)[1])
    sifted_var[nonnarows] <-sift.var
    return(sifted_var)
  }

  data.wide <- data %>% myspread(timevar,lnames) #ok
  wide_longivars <- grep(paste(paste("_",lnames,sep = ""),collapse = "|"),names(data.wide),value = TRUE)
  sub.data <- data.wide %>% dplyr::select(wide_longivars)
  sortedlist <- lapply(sub.data,function(x) cbind(rank(x,na.last = FALSE,ties.method = "first"),x))
  sdvec <- sapply(sortedlist,function(x) sd(x[,2],na.rm = TRUE))

  range <- ifelse(N%%2==0,N/2,(N-1)/2)
  sift.data.syn <-sapply(wide_longivars, function(x)
    neighbor.range.sift(sortedlist=sortedlist,varname = x,N=range,sdvec=sdvec[x]))

  data.wide[,wide_longivars]<-sift.data.syn

  data.long <- data.wide %>%
    gather(v, value, wide_longivars) %>%
    separate(v, c(timevar, "col")) %>%
    arrange(ID) %>%
    spread(col, value)
  data.long <- data.long %>% dplyr::select(names(data))

  return(list(data.long,data))
}
