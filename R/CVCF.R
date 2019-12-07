#' Last value carry forward and next value carry backward combination for initiating imputation.
#'
#' @param data A data frame in long format with missing value.
#' @param lnames A vector of the longitudinal varaible names.
#' @param timevar Name of the visit time variable.
#' @param ID Name of the ID variable.
#'
#' @export
CVCF<-function(data, lnames, timevar, ID){
  #data=dataname, lnames=longitudinal variable names, timevar= time variable name
  #Check if the inputs are valid?
  if(!"data.frame"%in% class(data)){stop("Invalid data input.")}
  if(class(lnames)!="character"){stop("Invalid longitudinal variable names.")}
  if(class(timevar)!="character"){stop("Invalid time variable names.")}
  #Calculate the mean of each time point
  ln <- length(lnames)
  time_means <- data %>% group_by(get(timevar)) %>%
    summarize_at(lnames,funs(mean(.,na.rm = TRUE)))
  names(time_means)[1]<-timevar
  #Obtain the missingness for each id and each timepoint each varaible.
  #When categorical is observed the rowSums would be non NA all the time.
  misind <- data %>% dplyr::select(ID,timevar,lnames) %>% mutate(rsum=rowSums(.[-c(1,2)])) %>%
    filter(is.na(rsum)) %>% dplyr::select(ID,lnames,timevar)
  misind$misvar <- apply(misind,1,function(x) paste(names(misind)[which(is.na(x))],collapse = ","))
  workdata <- data %>% dplyr::select(ID,lnames,timevar) %>%
    filter(ID %in% misind[,1]) %>% arrange(ID,get(timevar))
  #Check if the misvar is which one or multiple vectors
  #Slow version
  for(i in 1:nrow(misind)){
    misvarvec <- unlist(strsplit(misind$misvar[i],","))
    tempdata <- workdata %>% filter(ID %in% misind[i,1])
    tempdatab <- tempdata %>% filter(get(timevar)<misind[i,timevar])
    tempdataa <- tempdata %>% filter(get(timevar)>misind[i,timevar])
    for(j in misvarvec){
      if(dim(tempdatab %>% filter(!is.na(get(j))))[1]>0){
        tempdatab_sub <- tempdatab %>% filter(!is.na(get(j)))
        misind[i,which(names(misind)%in% j)] <- tempdatab_sub[nrow(tempdatab_sub),j]
      }else if(dim(tempdataa %>% filter(!is.na(get(j))))[1]>0){
        tempdataa_sub <- tempdataa %>% filter(!is.na(get(j)))
        misind[i,which(names(misind)%in% j)] <- tempdataa_sub[1,j]
      } else{
        misind[i,which(names(misind)%in% j)] <- as.numeric(time_means %>% filter(get(timevar)==misind[i,timevar]) %>% dplyr::select(j))
      }
    }
  }
  for(i in 1:nrow(misind)){
    var <- unlist(strsplit(as.character(misind[i,"misvar"]),","))
    varind <- which(names(data) %in% var)
    id <- as.character(misind[i,ID])
    t <- misind[i,timevar]
    data[which(as.vector(data[,ID]) == c(id) & as.vector(data[,timevar]) == c(t)),var] <- misind[i,var]
  }
  return(data)
}
