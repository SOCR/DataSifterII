
#1 sifter main function
###########################################################################
dataSifter<-function(level="indep",data,unstructured.names=NULL,subjID=NULL,batchsubj=1000,missingback=FALSE,usecore=2,col_preserve = 0.5, col_pct = 0.7,
                     nomissing=FALSE,k0=NA,k1=NA,k2=NA,k3=NA,k4=NA,maxiter=1){
  if(!level %in% c(NULL,"none","small","medium","large","indep")) stop("Invalid Level Specification")
  if(is.null(level) & is.na(k0+k1+k2+k3+k4)) stop("Invalid k Combination")

  if(!is.null(subjID))data<-data[,-which(names(data) %in%subjID)]
  originaldatacolnames<-colnames(data)

  #Thinning step
  data<-thinning(data,col_preserve, col_pct)[[1]]
  missinglist<-apply(data,1,function(x) which(is.na(x),arr.ind = TRUE))

  #Obfuscate date
  #unind<-which(names(data) %in%  unstructured.names)
  sddata <- sdDates(data)
  if(!is.null(sddata[[2]])){
    dateind <- sddata[[2]]
    data <- sddata[[1]]
    #data[,-unind] <- sddata[[1]]
    print(paste("Detected date column numbers:",dateind))
    Date <- data[,dateind]
    Datename <- names(data)[dateind]
    data <- data[,-dateind]
    #print(names(data))
    if(length(dateind)>1){
      Date <- as.data.frame(lapply(Date,function(x) Sifter_Dates(x)))
      names(Date) = Datename}
    else {Date<-Sifter_Dates(Date)}
  }
  #Separate into batches and do imputation if there are less than 1000 we don't separate by defult.
  if(nomissing==FALSE){
    unind<-which(names(data) %in%  unstructured.names)
    if(nrow(data)>batchsubj){
      data.sep<-sep(data,batchsubj)

      #First imputation
      data.imp<-firstImp(data.sep,unind,firstimputation = TRUE,maxiter=maxiter,cores=usecore)

      if(length(unind)!=0){data.imp<-cbind.data.frame(data[,unind],data.imp);names(data.imp)[1]=unstructured.names}
    } else{
      data.imp<-firstImp(data,unind,firstimputation = TRUE,maxiter=maxiter,cores = 1)
      if(length(unind)!=0){data.imp<-cbind.data.frame(data[,unind],data.imp);names(data.imp)[1]=unstructured.names}
    }
    #Make the original data into imputed compelete data.
    data<-data.imp
  }

  #Calculate distance
  unind<-which(names(data) %in%  unstructured.names)
  subjDist_sorted<-subjdist(data,unind)

  if(level=="none"&is.na(sum(k0+k1+k2+k3+k4))){
    print("No Obfuscation")
    return(data)
  }else if (level=="indep"&is.na(sum(k0+k1+k2+k3+k4))){
    rownumb<-nrow(data)
    colnumb<-ncol(data)
    #Obtain factor, Character, numeric lists
    listfac<-which(sapply(data, is.factor) == TRUE)
    listnum<-which(sapply(data, is.numeric) == TRUE)
    listchar<-which(sapply(data, is.character) == TRUE)
    #Build factor/char freq table
    if(length(listchar)>1){
      chartable<-plyr::llply(.data=data[,listchar],.fun=table)
    } else if(length(listchar)==1){chartable<-list(table(data[,listchar])); names(chartable)=names(data)[listchar]}
    if(length(listfac)>1){
      factable<-plyr::llply(.data=data[,listfac],.fun=table)
    } else if(length(listfac)==1){factable<-list(table(data[,listfac])); names(factable)=names(data)[listfac]}

    #Estimate num densities
    vardens<-lapply(data[,listnum],density)#create empirical distribution list
    siftedData<-data.frame(plyr::llply(.data=vardens,.fun=function(k) {sample(x=k$x,size =rownumb,prob=k$y/sum(k$y),replace=T)}))

    if(length(listfac)>0){
      siftedData.fac<-data.frame(plyr::llply(.data=factable,.fun=function(k) {sample(x=names(k),size =rownumb,prob=as.numeric(k/sum(k)),replace=T)}))
      siftedData<-cbind.data.frame(siftedData,siftedData.fac)
    }
    if(length(listchar)>0){
      siftedData.char<-data.frame(plyr::llply(.data=chartable,.fun=function(k) {sample(x=names(k),size =rownumb,prob=as.numeric(k/sum(k)),replace=T)}))
      siftedData<-cbind.data.frame(siftedData,siftedData.char)
    }
    siftedData<-siftedData[,originaldatacolnames]
    return(siftedData)
  } else{
    if(is.na(sum(k0+k1+k2+k3+k4))){
      if(level=="small"){k0=1;k1=0.05;k2=1;k3=0.1;k4=0.01}
      if(level=="medium"){k0=0;k1=0.25;k2=2;k3=0.6;k4=0.05}
      if(level=="large"){k0=0;k1=0.4;k2=5;k3=0.8;k4=0.2}
    }
    # #k0
    #Cutoff for distance
    distmin<-min(subjDist_sorted$Distance)
    dist_cut<-sd(subjDist_sorted$Distance)
    list_string<-which(names(data)%in%unstructured.names)
    cases_to_swap <- round(k4 * nrow(data))
    if (is.null(list_string)) {
      print("No unstructured features to swap")
      data_k0 <- data #no unstructrued-sifted data=original one
    } else {
      if (k0 == 0) {
        data_k0 <- data #no swap-sifted data=original one
        print("No unstructured features to swap")
      } else {
        data_k0 <- data
        if(length(list_string)!=0){ B <- data
        for (i in 1:dim(data)[1]) {
          swap_subset_1 <- subjDist_sorted[which(subjDist_sorted[, 1] == i)[1:cases_to_swap], ]#subjDist_sorted=distance between cases
          swap_subset_2 <- subjDist_sorted[which(subjDist_sorted[, 2] == i)[1:cases_to_swap], ]
          swap_subset_1<-swap_subset_1[complete.cases(swap_subset_1),]
          swap_subset_2<-swap_subset_2[complete.cases(swap_subset_2),]
          swap_subset <- rbind(swap_subset_1, swap_subset_2)
          swap_subset <- swap_subset[order(swap_subset$Distance), ]
          swap_subset<-swap_subset[which(swap_subset$Distance<=(distmin+2*dist_cut)),]
          num_case<-ifelse(dim(swap_subset)[1]<cases_to_swap,dim(swap_subset)[1],cases_to_swap)
          if(num_case>0){
            case_to_swap <- sample(swap_subset[1:num_case, 1], 1)#there is something seriously wrong in this statement
            #print(case_to_swap)
            b1 <- B[i, list_string]
            b2 <- B[case_to_swap, list_string]
            B[i, list_string] <- b2
            B[case_to_swap, list_string] <- b1
            data_k0[c(i, case_to_swap), list_string] <- B[c(i, case_to_swap), list_string]
          }
        }
        }
      }
      print("Swap unstructred done")
    }

    #k1k2
    list_string<-which(names(data_k0)%in% unstructured.names)#no unstructured result in NULL value
    #work begin
    list5 <- c(which(sapply(data_k0, is.numeric) == TRUE), which(sapply(data_k0, is.factor) == TRUE),which(sapply(data_k0, is.character) == TRUE))# list 5 are numerical/factorial/characteristic features
    list5_cols <- as.numeric(list5)
    if(length(list_string)!=0){list5_cols<-list5_cols[-which (list5_cols %in% list_string)]}
    if (k1 == 0) {
      print("No imputation necessary, steps k1 and k2 skipped")
      data_k2 <- data_k0
    } else {
      C <- data_k0[, list5_cols]
      for (j in 1:k2) {
        C <- missForest::prodNA(C, noNA = k1)#introduce NA values equal to k1%
        SepC<-sep(C,numsample = batchsubj)
        C<-firstImp(sepdata = SepC,unstructured.list = list_string,firstimputation = FALSE,cores = usecore,maxiter=maxiter)
      }
      data_k2 <- data_k0
      data_k2[, list5_cols] <- C
    }
    print("Artifical missingness and imputation done")

    cases_to_swap <- round(k4 * nrow(subjDist_sorted))
    if (k3 == 0) {
      data.obfusc <- data_k2
      list_feature_obfusc <- list5_cols
    } else {
      list_feature_obfusc <- sample(list5_cols, round(k3 * length(list5_cols)))# select features to obfuscate
      data.obfusc <- data_k2
      D <- data_k2
      for (i in 1:nrow(data.obfusc)) {
        print(i)
        #################################Swap step is causing problems#################################
        swap_subset_1 <- subjDist_sorted[which(subjDist_sorted[, 1] == i)[1:cases_to_swap], ]
        swap_subset_2 <- subjDist_sorted[which(subjDist_sorted[, 2] == i)[1:cases_to_swap], ]#pick distances that has ith case in it
        swap_subset <- rbind(swap_subset_1, swap_subset_2)
        swap_subset <- swap_subset[order(swap_subset$Distance), ]
        swap_subset<-swap_subset[which(swap_subset$Distance<=(distmin+2*dist_cut)),]
        num_case<-ifelse(dim(swap_subset)[1]<cases_to_swap,dim(swap_subset)[1],cases_to_swap)
        if(num_case>0){
          case_to_swap <- sample(swap_subset[1:num_case, 1], 1)
          b1 <- D[i, list_feature_obfusc]
          b2 <- D[case_to_swap, list_feature_obfusc]
          D[i, list_feature_obfusc] <- b2
          D[case_to_swap, list_feature_obfusc] <- b1
          data.obfusc[c(i, case_to_swap), list_feature_obfusc] <- D[c(i, case_to_swap), list_feature_obfusc]
        }
      }
    }
    print("Swapping structured data done")

    #return to orginal order of features
    if(!is.null(sddata[[2]])){
      data.obfusc <- cbind(data.obfusc,Date)
      if(length(sddata[[2]])==1){
        names(data.obfusc)[ncol(data.obfusc)]=Datename
      }
      data.obfusc <- cbind(data.obfusc,Date)}
    data.obfusc<-data.obfusc[,originaldatacolnames]
    #make originally NA cells still NA
    if(missingback==TRUE&length(missinglist)!=0){
      misrow<-sapply(missinglist,function(x)length(x)!=0)
      misrow<-which(misrow=="TRUE")
      nadata<-data.obfusc[misrow,]#pick rows that has missing

      for(i in 1:length(misrow)){
        nadata[i,c(missinglist[[misrow[i]]])]=rep(NA,length(missinglist[[misrow[i]]]))
      }
      data.obfusc[misrow,]<-nadata
    }
    return(data.obfusc)
  }#after eta!=0
}



#2.Imputation function
#############################################################################################


firstImp<-function(sepdata,unstructured.list=NULL,firstimputation=FALSE,cores=2,maxiter=1){
  numgroup<-length(sepdata)
  cl1<-parallel::makeCluster(cores)
  registerDoParallel(cl1)
  clusterEvalQ(cl1, c(library(parallel), library(doParallel), library(missForest))  )
  clusterExport(cl1, list("sepdata","unstructured.list","maxiter"),envir = environment())
  if(firstimputation==TRUE){
    #cl<-parallel::makeCluster(cores,type = "SOCK")
    if(!is.null(unstructured.list)){
      K<-foreach::foreach(i=1:numgroup,.combine=rbind.data.frame,.packages = c("missForest")) %dopar%{missForest(sepdata[[i]][,-c(unstructured.list,which(names(sepdata[[i]])=="ID"))],maxiter=maxiter)$ximp}
    }else{
      K<-foreach::foreach(i=1:numgroup,.combine=rbind.data.frame,.packages = c("missForest")) %dopar%{missForest(sepdata[[i]][,-c(which(names(sepdata[[i]])=="ID"))],maxiter=maxiter)$ximp}
    }


  } else {
    for(i in 1:numgroup){
      list.factor<-which(sapply(sepdata[[i]], is.factor) == TRUE)
      if(length(list.factor)==1){
        sepdata[[i]][,list.factor]<-sapply(sepdata[[i]][,list.factor],function(x) factor(x,levels = unique(x)))

      } else if (length(list.factor)>1){
        sepdata[[i]][,list.factor]<-do.call(cbind.data.frame,lapply(sepdata[[i]][,list.factor],function(x) factor(x,levels = unique(x))))
      }
    }
    K<-foreach::foreach(i=1:numgroup,.combine=rbind.data.frame,.packages = c("missForest")) %dopar%{
      missForest(sepdata[[i]][,-which(names(sepdata[[i]])=="ID")],maxiter=maxiter)$ximp}
  }
  parallel::stopCluster(cl1)
  ID<-numeric()
  for(i in 1:length(sepdata)){
    ID<-c(ID,sepdata[[i]]$ID)
  }
  K$ID<-ID
  K<-K[order(K$ID),]
  K<-K[,-which(names(K)=="ID")]
  return(K)
}

#3.Separate a big dataset into multiple batches.
#############################################################################################

sep<-function(data,numsample=1000){#change the name to batching

  #generate case ID
  data$ID<-1:nrow(data)
  #Do ranodm sampling without replacement with the case ID
  group<-list()
  IDs<-data$ID #IDs<-1:2002
  numgroup<-floor(nrow(data)/numsample)

  for(i in 1:(numgroup)){
    group[[i]]=sample(IDs,size = numsample,replace = FALSE)
    IDs<-IDs[! IDs %in% group[[i]]]
  }

  group[[numgroup]]<-c(IDs,group[[numgroup]])
  #subsetting the original dataset according to case ID gorups and assign names or save into a list
  sepdata<-list()
  for(i in 1:numgroup){
    sepdata[[i]]<-data[group[[i]],]
  }
  return(sepdata)
}

#4. Subject distances
#############################################################################################
subjdist<-function(completedata,unstructured.list=NULL){
  library(wordspace)
  library(cluster)
  if(!is.null(unstructured.list)){completedata[,unstructured.list]<-sapply(completedata[,unstructured.list],as.character)}
  numcols<-which(sapply(completedata,class) %in% c("numeric","integer"))
  catcols<-which(sapply(completedata,class) %in% c("character","factor"))
  if(!is.null(unstructured.list))catcols<-catcols[-which(catcols%in% unstructured.list)]
  numdist<-as.dist(dist.matrix(as.matrix(completedata[,numcols]),method="euclidean"))#include the scale before calculating it
  numdist<-(numdist-min(numdist))/(max(numdist)-min(numdist))#scale to 0-1
  comb_subjects = combn(1:nrow(completedata), 2)
  if(length(catcols)<=1){#ignored when single catagory data is available
    b = rbind(comb_subjects, numdist)
    subjDist = data.frame(t(b))
    names(subjDist) = c("Case j", "Case i", "Distance")
    subjDist_sorted <- subjDist[order(subjDist$Distance), ]
  } else {
    catdist<-daisy(completedata[,catcols],metric="gower")
    catdist<-(catdist-min(catdist))/(max(catdist)-min(catdist))
    b = rbind(comb_subjects, numdist,catdist)
    subjDist = data.frame(t(b))
    names(subjDist) = c("Case j", "Case i", "Numdist","catdist")
    subjDist$Distance<-subjDist$Numdist*length(numcols)/ncol(completedata)+subjDist$catdist*length(catcols)/ncol(completedata)
    subjDist<-subjDist[,-c(3,4)]
    subjDist_sorted <- subjDist[order(subjDist$Distance), ]
  }
  return(subjDist_sorted)
}

#5. Remove non-informative variables
#############################################################################################
thinning<-function(data,col_preserve=.5,col_pct=.7){

  #Preprocessing Delete features that are 70% or more missing under col_preserve
  leastcol=floor(col_preserve*ncol(data))
  misscol<-which(apply(data,2,function(x) sum(is.na(x))/nrow(data))>=col_pct)
  ordered.misscolpct<-sort(apply(data,2,function(x) sum(is.na(x))/nrow(data)),decreasing = TRUE,na.last = TRUE)

  if(length(misscol)==0){
    data.new<-data
  } else if(length(misscol)<=leastcol ) {
    data.new<-data[,-misscol]
  } else if(length(misscol)>leastcol){
    misscol<-which(apply(data,2, function(x) sum(is.na(x))/nrow(data))>=ordered.misscolpct[leastcol])
    data.new<-data[,-misscol]
  }

  return(list(data.new,misscol))# output=deleted data[[1]], deleted missing cols[[2]], deleted missing rows[[3]]
}

#6. Detect dates and cast into "%m/%d/%Y" (returns a dataset and indexes of the date variables)
#############################################################################################

sdDates <- function(data) {
  patterns = c('[0-9][0-9][0-9][0-9]/[0-9][0-9]/[0-9][0-9]',
               '[0-9][0-9]/[0-9][0-9]/[0-9][0-9][0-9][0-9]',
               '[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]',
               'January|February|March|April|May|June|July|August|September|October|Novembber|December+[ ][0-9][0-9],[0-9][0-9][0-9][0-9]')
  formatdates = c('%Y/%m/%d','%d/%m/%Y','%Y-%m-%d','%d%b%y',"%B %d, %Y")
  standardformat='%m/%d/%Y'
  dataSampleElement <- lapply(data,function(x) x[!is.na(x)][1:2])
  indDate <- c()
  for(i in 1:length(patterns)){
    for(j in 1:ncol(data)){
      if(grepl(patterns[i], dataSampleElement[[j]][1])){
        if(grepl(patterns[i], dataSampleElement[[j]][2])){
          indDate<-c(j,indDate)
          aux=as.Date(dataSampleElement[[j]][1],format=formatdates[i],origin="1970-01-01")}
        else{ warning(paste("Variable #",j,"might have an inconsistent date format."))
          aux=NA}
        if(!is.na(aux)){
          data[,j]<-as.Date(data[,j],format=formatdates[i])
        }}
    }
  }
  return(list(data,indDate))
}

#7. Sifting a standardized date column
#############################################################################################

Sifter_Dates<-function(date_col){
  n=length(date_col)
  #date_col<-as.Date(date_col)
  date_col<-data.frame(id=1:n,date_col)
  date <- date_col[order(date_col[,2]),2]
  id <- date_col[order(date_col[,2]),1]
  #Calculate the distance among all the cases
  datedist <- sapply(2:length(date),function(x) date[x]-date[x-1])
  datedist <- c(10^99,datedist,as.numeric(Sys.Date()-date[length(date)]))
  error <- sapply(1:length(date), function(x)
    floor(rtruncnorm(1,a=-min(datedist[x],datedist[x+1]),b =min(datedist[x],datedist[x+1]),
                     mean = 0,sd = 0.5*min(datedist[x],datedist[x+1]))))
  # Follows a normal distribution with sd= 0.5*min(dist[i],dist[i+1])
  date_col <- data.frame(id=id,sifterdate=as.Date(date+error,origin="1970-01-01"))
  return(date_col[order(date_col$id),2])
}

#8. Unstructed DTM
#############################################################################################
un_DTM <- function(unstrvar,nterms=20,wb=NULL){
  corp<-tm::Corpus(tm::VectorSource(unstrvar))
  corp<-tm::tm_map(corp,stripWhitespace)
  corp<-tm::tm_map(corp,removePunctuation)
  corp<-tm::tm_map(corp,removeNumbers)
  corp <- tm::tm_map(corp, removeWords, stopwords("english"))

  dtm<-tm::DocumentTermMatrix(corp)
  dtm1<-tm::removeSparseTerms(dtm,0.99)
  dtm1<-data.frame(as.matrix(dtm1))
  freqdtm<-sort(colSums(dtm1),decreasing = TRUE)

  if(!is.null(wb)){wb <- wb[which(wb %in% names(dtm1))]}

  if(is.null(wb)){
    topfreq<-names(freqdtm)[1:nterms]
  }else{
    topfreq<-union(names(freqdtm)[1:nterms],wb)}
  dtm1 <- dtm1 %>% select(topfreq)
  dtm1$TermCombo <- apply(dtm1,1,function(x) paste(names(dtm1)[which(x>0)],collapse = ","))
  return(dtm1)
}


#9. Longitudinal Data Obfuscation
#############################################################################################



#10. Image features obfuscation
#############################################################################################
