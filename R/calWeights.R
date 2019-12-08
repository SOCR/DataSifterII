#' Calculate sampling weights for raw data imputation.
#'
#' Use logistic regression or GLMM with LASSO regularity to model missingness for time-varying correlated data.
#' Output weights and full observation rows for each time-varying varaible.
#'
#' @param data Original data (numerical variables only) with missing values.
#' @param lcolnames A vector of longitudinal variables names.
#' @param w Type of sampling weights or missingness level. "peri" is to consider weights on subject level, which means any subjects with partial missing would be excluded from complete cases. "perij" is to consider weights on subject and time level.
#' Only subjects with all time points missing would be excluded from complete cases.
#' @param timevar The time variable or cluster varaible name.
#'
#' @param method Missingness model. If method = "param", the function utilize logistic regression ("peri") or GLMM ("perij") for missingness model.
#' If method = "nonparam", the function utilize random forest ("peri") for missingness model.
#'
#' @return A list contains two elements.
#' \itemize{
#'    \item Element 1 - Inversed probability of being observed.
#'    \item Element 2 - Row numbers of complete records in the original data.
#' }
#'#' @import dplyr glmmLasso glmnet doParallel lme4 missForest REEMtree MASS nlme
#'
#'@export
calWeights <- function(data,lcolnames,ID,w=c("peri","perij"),timevar,method="param"){
  #Input inspection
  if(class(lcolnames)!="character") stop("Invalid longitudinal column names.")
  if(!w %in% c("peri","perij")) stop("Weights is eighter per person(“peri”) or per person time(“perij”).")
  nlcol <- length(lcolnames)
  wlist <- list(nlcol)
  lInd <- which(names(data) %in% lcolnames)
  idInd <- which(names(data) == ID)
  #Build models for missingness per variable
  if(w=="peri"){
    #limit the data into one patient 1 line.
    w_rowidx <- lapply(lcolnames, function(x) {
      misid <- data[which(is.na(data[,x])),ID]
      #datasub <- data[!data[,ID] %in% unique(misid),]
      datasub <- data %>% group_by(ID) %>% filter(row_number()==1) %>%
        mutate(mis=!(ID %in% unique(misid))) %>% ungroup()%>%
        dplyr::select(-c(lcolnames,ID,timevar))
      if(method=="param"){
        modelmis <- cv.glmnet(x=data.matrix(datasub[,-which(names(datasub)=="mis")]),
                              y=as.numeric(datasub$mis), alpha=1,family = "binomial")
        modelmis <- glmnet(x=data.matrix(datasub[,-which(names(datasub)=="mis")]),
                           y=as.numeric(datasub$mis),lambda = modelmis$lambda.min)
        #print(modelmis$beta)
        testdata<-as.matrix(data[!data[,ID] %in% unique(misid),]%>%
                              dplyr::select(names(data)[names(data) %in% names(datasub)]))
        weights <- predict(modelmis,newx=data.matrix(testdata) ,type = "response")
        rownumb <- as.numeric(rownames(weights))
        } else {
        mtryint <- floor(sqrt(ncol(datasub)))
        mtryseq <- seq(2:min(mtryint+20,ncol(datasub)))
        selectmtry <- sapply(mtryseq, function(x) randomForest(factor(mis)~.,data = datasub,ntree=1000,mtry=x)$err.rate)
        mtryopt <- mtryseq[which.min(apply(selectmtry,2,function(x) min(x)))]
        modelmis <- randomForest(factor(mis)~.,data = datasub,ntree=1000,mtry=mtryopt)
        testdata<-as.matrix(data[!data[,ID] %in% unique(misid),]%>%
                              dplyr::select(names(data)[names(data) %in% names(datasub)]))
        weights <- predict(modelmis,newdata=testdata ,type = "prob")[,2]
        rownumb <- as.numeric(names(weights))
        }
      list(1/weights,rownumb)
    })
  } else {
    w_rowidx <- lapply(lcolnames, function(x) {
      misid <- data[which(is.na(data[,x])),ID]
      misfreqtb <- data.frame(ID=names(table(misid)),freq=as.numeric(table(misid)))
      datasub <- data
      mis <- ifelse(is.na(datasub[,x]),0,1)
      datasub <- cbind(datasub,mis) %>% dplyr::select(-lcolnames)
      datasub[,ID] <- as.factor(as.character(datasub[,ID]))
      #Here ID must equal to "ID"
        #This gives the first Y all 1's check why!!!
        modelmis <- filter.lambda(data=datasub,outcome = "mis",ID=ID,lambdaseq = seq(0,500,5),family = binomial(link = "logit"),startPQL = FALSE)[[1]]
        weights <- predict(modelmis,newdata = datasub[mis !=0,],type = "response")
      #print(weights)
      rownumb <- which(datasub$mis !=0)
      list(1/weights,rownumb)
    })
  }
  return(w_rowidx)
}
