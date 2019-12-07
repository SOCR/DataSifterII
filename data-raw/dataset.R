#Simulation data generate
size=100
noisevar=5
set.seed(4)
x1 <- rnorm(size,mean=1,sd=5)
x2 <- rnorm(size,mean=4,sd=10)
x3 <- rnorm(size,mean=2,sd=4)
x4 <- rnorm(size,mean=6,sd=20)
x5 <- rnorm(size,mean=10,sd=30)


ni <- sample(c(1:10),size,replace = TRUE)
visit <- unlist(sapply(ni,function(x) 1:x))
e <- rnorm(length(visit),0,2)
e1 <- rnorm(length(visit),0,2)
ID <- unlist(sapply(1:size,function(x) rep(x,ni[x])))
visitdata<- data.frame(ID,visit,e,e1)

#Generate bi and eij, j=1,2,3,4
b <- rnorm(size,0,1)
b1 <- rnorm(size,0,1)
ID <- c(1:size)

#Combining
sim1<-data.frame(ID,x1,x2,x3,x4,x5,b,b1)
sim1 <- sim1 %>% right_join(visitdata,by="ID")
sim1 <- sim1 %>% mutate(Y=10+20*x1-15*x2-6*x3+0.8*visit+b+e)
sim1 <- sim1 %>% mutate(K=15+12*Y+5*x4+10*x5+2*visit+b1+e1)


noise <- data.frame(sapply(1:noisevar,function(x) rnorm(dim(sim1)[1],x*2.73,10)))
names(noise) <- paste("x",c(6:(5+noisevar)),sep = "")

sim1 <- cbind(sim1,noise)
sim1 <- sim1 %>% dplyr::select(-c(b1,e1,e,b))

misindY<-numeric()
#Everyone has a different bi
ni <- sim1 %>% group_by(ID) %>% mutate(ni=n())
ni <- ni %>%  filter(row_number()==1) %>% dplyr::select(ID,ni) %>% mutate(b1=rnorm(1,0,2))
b1 <- do.call("rbind",lapply(1:nrow(ni), function(x) data.frame(b1=rep(ni$b1[x],ni$ni[x]),ID=rep(ni$ID[x],ni$ni[x]))))
b1 <- b1[order(sim1$ID),]

set.seed(4)
for(i in 1:dim(sim1)[1]){misindY[i]<-rbinom(1,1,1/(1+exp(-2+20*sim1$x1[i]+sim1$visit[i]-b1$b1[i])))}
#table(misindY)

misindK<-numeric()

set.seed(4)
for(i in 1:dim(sim1)[1]){misindK[i]<-rbinom(1,1,1/(1+exp(-1-0.5*sim1$x5[i]+5*sim1$visit[i]-b1$b1[i])))}
#table(misindK)

sim1$Y[(misindY==1)]<-rep(NA,length(which(misindY==1)))
sim1$K[(misindK==1)]<-rep(NA,length(which(misindK==1)))

write.csv(sim1, file = "~/Box Sync/Longitudinal Sifter/Rpackage/DataSifter2/data-raw/sim.csv",row.names = FALSE)
setwd("~/Box Sync/Longitudinal Sifter/Rpackage/DataSifter2/data-raw/")
sim <- read.csv("sim.csv",na.strings = c("",NA))
usethis::use_data(sim,overwrite = TRUE)
