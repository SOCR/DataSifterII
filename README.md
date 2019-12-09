# DataSifter II: Statistical Obfuscation of Sensitive Time-varying Correlated Data

## Authors
Nina Zhou, Lu Wang, Simeone Marino, Yi Zhao, and Ivo Dinov

## Example

### Load Package

Install both DataSifter I and II.
```{r}
library(devtools)
install_github("SOCR/DataSifterII")
install_github("SOCR/DataSifter")
library(DataSifterII)
library(DataSifter.lite)
```

### Create Sifted dataset

Obfuscate time-varying and time-invariant data seperately.

```{r}
data(sim)

#DS II on Time-varying data
set.seed(1234)
siftsim <- repSifter(data=sim,mispct = 0.2,lnames=c("Y","K"),
                     timevar = "visit",ID="ID")
sifttv <- siftsim[[1]] %>% dplyr::select(c("ID","visit","Y","K"))

#DS I on Time-invariant data
sim_cs <- sim %>% filter(visit==1) %>% dplyr::select(-c("visit","Y","K"))
sifter_sim_cs <- dataSifter(level="medium",data = sim_cs,subjID = "ID",nomissing = TRUE)

#Merging two together
siftsim_final <- left_join(siftsim[[1]] %>% dplyr::select(c("ID","visit","Y","K")),
                           cbind(ID=as.factor(sim_cs$ID),sifter_sim_cs),by="ID")

```
### Test Privacy and Utility

The generating models for Y and K are

E(Y)=10+20\*X<sub>1</sub>-15\*X<sub>2</sub>-6\*X<sub>3</sub>+0.8\*visit,

and E(K)=15+12\*Y+5\*X<sub>4</sub>+10\*X<sub>5</sub>+2\*visit.

```{r}
#Data Privacy
boxplot(pifv(sim,siftsim_final))

#Data utility
library(nlme)
#Y=10+20*x1-15*x2-6*x3+0.8*visit+b+e
Ymodel <- lme(Y~x1+x2+x3+visit,random = ~1|ID,data=siftsim_final)
summary(Ymodel)
intervals(Ymodel)
#K=15+12*Y+5*x4+10*x5+2*visit+b1+e1
Kmodel <- lme(K~Y+x4+x5+visit,random = ~1|ID,data=siftsim_final)
summary(Kmodel)
intervals(Kmodel)
```
