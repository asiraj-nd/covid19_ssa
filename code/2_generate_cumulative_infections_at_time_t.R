##########################################################################

" 
  This file is part of the Preprint 
  https://www.medrxiv.org/content/10.1101/2020.04.07.20053421v2.
  
  Written by Amir S Siraj
  February 2020
##########################################################################
 This script generates the median and CI  of cumulative infections 
  at a given time after epidemic trigger point for each scenario

##########################################################################
- Specific tasks

  - 
    
  - 
  - 
    
##########################################################################
"
#rm(list = ls())
####################
setwd("~/covid19_ssa/code")

####### load numerical functions 
source("0_numfunctions.R")

### read scenario setup file
scenall<- read.csv("../data/afr/afr_scens_short.csv", stringsAsFactors = FALSE)


#### the time t to generate cumulative estimates for
point_t = 366
iterlen=100

#### use the following codes alternatively to choose betwen
#### 1000 parameter values of Weibull or Gamma distrbuted 
incub=2

### read in simulation daily infections
if (incub==1) allsims.pack<-read.csv("../output/afr/all_simulations_AFR_weibull.csv", stringsAsFactors = FALSE)[,-(1:2)]
if (incub==2) allsims.pack<-read.csv("../output/afr/all_simulations_AFR_gamma.csv", stringsAsFactors = FALSE)[,-(1:2)]
allsims.pack<- allsims.pack[,1:point_t]

ndays<-ncol(allsims.pack)
x = 0:ndays

allpts<- NULL
for (scen in 1:nrow(scenall)) {
  rower<- (1:iterlen)+(scen-1)*iterlen
  realistic_r <- which(allsims.pack[rower,1] <  (scenall[scen,4]*2)) # remove doubling estimates - rare events
  dcases<- allsims.pack[rower,] [realistic_r,]
  cumcases<- t(apply(dcases,1, cumsum))
  dim(cumcases)
  thissim<- NULL
  #### generate median and 95%CI
  for (j in 1:ndays) {
    thisday<- cumcases[,j]
    thisquant<- c(quantile(thisday,p=c(0.025,.5,0.975), na.rm=TRUE))
    thissim<- cbind(thissim,c(thisquant))
  }
  
  allpts<-rbind(allpts,thissim[c(2,1,3),ndays])
  
}


dim(allpts)

if (incub==1) write.csv(cbind(scenall,allpts), "../output/afr/cum_inf_AFR_t_length_and_95CI_weibull.csv", row.names=F, quote=F)
if (incub==2) write.csv(cbind(scenall,allpts), "../output/afr/cum_inf_AFR_t_length_and_95CI_gamma.csv", row.names=F, quote=F)

