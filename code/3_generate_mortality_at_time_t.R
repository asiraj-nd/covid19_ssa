##########################################################################
" 
  This file is part of the Preprint https://doi.org/10.1101/2020.04.07.20053421.
  
  Written by Amir S Siraj
  February 2020
##########################################################################
 This script generates estimates of age specific mortality based on  
  China's CFR and Ethiopia's age pyramid

##########################################################################
- Specific tasks

  - 
    
  - 
  - 
    
##########################################################################
"
####################
setwd("~/covid19_ssa/code")

#China CDC report
#http://weekly.chinacdc.cn/en/article/id/e53946e2-c6c4-41e9-9a9b-fea8db1a8f51

eth_dd<- read.csv("../data/afr/eth_age_specific_pop_mortality.csv", stringsAsFactors = FALSE)[,-1]
scenall<- read.csv("../data/afr/afr_scens_short.csv", stringsAsFactors = FALSE)[,-1]

#### use the following codes alternatively to choose betwen 1 or 2
incub=2

if (incub==1) point_inf<- read.csv('../output/afr/cum_inf_AFR_t_length_and_95CI_weibull.csv')
if (incub==2) point_inf<- read.csv('../output/afr/cum_inf_AFR_t_length_and_95CI_gamma.csv')

head(point_inf)

#china
agemort<- eth_dd[,2]

dim(point_inf)
####### convert infection to age specific infection
###### assuming equal hazard at all age group

###### Median infections then convert to symtomatic
ageinf<- t(apply(cbind(point_inf[,10],t( eth_dd[,c(1,1)][, rep(1,nrow(scenall))])),1, function(tr){tr[1] * tr[-1]}))
sympinf<- ageinf * .82

mortall<- apply(t(apply(sympinf,1, function(tr){tr*agemort/100})),1,sum)

#####################  95%CI Lower bound 
ageinf<- t(apply(cbind(point_inf[,11],t( eth_dd[,c(1,1)][, rep(1,nrow(scenall))])),1, function(tr){tr[1] * tr[-1]}))
sympinf<- ageinf * .82

mortall<- cbind(mortall, apply(t(apply(sympinf,1, function(tr){tr*agemort/100})),1,sum))

##################### 95% CI Upper bound
ageinf<- t(apply(cbind(point_inf[,12],t( eth_dd[,c(1,1)][, rep(1,nrow(scenall))])),1, function(tr){tr[1] * tr[-1]}))
sympinf<- ageinf * .82
mortall<- cbind(mortall, apply(t(apply(sympinf,1, function(tr){tr*agemort/100})),1,sum))

mortall<- data.frame(mortall)

names(mortall)<- c("median_mort","95CI_lower","95CI_upper")

if (incub==1) write.csv(round(mortall,0), "../output/afr/mort_at_time_t_all_scenarios_weibull.csv", row.names=FALSE, quote=F)
if (incub==2) write.csv(round(mortall,0), "../output/afr/mort_at_time_t_all_scenarios_gamma.csv", row.names=FALSE, quote=F)

