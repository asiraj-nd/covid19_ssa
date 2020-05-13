##########################################################################

" 
  This file is part of the Preprint 
  https://www.medrxiv.org/content/10.1101/2020.04.07.20053421v2.
  
  Written by Amir S Siraj
  February 2020
##########################################################################
  This script sets up random variables for parameter based on distributions, 
  assigns fixed parameters, scenarios, calls functions to populates all 
  pomp c++ snippets, creates the pomp object, sets up the iterative process 
  for simulation and saves the simulation results.

##########################################################################
- Specific tasks

  - 
    
  - 
  - 
    
##########################################################################
"
#rm(list = ls())
options(digits=10)

####################
setwd("~/covid19_ssa/code")

##########################################################################
####### load the pomp snippet setup functions, 
source("0_snippet_setup.R")

####### load numerical functions 
source("0_numfunctions.R")


#### load packages
library(inline)
library(pomp)
library(plyr)

# call function to populate the pomp process function ============================
rproc.ex1 = declare.rproc()


# call function to populate the pomp measurment function ==========================
rmeas.ex1 = declare.rmeas()

set.seed(123)

# infectious period fixed @5 days per Davies 2020
gamma.vec = rep(1/5,100)  

#### use the following codes alternatively to choose betwen
#### 1000 parameter values of Weibull or Gamma distrbuted 
incub=1

# Weibull incubation period per Lauer et al, 2020
if (incub==1) {
  xi.vec<- 1/rweibull(1000, shape=2.45, scale= 6.26)  
  weibull_f(2.45, 6.26)  # display the pdf  
}

# Gamma incubation period per Davies et al, 2020
if (incub==2) {
  xi.vec<- 1/rgamma(1000, shape=4, scale= 1.375) 
  gamma_f(4, 1.375) # display the pdf
}

##### this segment is now obsolute - left here for coding convenience
phi.vec<- xi.vec[1:100]*2         
phi.vec[which(phi.vec < gamma.vec)]<- gamma.vec[which(phi.vec < gamma.vec)]

######### set up the loop to setup pomp object and simulate based on scenarios

# length of time-series output in days
given.time = 548
given.start.time = 0

###### read scenario parameters file
scenall<- read.csv("../data/afr/afr_scens_short.csv")

###### simulation collection dataframe
allsims.pack  <- NULL
#scen=1
for (scen in 1:nrow(scenball)) {

    # Setup Pomp object's initial condition for state variables ====
    E.init = ceiling(scenall[scen,4]/2-0.01)
    I.init = floor(scenall[scen,4]/2+0.01)
    R.init = 0
    N.init= scenall[scen,2]
    C.init = scenall[scen,4]
    
    # declare pomp data object  ====================================
    data.ex1 = declare.data(given.time)
    
    mode(data.ex1)

    # initialize states  ============================================
    init.ex1 = declare.init(N.init, E.init, I.init, R.init, C.init)


    # declare pomp process object ===================================
    rproc.ex1 = declare.rproc()

    # create the pomp object ========================================
    ptm = proc.time()
      pomp.obj = declare.pomp(data.ex1, given.start.time, rproc.ex1, rmeas.ex1, init.ex1)
    print(proc.time()-ptm)

    set.seed(123)
    
    ##### run 100 simulations    
    for (iter in 1:100){

      given.gamma = gamma.vec[iter]
      given.xi = xi.vec[iter]   # obsolute
      given.phi = phi.vec[iter] # obsolute
      given.xi.vec = xi.vec
      given.theta = scenall[scen,6]
      given.R0 = scenall[scen,3]
      given.omega = scenall[scen,7]  ## effectiveness of facemask
      given.tau <-  scenall[scen,8] ## coverage of facemask
      given.epsilon = scenall[scen,5] ## effectivness of SD

      if(given.phi<given.gamma) given.phi = given.gamma  # obsolute - left for convenience
      
      #### beta parameter (fixed part of the force of infection)  
      given.beta = given.gamma * given.R0 * ( 1- given.tau * given.omega) *  (1-given.epsilon)

      print(given.beta)

      thisYCs<- NULL
      given.nsim = 1 ##### one simulation per iteration 

      #c('beta', 'phi', 'xi', 'tau', 'gamma')

      # simulate with this set of params ==================================
      # and those that radonmly selected in the process model =============
      sims1 = simulate.pomp(
      pomp.obj, given.beta, given.theta, given.phi, given.xi, given.tau, given.gamma, given.xi.vec, given.nsim)
      simdata = as.data.frame(sims1$sim)
      
      #### collect simulated new cases 
      allsims.pack<- rbind(allsims.pack,c(scen,iter,diff(c(C.init,simdata$YC))))
      print(c(scen,iter))
    }  ## iteration (iter of 100)
} # scenarion (scen of nrow(scenall))

dim(allsims.pack)

if (incub==1) write.csv(allsims.pack,"../output/afr/all_simulations_AFR_weibull.csv", row.names=F, quote=F)
if (incub==2) write.csv(allsims.pack,"../output/afr/all_simulations_AFR_gamma.csv", row.names=F, quote=F)
