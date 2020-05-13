##########################################################################

" 
  This file is part of the Preprint 
  https://www.medrxiv.org/content/10.1101/2020.04.07.20053421v2.
  
  Written by Amir S Siraj
  February 2020
  #================================================================
  numerical functions

"
## linear interpolate
interpolateCols <- function(vec, nSubsteps)
{
  matInterp <- rep(0,length(vec)+(length(vec)-1)*nSubsteps)
  
  for(i in 1:(length(vec)-1)) {
    interpVals <- approx(vec[i:(i+1)], n=(nSubsteps+2))
    lilloc<-i+nSubsteps*(i-1)
    matInterp[lilloc:(lilloc+nSubsteps+1)] <- interpVals$y
  }
  return(matInterp)
}


#### n time-unit moving average

ma <- function(x, n=2,parcial=TRUE){
  res = x #set the first values
  
  if (parcial==TRUE){
    for(i in 1:length(x)){
      t<-max(i-n+1,1)
      res[i] = mean(x[t:i])
    }
    res
    
  }else{
    for(i in 1:length(x)){
      t<-max(i-n+1,1)
      res[i] = mean(x[t:i])
    }
    res[-c(seq(1,n-1,1))] #remove the n-1 first,i.e., res[c(-3,-4,...)]
  }
}


# Weibull distribution with given parameters

weibull_f<- function (shape, scale) {
  
  
  #### the shape parameter (k)
  k = shape
  
  ### scale parameter (lambda)
  lambda= scale
  
  tmax=40
  tby=0.1
  x = seq(0,tmax,tby)
  
  ##### Density function math
  fx = k/ lambda * (x/lambda)^(k-1) * exp(-(x/lambda)^k)
  
  plot(x, fx, type='l')
  
  ##### Density function using R-function
  y<- dweibull(x, shape=k, scale= lambda)
  
  layout(matrix(1:2, nrow=1))
  par(mar=c(2,2,2,2))
  plot(x,fx, type='l', xlab="days", ylab="density", main="Math function")
  plot(x,y, type='l', xlab="days", ylab="density", main="R function")
  
  return (list(mean=lambda* gamma(1+1/k),  sd= sqrt( 
    lambda^2 *  ( gamma(1+2/k) - (gamma(1+1/k))^2 )
  ))) 
}



#Gamma distribution with given parameters

gamma_f<- function (shape, scale) {

  #### the shape parameter (k)
  alpha = shape
  
  ### scale parameter (lambda)
  lambda= scale
  
  tmax=40
  tby=0.1
  x = seq(0,tmax,tby)
  
  ##### Density function math
  fx =  (1/lambda^alpha) * x ^(alpha-1)*exp(-1/lambda*x) / gamma(alpha)
  
  plot(x, fx, type='l')
  
  ##### Density function using R-function
  y<- dgamma(x, shape=alpha, scale=lambda)
  
  layout(matrix(1:2, nrow=1))
  par(mar=c(2,2,2,2))
  plot(x,fx, type='l', xlab="days", ylab="density", main="Math function")
  plot(x,y, type='l', xlab="days", ylab="density", main="R function")
  
  return (list(mean=alpha * lambda,  sd= sqrt(alpha* lambda^2)
    )) 
}
