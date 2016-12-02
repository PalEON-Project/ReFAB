
#low biomass might cause bias in 2nd derivative
#plot the squared second derivative across sites
#we don't have covariates so we don't want to think about ind. sites too much
#if a site is linearly increasing then it will appear stable
#map the first or second derivatives over time? sum(abs(first derivative [1:T]))

#ask andria about climate proxies

biomass <- function(x){
  #4*x #linear
  #rnorm(length(x),4*x,1) #linear with random effects
  #4*x^2 #parabolic
  #100/(1+exp(-.1*(x-50))) #logistic
  abs(rnorm(length(x),1000/(1+exp(-.1*(x-50))),1000) )#logistic with random effects
  #floor(x/3.14 +.5) + sin(x^2) + 50 #wiggly with positive trend
  #rnorm(length(x),(floor(x/3.14 +.5) + sin(x^2) + 50),10)#wiggly with positive trend and random effects
}

second.deriv <- function(x,h){
  (biomass(x + h) - 2*biomass(x) + biomass(x-h))/(h^2)
}

x<-seq(1,100,1)
h <- 1 #changes scale of the sum of the second derivative
par(mfrow=c(3,1))
plot(x,biomass(x),main="Function",typ="l")
plot(diff(biomass(x))/diff(x),main="First Deriv",typ='l')
abline(h=0,col="red",lwd=2)
mtext(paste('sum = ',sum(diff(biomass(x))/diff(x))),side=3)
plot(second.deriv(x,h),main="Second Deriv",typ='l')
abline(h=0,col="red",lwd=2)
mtext(paste('sum = ',sum(second.deriv(x,h))),side=3)
