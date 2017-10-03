library(nimble) #, lib.loc='/tmp/nim063'
library(splines)
library(maps)
library(plyr)
library(oce)
library(RCurl)

load("~/babystepps/Data/calibration.data.Rdata") #from get.data.R
load("~/babystepps/Data/prediction.data.Rdata") #from get.data.R
load(file = '~/babystepps/Data/nimble.betas_1_22016-12-02.Rdata')

source("~/babySTEPPS/Workflow Code/utils/bs_nimble.R")

Z.knots = bs(biomass,intercept=TRUE,df=5)

beta1.est.real = matrix(colMeans(samples.mixed[100:nrow(samples.mixed),1:105]),ncol(Z.knots),ncol(Y))
beta2.est.real = matrix(colMeans(samples.mixed[100:nrow(samples.mixed),106:210]),ncol(Z.knots),ncol(Y))

#plots a confidence interval around an x-y plot (e.g. a timeseries)
ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}

dbetabin <- nimbleFunction(
  run = function(x = double(0), alpha = double(0), beta = double(0), size = double(0), 
                 log = integer(0, default = 0)) {
    returnType(double(0))
    logProb <- lgamma(size+1) - lgamma(x+1) - lgamma(size - x + 1) +
      lgamma(alpha + beta) - lgamma(alpha) - lgamma(beta) +
      lgamma(x + alpha) + lgamma(size - x + beta) - lgamma(size + alpha + beta)
    if(log) return(logProb)
    else return(exp(logProb))
  })

rbetabin <- nimbleFunction(
  run = function(n = integer(0), alpha = double(0), beta = double(0), size = double(0)) {
    returnType(double(0))
    if(n != 1) print("rbetabin only allows n = 1; using n = 1.")
    p <- rbeta(1, alpha, beta)
    return(rbinom(1, size = size, prob = p))
  })

order1 <- TRUE # for first order model; set to FALSE for 2nd order

pred_code <- nimbleCode({
  
  sigma ~ dunif(0,50) #GELMAN PAPER #5

  if(order1) {
    b[1] ~ dunif(0,145)
    for(t in 2:T){
      b[t] ~ T(dnorm(b[t-1],1/sigma^2),0,145)
    }
  } else  {
    b[1] ~ dunif(0,145)
    b[2] ~ dunif(0,145)
    for(t in 3:T){
      b[t] ~ dnorm(2*b[t-1] - b[t-2],1/sigma^2)
    }
  }
  
  for(t in 1:T){
    Zb[t,1:5] <- bs_nimble(b[t], u[1:3], N0[1:2], N1[1:3], N2[1:4], N3[1:5])
  }
  
  for(i in 1:I){
    for(t in 1:T){
      phi.first[t,i] <- sum(Zb[t,1:5] %*% beta[1:5,i])
      phi.first1[t,i] <- sum(Zb[t,1:5] %*% beta1[1:5,i])
    }
  }
  
  for(t in 1:T){
    for(i in 1:I){
      exp.phi[t,i] <- exp(phi.first[t,i])
      exp.phi1[t,i] <- exp(phi.first1[t,i])
    }
  }
  
  for(j in 1:J){
    Y[j, 1] ~ dbetabin(exp.phi[age.index[j,1], 1], exp.phi1[age.index[j,1], 1], n[j])
    for(i in 2:(I-1)){
      Y[j, i] ~ dbetabin(exp.phi[age.index[j,1], i], exp.phi1[age.index[j,1], i], n[j] - sum(Y[j,1:(i-1)]))
    }
  }
  
})

#Here's where you pick which lake you want to run
site_number = unique(x.meta[x.meta$site.name=='Cub Lake',1])
ten.count.use = ten.count[which(x.meta$site.id==site_number),]

Y = as.matrix(ten.count.use)

sample.ages <- x.meta[x.meta[,1]==site_number,]$age_bacon
age.bins <- seq(0,10000,100)
age.index <- as.matrix(as.numeric(cut(sample.ages,breaks = age.bins,labels=seq(1:(length(age.bins)-1)))))

tmp <- data.frame(cbind(age.index, Y))
names(tmp)[1] <- 'age.index'

Y2 <- aggregate(tmp, by = list(tmp$age.index), FUN = sum)

Y <- as.matrix(Y2[ , -c(1,2)])
age.index = as.matrix(Y2[,1])

Z.knots = Z
T = length(age.bins)-1
I = ncol(Y)
K = ncol(Z.knots)
J = length(age.index)
n = rowSums(Y)
Zb = matrix(NA,T,K)
phi.first = matrix(NA,T,I); exp.phi = phi.first
#beta.est = matrix(colMeans(samples1[100:nrow(samples1),]),K,I) 
new.biomass = seq(1,200,1)
Z.new = matrix(0,nrow=length(new.biomass),ncol=K)

u<-c(rep(attr(Z.knots,"Boundary.knots")[1],1),attr(Z.knots,"knots"),rep(attr(Z.knots,"Boundary.knots")[2],1))

data.pred = list(Y = Y)

constants.pred = list(beta = beta1.est.real, beta1 = beta2.est.real, I = I, J = J,
                      T = T, n = n, u = u, N0 = rep(0, (length(u)-1)), 
                      N1 = rep(0, (length(u))), N2 = rep(0, (length(u)+1)),
                      N3 = rep(0, (length(u)+2)), age.index = age.index)

inits.pred = list(b = rep(10, T),sigma = 4.5)#logb = matrix(log(10),1,T) #b = matrix(10,1,T), 

dimensions.pred = list(exp.phi = c(T,I), exp.phi1 = c(T,I), phi.first = c(T,I),
                       phi.first1 = c(T,I), Zb = dim(Zb), Y = dim(Y))

set.seed(0)

source('~/babySTEPPS/Workflow Code/samplers/samplers.R')
model_pred <- nimbleModel(pred_code, inits = inits.pred, constants = constants.pred,
                          data = data.pred, dimensions = dimensions.pred)
spec.pred <- configureMCMC(model_pred, thin = 10, print = FALSE,
                           control = list(log=TRUE))#,control = list(log=TRUE)


smp <- spec.pred$getSamplers()
for(i in 1:length(smp)) {
  if(smp[[i]]$name == 'RW' && smp[[i]]$target != 'sigma') {
    spec.pred$removeSamplers(smp[[i]]$target)
    spec.pred$addSampler(smp[[i]]$target, type = 'RWt_trunc', control = list(log=TRUE, range = c(0,145)))
    spec.pred$addSampler(smp[[i]]$target, type = 'jointb', control = list(log = TRUE, range = c(0,145), weights = c(.7,.2)))  # this seems to help avoid getting stuck at low-lik values early in chain and leads to higher ESS, but sampling does take longer... 
  }
}

spec.pred$addMonitors(c("b")) 
Rmcmc.pred <- buildMCMC(spec.pred)
cm <- compileNimble(model_pred, Rmcmc.pred)

# don't initialize all b's at same value as that can lead to samples for sigma being driven to be very small, at least for a while
b1 <- rnorm(T, 25, 10)
b2 <- rnorm(T, 75, 10)
b3 <- rnorm(T, 125, 10)
b1[b1 < 0] <- 2
b3[b3 > 145] <- 144

samplesList <- runMCMC(mcmc = cm$Rmcmc.pred, niter = 100, nchains = 3,
                       inits = list(list(b = b1, sigma = 4.5),
                                    list(b = b2, sigma = 4.5),
                                    list(b = b3, sigma = 4.5)))


stop()

samplesList <- runMCMC(mcmc = cm$Rmcmc.pred, niter = 10000, nchains = 1,
                       inits = list(b = b3, sigma = 4.5))

sl = samplesList[[3]]

par(mfrow=c(2,2))
ts.plot(sl[,88])
ts.plot(sl[,91])
ts.plot(sl[,94])
bvals <- sl[500:1000,1:100]
bs <- colMeans(bvals)
ts.plot(bs)
lines(1:100, apply(bvals, 2, quantile, .025), col = 'red')
lines(1:100, apply(bvals, 2, quantile, .975), col = 'red')


# basic likelihood plot (for biomass at time 'tt'
if(F) {
  tt = 86
  dps = cm$model_pred$getDependencies(paste0('b[', tt, ']'))
  ind = which(age.index == tt)
  bs = seq(10, 145, by = 1)
  i=1; out = rep(0, length(bs))
  for( b in bs ) {
    cm$model_pred$b[tt] = b
    cm$model_pred$calculate(dps)
    out[i] = cm$model_pred$calculate(paste0('Y[', ind, ', 1:20]'))
    i = i + 1
  }
  
  plot(bs, exp(out), type = 'l')
}
