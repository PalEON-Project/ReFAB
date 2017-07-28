library(nimble)
library(splines)
library(maps)
library(plyr)
library(oce)
library(RCurl)

setwd("/Users/paleolab/babySTEPPS/")

load('stick.break.like.start.Rdata')

pred_code <- nimbleCode({
  
  sigma ~ dunif(0,5) #GELMAN PAPER
  
  # b[1,1] ~ dunif(0,145)
  
  # for(t in 2:T){
  # b[1,t] ~ T(dlnorm(log(b[1,t-1]),1/sigma^2),0,145)
  # }
  
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
  
  for(t in 1:T){
    p.true[t,1] ~ dbeta(exp.phi[t,1],exp.phi1[t,1])
    p.rel[t,1] <- p.true[t,1]
    
    for(i in 2:(I-1)){
      p.rel[t,i]  ~ dbeta(exp.phi[t,i],exp.phi1[t,i]) 
      p.true[t,i] <-  p.rel[t,i] * (1 - sum(p.true[t,1:(i-1)]))
    }	
    p.true[t,21] <- 1 - sum(p.true[t,1:20])
  }    
  
  for(j in 1:J){
    Y[j,] ~ dmulti(p.true[age.index[j,1],],n[j])
  }
  
})

site_number = unique(x.meta[x.meta$site.name=='Cub Lake',1])
ten.count.use = ten.count[which(x.meta$site.id==site_number),]

Y = as.matrix(ten.count.use)

sample.ages <- x.meta[x.meta[,1]==site_number,]$age_bacon
age.bins <- seq(0,10000,100)
age.index <- as.matrix(as.numeric(cut(sample.ages,breaks = age.bins,labels=seq(1:(length(age.bins)-1)))))

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

data.pred = list(Y = Y, b = rep(100,T))

constants.pred = list(beta = beta1.est.real, beta1 = beta2.est.real, I = I, J = J,
                      T = T, n = n, u = u, N0 = rep(0, (length(u)-1)), 
                      N1 = rep(0, (length(u))), N2 = rep(0, (length(u)+1)),
                      N3 = rep(0, (length(u)+2)), age.index = age.index)

inits.pred = list(sigma = 4.5) #

dimensions.pred = list(exp.phi = c(T,I), exp.phi1 = c(T,I), phi.first = c(T,I),
                       phi.first1 = c(T,I), Zb = dim(Zb), Y = dim(Y),
                       p.rel = c(T,I), p.true = c(T,I)) #  b = dim(inits.pred$b),

set.seed(0)

model_pred <- nimbleModel(pred_code, inits = inits.pred, constants = constants.pred,
                          data = data.pred, dimensions = dimensions.pred)

cm <- compileNimble(model_pred)

##Your challenge is to rewrite this so it handles the different pollen samples at different times cleanly. You should be able to do all the calcs for all times in a single site at once.

calclik <- nimbleFunction(
  setup = function(model, calcnode, parentnodes) {
    deps = model$getDependencies(as.vector(parentnodes))
  },
  run = function(biomasses=double(1), nSims=double()) {
    n <- length(biomasses)
    returnType(double(2))
    out = matrix(NA,n,100)
    for(t in 1:100){
      for(i in 1:n) {
        #set.seed(0)
        model$b[t] <<- biomasses[i]
        model$calculate() # so deterministic dependencies of b are updated
        
        out[i,t] = 0
        for(j in 1:nSims) {
          model$simulate(parentnodes[t]) # simulate p.rel and fill in resulting p.true
          model$calculate(deps) # fill in deterministic dependencies of parentnodes and calculate densities for dependencies (which should include calcnode)
          out[i,t] <- out[i,t] + exp(model$getLogProb(calcnode[t])) # add density values across iterations
        }
        out[i,t] = out[i,t]/nSims
      }
      
    }
    return(out)
  })

parentnodes <- cbind(paste0('shape1[',age.index,',1:20]'),paste0('shape2[',age.index,',1:20]'))
calcnode <- paste0('Y[',1:length(age.index),',1:20]')
rcalclik <- calclik(model_pred, calcnode = calcnode[1:7], parentnodes = parentnodes[1:7,])

ccalclik <- compileNimble(rcalclik, project = model_pred, showCompilerOutput = TRUE)
bvals <- seq(1,150, length = 50)
tmp = ccalclik$run(bvals, 50000)
plot(bvals, log(tmp))

# alternative that tries to reset seed and do calculation one biomass at a time, but for some reason this doesn't smooth things out...
for( i in seq_along(bvals)) {
  set.seed(0)
  tmp[i] = ccalclik$run(bvals[i], 50000)
}


calclik <- nimbleFunction(
  setup = function(model, calcnode, parentnodes, timebin) {
    deps = model$getDependencies(parentnodes)
  },
  run = function(biomasses=double(1), nSims=double()) {
    n <- length(biomasses)
    returnType(double(1))
    out = numeric(n)
    for(i in 1:n) {
      model$b[timebin] <<- biomasses[i]
      model$calculate()# so deterministic dependencies of b are updated
      
      out[i] = 0
      for(j in 1:nSims) {
        model$simulate(parentnodes) # simulate p.rel and fill in resulting p.true
        model$calculate(deps) # fill in deterministic dependencies of parentnodes and calculate densities for dependencies (which should include calcnode)
        out[i] <- out[i] + exp(model$getLogProb(calcnode)) # add density values across iterations
      }
      out[i] = out[i]/nSims
    }
    return(out)
  })

bvals <- seq(1,150, length = 100)
tmp.mat5 <-matrix(NA,length(age.index),length(bvals))

for(i in 46:52){
  rcalclik <- calclik(model = model_pred,
                      calcnode = paste0('Y[',i,',1:21]'),
                      parentnodes = c(paste0('shape1[',age.index[i],',1:20]'),
                                      paste0('shape2[',age.index[i],',1:20]')),
                      timebin = age.index[i])
  ccalclik <- compileNimble(rcalclik, project = model_pred)
  
  tmp.mat5[i,] = ccalclik$run(biomasses = bvals, nSims = 100000)
  
  a = log(tmp.mat5[i,])
  a[is.na(a)] = 0
  plot(bvals, exp(a - max(a))/-sum(a),typ='l')
  
  plot(bvals, a, main=age.index[i],pch=19)
}

save(tmp.mat4,file='tmp.mat4.Rdata')

pdf('cub.likes.pdf')
par(mfrow=c(3,3))
for(i in 1:25){
  plot(bvals, log(tmp.mat4[i,]),main=age.index[i],typ='l')
}
dev.off()

a = log(tmp.mat5[i,])
plot(bvals, exp(a - max(a))/-sum(a),typ='l')

plot(tmp.mat5[i,]-max(tmp.mat5[i,]))
## 89, 96, 99 fit non temporal model. no smoothing over time.
## Cub Lake data and code


