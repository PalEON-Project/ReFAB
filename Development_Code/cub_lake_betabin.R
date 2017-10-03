library(nimble)
library(splines)
library(maps)
library(plyr)
library(oce)
library(RCurl)

data.dir = "/Users/paleolab/babySTEPPS/Data/"
fig.dir = "/Users/paleolab/babySTEPPS/Figures/"
model.dir = "/Users/paleolab/babySTEPPS/Code/"

setwd("/Users/paleolab/babySTEPPS/")


load("add.bacon2.Rdata")
#load("nimble.betas2016-11-22.Rdata") #load("2016-05-31nimble.betas.Rdata") 
source(paste0(model.dir,"bs_nimble.R"))
load(file = 'nimble.betas_1_22016-12-02.Rdata')

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


pred_code <- nimbleCode({
  
  sigma ~ dunif(0,50) #GELMAN PAPER #5
  
  # Notes -- log biomass evolution
  #logb[1,1] ~ dunif(log(1),log(145)) #this makes it hard to jump modes because it's on a different scale
  # for(t in 2:T){
  #   logb[1,t] ~ T(dnorm(logb[1,t-1]),1/sigma^2),-Inf,log(145))
  # }
  
  #First Order -- biomass evolution
  b[1] ~ dunif(0,145)
  for(t in 2:T){
    b[t] ~ T(dnorm(b[t-1],1/sigma^2),0,145)
  }
  # #Second Order - What we decided we're going to use
  # b[1] ~ dunif(0,145)
  # b[2] ~ T(dnorm(b[1],1/sigma^2),0,145)
  # for(t in 3:T){
  #   b[t] ~ T(dnorm((2*b[t-1] - b[t-2]),1/sigma^2),0,145)
  # }
  
  #Second Order - Ann Version
  # b[1,1] ~ dunif(0,145)
  # b[1,2] ~ dunif(0,145)
  # for(t in 3:T){
  #   b[1,t] ~ T(dlnorm(log(2*b[1,t-1] - b[1,t-2]),1/sigma^2),0,145)
  # }
  
  #Second Order - Chris Version
  # logb[1,1] ~ dunif(log(2),log(145))
  # logb[1,2] ~ dunif(log(2),log(145))
  # for(t in 3:T){
  # logb[1,t] ~ T(dnorm(2*logb[1,t-1] - logb[1,t-2], 1/sigma^2), log(2), log(145))
  # }
  # 
  # for(t in 1:T){
  #   b[1,t] <- exp(logb[1,t])
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

  for(j in 1:J){
      Y[j, 1] ~ dbetabin(exp.phi[age.index[j,1], 1], exp.phi1[age.index[j,1], 1], n[j])
      for(i in 2:(I-1)){
        Y[j, i] ~ dbetabin(exp.phi[age.index[j,1], i], exp.phi1[age.index[j,1], i], n[j] - sum(Y[j,1:(i-1)]))
        
      }
  }
 
})

u<-c(rep(attr(Z,"Boundary.knots")[1],1),attr(Z,"knots"),rep(attr(Z,"Boundary.knots")[2],1))

x = new.pol1[new.pol1$age_bacon>=200,]
x = x[x$age_bacon<=10000,]

x.meta = x[,c('site.id','lat',"long","dataset.id","site.name","age_bacon")]

trees <- c("PINUSX","ALNUSX","JUGLANSX","ACERX","CUPRESSA","FRAXINUX","FAGUS","CYPERACE","LARIXPSEU","TSUGAX","QUERCUS","TILIA","BETULA","PICEAX","OSTRYCAR","ULMUS","ABIES","POPULUS")
other.trees <- c("TAXUS","NYSSA","CASTANEA","PLATANUS","SALIX","LIQUIDAM")
ten.count = matrix(0,nrow(x),length(trees)+3)
prairie <- c("ARTEMISIA","ASTERX","POACEAE","AMBROSIA","CHENOAMX","CORYLUS")
ten.count[,1] <- unlist(rowSums(x[,prairie]))
ten.count[,2] <- unlist(rowSums(x[,other.trees]))
ten.count[,3:(length(trees)+2)] <- as.matrix(x[,trees])
ten.count[,(length(trees)+3)] <- rowSums(x[,20:99]) - rowSums(ten.count)
colnames(ten.count)<-c("prairie","other trees",trees,"other herbs")

ten.count.save = ten.count
ten.count = round(ten.count.save)

counts <- Y[,rev(order(colMeans(Y)))]

ten.count <- ten.count[,colnames(counts)]

samples.pred.save = list()
biomassCI=list()

for(s in 1:length(unique(x.meta[,1]))){ #length(unique(x.meta[,1]))
  print(paste('working on number',s,'of 183',(s/183)*100,'% complete'))
  site_number = unique(x.meta[,1])[s]

  #site_number = unique(x.meta[x.meta$site.name=='Cub Lake',1])
ten.count.use = ten.count[which(x.meta$site.id==site_number),]

if(length(ten.count.use)>25*20&	min(x.meta[x.meta[,1]== site_number,]$age_bacon)<1000 & 
   max(x.meta[x.meta[,1]== site_number,]$age_bacon)>9000){
  
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

source('~/Downloads/newsampler/sampler_RWt.R')
source('~/Downloads/newsampler/sampler_jointb.R')
model_pred <- nimbleModel(pred_code, inits = inits.pred, constants = constants.pred,
                          data = data.pred, dimensions = dimensions.pred)
spec.pred <- configureMCMC(model_pred, thin = 10, print = FALSE,
                           control = list(log=TRUE))#,control = list(log=TRUE)


smp <- spec.pred$getSamplers()
for(i in 1:length(smp)) {
    if(smp[[i]]$name == 'RW' && smp[[i]]$target != 'sigma') {
        spec.pred$removeSamplers(smp[[i]]$target)
        spec.pred$addSampler(smp[[i]]$target, type = 'RWt', control = list(log=TRUE))
        spec.pred$addSampler(smp[[i]]$target, type = 'jointb', control = list(log = TRUE, weights = c(.7,.2)))  # this seems to help avoid getting stuck at low-lik values early in chain and leads to higher ESS, but sampling does take longer... 
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
samplesList <- runMCMC(mcmc = cm$Rmcmc.pred, niter = 5000, nchains = 3,
                       inits = list(list(b = b1, sigma = 4.5),
                                    list(b = b2, sigma = 4.5),
                                    list(b = b3, sigma = 4.5)))

samples.pred.save[[s]]<-samplesList

# rug(sample.ages/100)
dyn.unload(model_pred$nimbleProject$cppProjects[[1]]$getSOName())
biomassCI[[s]] <- apply(as.matrix(samplesList[[1]][,1:100],samplesList[[2]][,1:100],samplesList[[3]][,1:100]),2,quantile,c(0.025,0.5,0.975))
#samples.pred.good.inits.save[[s]] <- samples.pred.good.inits

}
save(biomassCI, file="biomass.CI13.Rdata")
save(samples.pred.save, file="samples.pred.save5.Rdata")
}

# basic likelihood plot (for biomass at time 'tt'
if(F) {
    tt = 96
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

pdf(paste0(Sys.Date(),unique(x.meta[x.meta[,1]==site_number,'site.name']),'.trace.pdf'))
bs = seq(10, 145, by = 1)
out.keep <- matrix(NA,100,length(bs))
par(mfrow=c(3,2))
for(t in age.index){ #
  
  tt = t
  dps = cm$model_pred$getDependencies(paste0('b[', tt, ']'))
  ind = which(age.index == tt)
  i=1; out = rep(0, length(bs))
  for( b in bs ) {
    cm$model_pred$b[tt] = b
    cm$model_pred$calculate(dps)
    out[i] = cm$model_pred$calculate(paste0('Y[', ind, ', 1:20]'))
    i = i + 1
  }
  
  out.keep[t,] <- out 
  
  plot(bs, exp(out), type = 'l',main=t)
  plot(samplesList[[1]][,t],typ='l',main=t,ylim=c(0,150),col=rainbow(3,alpha = 1)[1])
  points(samplesList[[2]][,t],typ='l',col=rainbow(3,alpha = 0.6)[2])
  points(samplesList[[3]][,t],typ='l',col=rainbow(3,alpha = 0.6)[3])
  
}
for(i in 42){
  plot(samplesList[[1]][,i],typ='l',main='SIGMA',
       ylim=c(0,max(c(samplesList[[1]][,i],samplesList[[2]][,i],samplesList[[3]][,i]))),col=rainbow(3,alpha = 1)[1])
  points(samplesList[[2]][,i],typ='l',col=rainbow(3,alpha = 0.6)[2])
  points(samplesList[[3]][,i],typ='l',col=rainbow(3,alpha = 0.6)[3])
}
dev.off()

