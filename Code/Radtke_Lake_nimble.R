library(nimble)
library(splines)
library(maps)
library(plyr)
library(oce)
library(RCurl)

load("add.bacon.Rdata")
load("2016-05-31nimble.betas.Rdata")
source("/Users/paleolab/babySTEPPS/Code/bs_nimble.R")

# set up the "d" function for the distribution
ddirchmulti <- nimbleFunction(
  run = function(x = double(1), alpha = double(1), size = double(0), log_value = integer(0)){
  returnType(double(0))
  logProb <- lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum(lgamma(alpha + x)) - lgamma(sum(alpha) + size)
    
    if(log_value) {
      return(logProb)
    } else {
      return(exp(logProb))
    }
    
   }
)

# set up the "r" function
rdirchmulti <- nimbleFunction(
run = function(n = integer(0), alpha = double(1), size = double(0)) {
returnType(double(1))
if(n != 1) nimPrint("rdirchmulti only allows n = 1; using n = 1.")
p <- rdirch(1, alpha)
return(rmulti(1, size = size, prob = p))
})

# tell NIMBLE about the newly available distribution
registerDistributions(list(ddirchmulti = list(BUGSdist = "ddirchmulti(alpha, size)", 
  types = c('value = double(1)', 'alpha = double(1)'))))
    
#####
##### BUGS CODE #####
#####
    
pred_code <- nimbleCode({
	
  sigma ~ dunif(0,5) #GELMAN PAPER
  
  b[1,1] ~ dunif(0,157)
  
  for(t in 2:T){
  	   b[1,t] ~ T(dlnorm(log(b[1,t-1]),1/sigma^2),0,157)
    }
 
  	for(t in 1:T){
      Zb[t,1:5] <- bs_nimble(b[1,t], u[1:3], N0[1:2], N1[1:3], N2[1:4], N3[1:5])
    }

  for(i in 1:I){
  	  for(t in 1:T){
  		phi.first[t,i] <- sum(Zb[t,1:5] %*% beta[1:5,i])
  		exp.phi[t,i] <- exp(phi.first[t,i])
  	  }
  }

  	for(j in 1:J){
        Y[j,1:20] ~ ddirchmulti(exp.phi[age.index[j,1], 1:20], n[j])
     }
     
  
})

#####
##### SETUP #####
#####

x = new.pol1[new.pol1$age_bacon>=200,]
x = x[x$age_bacon<=10000,]

x.meta = x[,c('SiteID','LatitudeNorth',"LongitudeWest","dataset.id","sitename","age_bacon")]
ten.count = x[,7:26]

ten.count.save = ten.count
ten.count = round(ten.count.save)

site_number = unique(x.meta[x.meta$sitename=='RadtkeLake',1])
ten.count.use = ten.count[which(x.meta$SiteID==site_number),]

Y = as.matrix(ten.count.use)

sample.ages <- x.meta[x.meta[,1]==site_number,]$age_bacon
age.bins <- seq(0,10000,100)
age.index <- as.matrix(as.numeric(cut(sample.ages,breaks = age.bins,labels=seq(1:(length(age.bins)-1)))))

Z.knots = Z.new
T = length(age.bins)-1
I = ncol(Y)
K = ncol(Z.knots)
J = length(age.index)
n = rowSums(Y)
Zb = matrix(NA,T,K)
phi.first = matrix(NA,T,I); exp.phi = phi.first
beta.est = matrix(colMeans(samples1[300:nrow(samples1),]),K,I) 
new.biomass = seq(1,200,1)
Z.new = matrix(0,nrow=length(new.biomass),ncol=K)

u<-c(rep(attr(Z.knots,"Boundary.knots")[1],1),attr(Z.knots,"knots"),rep(attr(Z.knots,"Boundary.knots")[2],1))

data.pred = list(Y = Y)

constants.pred = list(beta = beta.est, I = I, J=J, T=T, n = n, u = u, N0 = rep(0, (length(u)-1)),N1 = rep(0, (length(u))), N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)),age.index = age.index)

inits.pred = list(b=matrix(100,1,T), sigma = 4.5)

dimensions.pred = list(exp.phi = dim(exp.phi), phi.first = dim(exp.phi), Zb = dim(Zb), beta = dim(beta.est), Y = dim(Y),  b = dim(inits.pred$b))

#####
##### MODEL FUNCTIONS #####
#####

set.seed(0)
model_pred <- nimbleModel(pred_code, inits = inits.pred, constants = constants.pred, data = data.pred, dimensions = dimensions.pred)

spec.pred <- configureMCMC(model_pred, thin = 10, print = FALSE)
spec.pred$getSamplers()
spec.pred$addMonitors(c('b')) 
Rmcmc.pred <- buildMCMC(spec.pred)
cm <- compileNimble(model_pred)

Cmcmc.pred <- compileNimble(Rmcmc.pred, project = model_pred) 
Cmcmc.pred$run(20000)

samples.pred <- as.matrix(Cmcmc.pred$mvSamples)

#####
##### TRACE PLOTS #####
#####

par(mfrow=c(3,1))
for(i in 40:42) {
	plot(samples.pred[,i],typ='l',ylab=i)
	abline(h=157)
	}
	
#####
##### POLLEN PLOTS #####
#####	
 par(mfrow=c(3,3))
 for(i in 1:20) { 
 	plot(sample.ages,Y[,i],main=colnames(Y)[i])
 	 abline(v=3000)
 	 }
 	 
#####
##### TIME SERIES #####
#####
plot(apply(samples.pred[500:2000,1:99],2,quantile,c(0.025,0.5,0.975))[2,],typ='l',ylim=c(0,200))
rug(sample.ages/100)

dyn.unload(model_pred$nimbleProject$cppProjects[[1]]$getSOName())
