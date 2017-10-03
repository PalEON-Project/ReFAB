setwd("/Users/paleolab/babySTEPPS/")

library(nimble)
library(splines)

load("add.bacon2.Rdata")
source(paste0('~/babystepps/Code/',"bs_nimble.R"))

if(length(which(colnames(Y)=='PINUSX'))>0){
	Y = Y[,-3]
	counts = counts[,-3]
}

code <- nimbleCode({
  
  #delta ~ dunif(0,1000)
  
  for(r in 1:R){ 
    for(i in 1:I){
      beta[r,i] ~ dnorm(0,.04)
    }
  }
  
  phi.first[,] <- Z[,]%*%beta[,]
  
  for(j in 1:J){
    for(i in 1:I){
      exp.phi[j,i] <- exp(phi.first[j,i])
    }
    row.sums[j] <- sum(exp.phi[j,])
    # 	  for(i in 1:10){
    # 	    phi[j,i] <- exp.phi[j,i]/(row.sums[j])
    # 	  }
  }
  
  for(j in 1:J){
    p[j,] ~ ddirch(exp.phi[j,]) #http://sourceforge.net/p/mcmc-jags/discussion/610037/thread/c21ef62a
    Y[j,] ~ dmulti(size = n[j], prob = p[j,])
  }
  
})

if(FALSE){
Z.knots = matrix(0,nrow=length(biomass),ncol=(length(u)+2));
Z.knots = bs(biomass,intercept=TRUE,df=5)
# 0 to 150 grid points rows are biomass and each column is basis function
# 5 basis functions 
#plot emp props of key taxa based on sampling dates by site
#get a sense of the raw data how long it took for the transition

beta = matrix(NA,ncol(Z.knots),ncol(Y))
p = matrix(NA,nrow(Y),ncol(Y)) ; phi = p
J = nrow(Y)

data = list(Y = as.matrix(counts) ,  Z =  Z.knots)

constants = list(n = rowSums(counts), R = ncol(Z.knots), I = ncol(Y), J = nrow(Y))

inits = list(beta = matrix(1,ncol(Z.knots),ncol(Y)),p = matrix(1/20,nrow(Y),ncol(Y)))

dimensions = list(exp.phi = dim(phi), phi.first = dim(phi), Z = dim(Z.knots), beta = dim(beta), p = dim(p), Y = dim(counts), n = nrow(Y))

# in BUGS code, to calculate the vector of basis matrix values for a given biomass, pass that biomass in as 'u_given', pass in the vector of u values for the knots and pass in N0,N1,N2,N3 of correct length - you can do this simply by providing N0,N1,N2,N3 as part of the 'constants' argument given to the 'nimbleModel' function

model <- nimbleModel(code, inits = inits, constants = constants, data = data, dimensions = dimensions)

# compiled version of the model
Cmodel <- compileNimble(model)

# set up MCMC
spec <- configureMCMC(model, thin = 10, print = TRUE)
spec$addMonitors(c('beta','row.sums')) 

# set up monitoring of whatever
# model variables you want posterior samples for - by default, top level
# parameters are already included, so 'mu' in the above example would by
# default be monitored. 'psi' and 'theta' are just for illustration -
# obviously they are not part of my toy model above

# create MCMC algorithm for the model
Rmcmc <- buildMCMC(spec)

# compiled version of the MCMC
Cmcmc <- compileNimble(Rmcmc, project = model)

# run MCMC for 5000 iterations
Cmcmc$run(2000)
samples1 <- as.matrix(Cmcmc$mvSamples)

beta.est.real = matrix(colMeans(samples1[100:nrow(samples1),grep('beta',colnames(samples1))]),ncol(Z.knots),ncol(Y))

plot.betas <- as.matrix(exp(Z.knots%*%beta.est.real)/rowSums(exp(Z.knots%*%beta.est.real)))

pdf(paste0(fig.dir,paste0("scatter.plus.betas",Sys.Date(),".pdf")))
par(mfrow=c(3,3))
for(i in 1:ncol(counts)){
	plot(biomass,counts[,i]/total_counts,pch=19,cex=.7,col='black',ylab="Pollen Proportions",main=colnames(counts)[i],xlab="Biomass")
	points(biomass,plot.betas[,i],col='red',pch=19,cex=1)
	abline(v=u)
  }
dev.off()

row.sums.est = colMeans(samples1[100:nrow(samples1),grep('row.sums',colnames(samples1))])
pdf('calib.lik.pdf')
plot(biomass,row.sums.est,pch=19)
dev.off()
for(i in 1:ncol(counts)){
	plot(samples1[,i],type="l")
  }

save(samples1,file = paste0("nimble.betas",Sys.Date(),".Rdata"))
}
load("nimble.betas.Rdata")

# samples has rows as iterations and columns as variables, you'll have
# to manually remove a burn-in period
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

pred_code <- nimbleCode({
   for(j in 1:J){
    b[j] ~ dunif(0,156)
 
    Zb[j,1:5] <- bs_nimble(b[j], u[1:3], N0[1:2], N1[1:3], N2[1:4], N3[1:5])
    }

  for(i in 1:I){
  	for(j in 1:J){
  		phi.first[j,i] <- sum(Zb[j,1:5] %*% beta[1:5,i])
  	}
  }
  
  for(j in 1:J){
    for(i in 1:I){
      exp.phi[j,i] <- exp(phi.first[j,i])
    }
  }
  
  for(j in 1:J){
   Y[j,] ~ ddirchmulti(exp.phi[j,],n[j])
  }
  
})

J = nrow(counts)#nrow(Z.new)
Zb = matrix(NA,J,ncol(Z.knots))
DFS = ncol(Zb)
phi = matrix(NA,J,ncol(counts)); phi.first = phi;
#beta.est.sim = matrix(summary(csamp.sim.cal)$statistics[,1],4,ncol(Y))
#load("beta.samps.Rdata")
beta.est = matrix(colMeans(samples1[100:nrow(samples1),]),ncol(Z.knots),ncol(Y)) #matrix(summary(csamp.real.cal)$statistics[,1],ncol(Z.knots),ncol(Y))# 
new.biomass = seq(1,u[length(u)],1)
#Z.new = bs(new.biomass,intercept=TRUE,knots = attr(Z.knots,"knots"),Boundary.knots = attr(Z.knots,"Boundary.knots"))
Z.new = matrix(0,nrow=length(new.biomass),ncol=5)
u <- u #should have defined knots in calibration
u<-c(rep(attr(Z.knots,"Boundary.knots")[1],1),attr(Z.knots,"knots"),rep(attr(Z.knots,"Boundary.knots")[2],1))

for(i in 1:length(new.biomass)){
    u_given <- new.biomass[i]
	Z.new[i,] = bs_nimble(u_given, u=u, N0 = rep(0, (length(u)-1)),N1 = rep(0, (length(u))), N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))
}

data.pred = list(Y = as.matrix(counts))

constants.pred = list(beta = beta.est, I = ncol(counts), DFS = DFS, J = J, n = rowSums(counts),  Z =  Z.new, u = u, N0 = rep(0, (length(u)-1)),N1 = rep(0, (length(u))), N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))

inits.pred = list(b=rep(100,J))

dimensions.pred = list(exp.phi = dim(phi), phi.first = dim(phi), Zb = dim(Zb), beta = dim(beta.est), Y = dim(counts))

model_pred <- nimbleModel(pred_code, inits = inits.pred, constants = constants.pred, data = data.pred, dimensions = dimensions.pred)

cm <- compileNimble(model_pred)
vals <- 1:157
outLik = outPost = matrix(NA, 157, J)
for(j in 1:J){
	calcNodes <-  cm$getDependencies(paste0('b[',j,']'))
for(val in vals) {
    cm$b[j] <- val
    outPost[val,j] = calculate(cm,calcNodes)# cm$calculate(calcNodes)
    # likelihood portion
    outLik[val,j] =  calculate(cm,calcNodes[grep("Y", calcNodes)]) # cm$calculate(calcNodes[45])  #
}	
}

pdf('calib_dataset_max_liks.pdf')
par(mfrow=c(3,3))
   for(j in 1:J){
        plot(vals,exp(outLik[,j] - max(outLik[,j]))/-sum(outLik[,j])
        ,typ='l',ylab=NA,main=w[j])
    }
dev.off()

#Cmodel.pred <- compileNimble(model_pred) #opens browser...
spec.pred <- configureMCMC(model_pred, thin = 10, print = TRUE)
spec.pred$addMonitors(c('b')) 
Rmcmc.pred <- buildMCMC(spec.pred)

Cmcmc.pred <- compileNimble(Rmcmc.pred, project = model_pred) #Error in cModel$.nodeValPointers_byGID : $ operator not defined for this S4 class

ptm <- proc.time()
Cmcmc.pred$run(5000)
samples.pred <- as.matrix(Cmcmc.pred$mvSamples)
proc.time() - ptm

vals <- 1:145
outLik = outPost = matrix(NA, 145, J)
for(j in 1:J){
	calcNodes <-  cm$getDependencies(paste0('b[',j,']'))
for(val in vals) {
    cm$b[j] <- val
    outPost[val,j] = calculate(cm,calcNodes)# cm$calculate(calcNodes)
    # likelihood portion
    outLik[val,j] =  calculate(cm,calcNodes[grep("Y", calcNodes)]) # cm$calculate(calcNodes[45])  #
}	
}

pdf('calib_dataset_max_liks.pdf')
par(mfrow=c(3,3))
   for(j in 1:J){
        plot(vals,exp(outLik[,j] - max(outLik[,j]))/-sum(outLik[,j])
        ,typ='l',ylab=NA,main=w[j])
    }
dev.off()

par(mfrow=c(2,2))
for(i in 1:ncol(samples.pred)){
	plot(samples.pred[,i])
}

# #samples.pred1<-samples.pred
load("samples.pred1.Rdata")
save(samples.pred,file="twothirds.pred.Rdata")

pdf(paste0("pred_validation",Sys.Date(),".pdf"))
plot(biomass,colMeans(samples.pred[100:nrow(samples.pred),]),xlim=c(0,200),ylim=c(0,200),pch=19,xlab="True Biomass",ylab="Predicted Mean Biomass")
abline(a=0,b=1)
abline(lm(biomass~colMeans(samples.pred[100:nrow(samples.pred),])+0),lty=2)
mtext(paste("r-squared",summary(lm(biomass~colMeans(samples.pred[100:nrow(samples.pred),])+0))$r.squared))
dev.off()
