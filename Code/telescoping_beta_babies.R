setwd("/Users/paleolab/babySTEPPS/")

library(nimble)
library(splines)

load("add.bacon2.Rdata")
model.dir <- c('/Users/paleolab/babySTEPPS/Code/')
source(paste0(model.dir,"bs_nimble.R"))

code <- nimbleCode({
  
  for(r in 1:R){ 
    for(i in 1:I){
      beta[r,i] ~ dnorm(0,.04)
      beta.pine[r,i] ~ dnorm(0,.04)
    }  
  }
  
  phi.first[,] <- Z[,]%*%beta[,]
  pine.phi[,] <- Z[,]%*%beta.pine[,]
    
  for(j in 1:J){
    for(i in 1:I){
      exp.phi[j,i] <- exp(phi.first[j,i])
      exp.pine.phi[j,i] <- exp(pine.phi[j,i])
    }
    row.sums[j] <- sum(exp.phi[j,])
  }
     
    for(j in 1:J){
    	p.true[j,1] ~ dbeta(exp.phi[j,1],exp.pine.phi[j,1])
    	p.sum[j,1] <- 1
    	
    for(i in 2:(I-1)){
        p.rel[j,i]  ~ dbeta(exp.phi[j,i],exp.pine.phi[j,i])
        p.sum[j,i] <- sum(p.true[j,1:(i-1)])
        p.true[j,i] <-  p.rel[j,i] * (1 - p.sum[j,i])
    }	
       p.sum[j,21] <- sum(p.true[j,1:20])
       p.true[j,21] <- 1 - p.sum[j,21]
    }  
       
  for(j in 1:J){
    Y[j,] ~ dmulti(size = n[j], prob = p.true[j,])
  }
  
})

counts <- Y[,rev(order(colMeans(Y)))]

Z.knots = matrix(0,nrow=length(biomass),ncol=(length(u)+2));
Z.knots = bs(biomass,intercept=TRUE,df=5)
# 0 to 150 grid points rows are biomass and each column is basis function
# 5 basis functions 
#plot emp props of key taxa based on sampling dates by site
#get a sense of the raw data how long it took for the transition

beta = matrix(NA,ncol(Z.knots),ncol(Y))
beta.pine = matrix(NA,ncol(Z.knots),1)
p = matrix(NA,nrow(Y),ncol(Y)) ; phi = p
J = nrow(Y)

data = list(Y = as.matrix(counts) ,  Z =  Z.knots)

constants = list(n = rowSums(counts), R = ncol(Z.knots), I = ncol(Y), J = nrow(Y))

inits = list(beta = matrix(1, ncol(Z.knots), ncol(Y)), beta.pine = matrix(1, ncol(Z.knots), ncol(Y)), p.rel = matrix(1/21, nrow(Y), ncol(Y)))

dimensions = list(exp.phi = dim(phi), exp.pine.phi = dim(phi), phi.first = dim(phi), pine.phi = dim(phi), Z = dim(Z.knots), beta = dim(beta), beta.pine = dim(beta), p.true = dim(p), Y = dim(counts), n = nrow(Y), pine.dirch = dim(phi),p.rel = dim(p),p.sum = dim(p))

# in BUGS code, to calculate the vector of basis matrix values for a given biomass, pass that biomass in as 'u_given', pass in the vector of u values for the knots and pass in N0,N1,N2,N3 of correct length - you can do this simply by providing N0,N1,N2,N3 as part of the 'constants' argument given to the 'nimbleModel' function

model <- nimbleModel(code, inits = inits, constants = constants, data = data, dimensions = dimensions)

# compiled version of the model
Cmodel <- compileNimble(model)

# set up MCMC
spec <- configureMCMC(model, thin = 10, print = TRUE)
spec$addMonitors(c('beta','beta.pine','p.true','p.rel')) 

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
samples.mixed <- as.matrix(Cmcmc$mvSamples)

par(mfrow=c(3,3))
for(i in grep('p.rel',colnames(samples.mixed))[1:10]){
	plot(samples.mixed[,i],main=colnames(samples.mixed)[i],type="l")
  }

beta1.est.real = matrix(colMeans(samples.mixed[100:nrow(samples.mixed),1:105]),ncol(Z.knots),ncol(Y))
beta2.est.real = matrix(colMeans(samples.mixed[100:nrow(samples.mixed),106:210]),ncol(Z.knots),ncol(Y))

plot.betas1 <- as.matrix(exp(Z.knots%*%beta1.est.real)/rowSums(exp(Z.knots%*%beta1.est.real)))
plot.betas2 <- as.matrix(exp(Z.knots%*%beta2.est.real)/rowSums(exp(Z.knots%*%beta2.est.real)))


beta1.est.real <- array(NA, dim=c(ncol(Z.knots),ncol(Y),200))
beta2.est.real <- array(NA, dim=c(ncol(Z.knots),ncol(Y),200))
plot.betas1 <- array(NA, dim=c(length(biomass),ncol(Y),200))
plot.betas2 <- array(NA, dim=c(length(biomass),ncol(Y),200))

for(i in 1:200){
	beta1.est.real[,,i] = matrix(samples.mixed[i,1:105],ncol(Z.knots),ncol(Y))
    beta2.est.real[,,i] = matrix(samples.mixed[i,106:210],ncol(Z.knots),ncol(Y))

    plot.betas1[,,i] <- as.matrix(exp(Z.knots%*%beta1.est.real[,,i])/rowSums(exp(Z.knots%*%beta1.est.real[,,i])))
    plot.betas2[,,i] <- as.matrix(exp(Z.knots%*%beta2.est.real[,,i])/rowSums(exp(Z.knots%*%beta2.est.real[,,i])))
    
}

prop.quants <- array(NA,dim=c(21,101,3))

for(s in 1:21){
	for(b in 1:101){
		prop.quants[s,b,] <- quantile(plot.betas1[b,s,],c(.025,.5,.975))
	}
}

par(mfrow=c(3,3))
for(i in 1:ncol(counts)){
    plot(biomass,prop.quants[i,,1],ylim=c(0,1))
    for(r in 2:3){
    	points(biomass,prop.quants[i,,r])
    }
}

matrix(samples.mixed[i,grep('p.new[1,',colnames(samples.mixed))],J,21,byrow=TRUE)

pdf(paste0(fig.dir,paste0("tele.betas",Sys.Date(),".pdf")))
par(mfrow=c(3,3))
for(i in 1:ncol(counts)){
	plot(biomass,counts[,i]/total_counts,pch=19,cex=.7,col='black',ylab="Pollen Proportions",main=colnames(counts)[i],xlab="Biomass",ylim=c(0,max(c(plot.betas1[,i],plot.betas2[,i],counts[,i]/total_counts))))
	points(biomass,plot.betas1[,i],col='red',pch=19,cex=1)
	points(biomass,plot.betas2[,i],col='blue',pch=19,cex=1)
	abline(v=u)
  }
dev.off()

grep('p',colnames(samples.mixed))

par(mfrow=c(3,3))
for(i in grep('p.new',colnames(samples.mixed))){
	plot(samples.mixed[,i],main=colnames(samples.mixed)[i],type="l")
  }

save(samples.mixed,file = paste0("nimble.betas_1_2",Sys.Date(),".Rdata"))

load("nimble.betas_1_2.Rdata")

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
  		phi.first1[j,i] <- sum(Zb[j,1:5] %*% beta1[1:5,i])
  	}
  }
  
  for(j in 1:J){
    for(i in 1:I){
      exp.phi[j,i] <- exp(phi.first[j,i])
      exp.phi1[j,i] <- exp(phi.first1[j,i])
    }
  }
  
   for(j in 1:J){
    	p[j,1] ~ dbeta(exp.phi[j,1],exp.phi1[j,1])
    	
    for(i in 2:I){
        p.new[j,i]  ~ dbeta(exp.phi[j,i],exp.phi1[j,i])
        p[j,i] <-  p.new[j,i] * (1 - p[j,i-1]) 
    }	
    }    
    
  for(j in 1:J){
   Y[j,] ~ dmulti(p[j,],n[j])
  }
  
})

J = nrow(counts)#nrow(Z.new)
Zb = matrix(NA,J,ncol(Z.knots))
DFS = ncol(Zb)
phi = matrix(NA,J,ncol(counts)); phi.first = phi;
#beta.est.sim = matrix(summary(csamp.sim.cal)$statistics[,1],4,ncol(Y))
#load("beta.samps.Rdata")
beta.est = matrix(colMeans(samples.mixed[100:nrow(samples.mixed),]),ncol(Z.knots),ncol(Y)) #matrix(summary(csamp.real.cal)$statistics[,1],ncol(Z.knots),ncol(Y))# 
new.biomass = seq(1,u[length(u)],1)
#Z.new = bs(new.biomass,intercept=TRUE,knots = attr(Z.knots,"knots"),Boundary.knots = attr(Z.knots,"Boundary.knots"))
Z.new = matrix(0,nrow=length(new.biomass),ncol=5)
u <- u #should have defined knots in calibration
u<-c(rep(attr(Z.knots,"Boundary.knots")[1],1),attr(Z.knots,"knots"),rep(attr(Z.knots,"Boundary.knots")[2],1))

ro = mean(samples.mixed[,'ro'])

for(i in 1:length(new.biomass)){
    u_given <- new.biomass[i]
	Z.new[i,] = bs_nimble(u_given, u=u, N0 = rep(0, (length(u)-1)),N1 = rep(0, (length(u))), N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))
}

data.pred = list(Y = as.matrix(counts))

constants.pred = list(beta = beta1.est.real, beta1 = beta2.est.real, I = ncol(counts), DFS = DFS, J = J, n = rowSums(counts), u = u, N0 = rep(0, (length(u)-1)),N1 = rep(0, (length(u))), N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))

inits.pred = list(b=rep(100,J))

dimensions.pred = list(exp.phi = dim(phi), exp.phi1 = dim(phi), phi.first = dim(phi), phi.first1 = dim(phi), Zb = dim(Zb), beta = dim(beta.est), beta1 = dim(beta.est), Y = dim(counts), p = dim(phi))

model_pred <- nimbleModel(pred_code, inits = inits.pred, constants = constants.pred, data = data.pred, dimensions = dimensions.pred)

#Cmodel.pred <- compileNimble(model_pred) #opens browser...
spec.pred <- configureMCMC(model_pred, thin = 10, print = TRUE)
spec.pred$addMonitors(c('b')) 
Rmcmc.pred <- buildMCMC(spec.pred)
cm <- compileNimble(model_pred)
Cmcmc.pred <- compileNimble(Rmcmc.pred, project = model_pred) #Error in cModel$.nodeValPointers_byGID : $ operator not defined for this S4 class

ptm <- proc.time()
Cmcmc.pred$run(5000)
samples.pred <- as.matrix(Cmcmc.pred$mvSamples)
proc.time() - ptm

par(mfrow=c(2,2))
for(i in 1:ncol(samples.pred)){
	plot(samples.pred[,i])
}

# #samples.pred1<-samples.pred
load("samples.pred1.Rdata")
save(samples.pred,file="twothirds.pred.Rdata")

pdf(paste0("pred_validation_two_betas",Sys.Date(),".pdf"))

plot(biomass,colMeans(samples.pred[100:nrow(samples.pred),grep('b',colnames(samples.pred))]),xlim=c(0,200),ylim=c(0,200),pch=19,xlab="True Biomass",ylab="Predicted Mean Biomass")
abline(a=0,b=1)
abline(lm(biomass~colMeans(samples.pred[100:nrow(samples.pred),grep('b',colnames(samples.pred))])+0),lty=2)
mtext(paste("r-squared",summary(lm(biomass~colMeans(samples.pred[100:nrow(samples.pred),grep('b',colnames(samples.pred))])+0))$r.squared))

dev.off()
