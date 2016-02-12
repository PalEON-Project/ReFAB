setwd("/Users/paleolab/babySTEPPS/nimbleRuns/")

library(nimble)
library(splines)

load("min.list.june24.Rdata")

bs_nimble <-  nimbleFunction(
    run = function(u_given = double(0), u = double(1), N0 = double(1),
                    N1 = double(1), N2 = double(1), N3 = double(1)) {
        returnType(double(1))
        
        for(i in 1:7)
            N0[i] = 0
        for(i in 1:7)
            N1[i] = 0
        for(i in 1:6)
            N2[i] = 0
        for(i in 1:5)
            N3[i] = 0

        if(u_given < u[5]){
            N0[4] = 1
        }  else {
            N0[5] = 1
        }

        for(i in 1:7){
            p = 1
            if(N0[i]==0 & N0[i+1]==0){
                N1[i] = 0
            }
            if(N0[i]!=0 & N0[i+1]==0){
                N1[i] = ((u_given - u[i]) / (u[i+p] - u[i])) * N0[i]
            }
            if(N0[i]==0 & N0[i+1]!=0){
                N1[i] = ((u[i+p+1] - u_given) / (u[i+p+1] - u[i+1]) )* N0[i+1]
            }

        }

        for(i in 1:6){
            p = 2
            if(N1[i]==0 & N1[i+1]==0){
                N2[i] = 0
            }
            if(N1[i]!=0 & N1[i+1]!=0){
                N2[i] = ((u_given - u[i]) / (u[i+p] - u[i])) * N1[i] +
                    ((u[i+p+1] - u_given) / (u[i+p+1] - u[i+1]) )* N1[i+1]
            }
            if(N1[i]!=0 & N1[i+1]==0){
                N2[i] = ((u_given - u[i]) / (u[i+p] - u[i])) * N1[i]
            }
            if(N1[i]==0 & N1[i+1]!=0){
                N2[i] = ((u[i+p+1] - u_given) / (u[i+p+1] - u[i+1]) )* N1[i+1]
            }
            
        }

        for(i in 1:5){
            p = 3
            if(N2[i]==0 & N2[i+1]==0){
                N3[i] = 0
            }
            if(N2[i]!=0 & N2[i+1]!=0){
                N3[i] = ((u_given - u[i]) / (u[i+p] - u[i])) * N2[i] +
                    ((u[i+p+1] - u_given) / (u[i+p+1] - u[i+1]) )* N2[i+1]
            }
            if(N2[i]!=0 & N2[i+1]==0){
                N3[i] = ((u_given - u[i]) / (u[i+p] - u[i])) * N2[i]
            }
            if(N2[i]==0 & N2[i+1]!=0){
                N3[i] = ((u[i+p+1] - u_given) / (u[i+p+1] - u[i+1]) )* N2[i+1]
            }
            
        }

        return(N3)
    })


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

#z.seq1 = seq(1,400, length=5)
#Z.knots = bs(biomass,degree = 3,knots=quantile(biomass,c(.5)),intercept=TRUE)

Z.knots = matrix(0,nrow=length(biomass),ncol=5)
u <- c(rep(.0001606565,4), quantile(biomass,c(.5)), rep(190.0118,4))

for(i in 1:length(biomass)){
    u_given <- biomass[i]
	Z.knots[i,] = bs_nimble(u_given, u, rep(0, 8), rep(0, 7), rep(0, 6), rep(0, 5))
}


beta = matrix(NA,ncol(Z.knots),ncol(Y))
p = matrix(NA,nrow(Y),ncol(Y)) ; phi = p
J = nrow(Y)

data = list(Y = counts ,  Z =  Z.knots)

constants = list(n = rowSums(counts), R = ncol(Z.knots), I = ncol(Y), J = nrow(Y))

inits = list(beta = matrix(1,ncol(Z.knots),ncol(Y)),p = matrix(1/20,nrow(Y),ncol(Y)))

dimensions = list(exp.phi = dim(phi), phi.first = dim(phi), Z = dim(Z.knots), beta = dim(beta), p = dim(p), Y = dim(counts), n = nrow(Y))

# in BUGS code, to calculate the vector of basis matrix values for a given biomass, pass that biomass in as 'u_given', pass in the vector of u values for the knots and pass in N0,N1,N2,N3 of correct length - you can do this simply by providing N0,N1,N2,N3 as part of the 'constants' argument given to the 'nimbleModel' function


if(FALSE){
model <- nimbleModel(code, inits = inits, constants = constants, data = data, dimensions = dimensions)

# compiled version of the model
Cmodel <- compileNimble(model)

# set up MCMC
spec <- configureMCMC(model, thin = 10, print = TRUE)
spec$addMonitors(c('beta')) 

# set up monitoring of whatever
# model variables you want posterior samples for - by default, top level
# parameters are already included, so 'mu' in the above example would by
# default be monitored. 'psi' and 'theta' are just for illustration -
# obviously they are not part of my toy model above

# create MCMC algorithm for the model
Rmcmc <- buildMCMC(spec)
# compiled version of the MCMC
Cmcmc <- compileNimble(Rmcmc, project = model)
# 
# run MCMC for 5000 iterations
Cmcmc$run(5000)
samples1 <- as.matrix(Cmcmc$mvSamples)
}

# save(samples,file = "nimble.betas.Rdata")
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
  
  sigma ~ dunif(0,1) #GELMAN PAPER
  
  for(p in 1:P){
  	for(t in 1:T){
  	   log(b[p,t]) ~ dnorm(log(b[p,t-1]),1/sigma^2)
    }
  }
  
  
  for(p in 1:P){
  	for(t in 1:T){
    Zb[p,t,1:5] <- bs_nimble(b[p,t], u[1:8], N0[1:8], N1[1:7], N2[1:6], N3[1:5])
    }
    }

  for(i in 1:I){
  	for(p in 1:P){
  	for(t in 1:T){
  		phi.first[p,t,i] <- sum(Zb[p,t,1:5] %*% beta[1:5,i])
  	}
  }
  }
 # phi.first[,] <- Zb[,]%*%beta[,]
  
  for(p in 1:P){
  	for(t in 1:T){
    for(i in 1:I){
      exp.phi[p,t,i] <- exp(phi.first[p,t,i])
    }
  }
  }
  
  for(p in 1:P){
  	for(t in 1:T){
   # for(i in 1:10){
   #  phi[j,i] <- exp.phi[j,i]/sum(exp.phi[j,])
   #} 
   #p[j,] ~ ddirch(exp.phi[j,]) #http://sourceforge.net/p/mcmc-jags/discussion/610037/thread/c21ef62a
   #Y[j,] ~ dmulti(prob = p[j,], size = n[j])
   Y[j,] ~ ddirchmulti(exp.phi[p,t,],n[p,t])
  }
  }
  
})


J = nrow(counts)#nrow(Z.new)
Zb = matrix(NA,J,ncol(Z.knots))
DFS = ncol(Zb)
phi = matrix(NA,J,ncol(counts)); phi.first = phi;
#beta.est.sim = matrix(summary(csamp.sim.cal)$statistics[,1],4,ncol(Y))
#load("beta.samps.Rdata")
beta.est = matrix(colMeans(samples[100:nrow(samples),]),ncol(Z.knots),ncol(Y)) #matrix(summary(csamp.real.cal)$statistics[,1],ncol(Z.knots),ncol(Y))# 
new.biomass = seq(1,190,1)
#Z.new = bs(new.biomass,intercept=TRUE,knots = attr(Z.knots,"knots"),Boundary.knots = attr(Z.knots,"Boundary.knots"))
Z.new = matrix(0,nrow=length(new.biomass),ncol=5)
u <- c(rep(0,4), 26, rep(190,4))

for(i in 1:length(new.biomass)){
    u_given <- new.biomass[i]
	Z.new[i,] = bs_nimble(u_given, u, rep(0, 8), rep(0, 7), rep(0, 6), rep(0, 5))
}

data.pred = list(Y = counts)

constants.pred = list(beta = beta.est, I = ncol(counts), DFS = DFS, J = J, n = rowSums(counts),  Z =  Z.new, u = u, N0 = rep(0, 8), N1 = rep(0, 7), N2 = rep(0, 6), N3 = rep(0, 5))

inits.pred = list(b=rep(100,J))

dimensions.pred = list(exp.phi = dim(phi), phi.first = dim(phi), Zb = dim(Zb), beta = dim(beta.est), Y = dim(counts))

model_pred <- nimbleModel(pred_code, inits = inits.pred, constants = constants.pred, data = data.pred, dimensions = dimensions.pred)

#Cmodel.pred <- compileNimble(model_pred) #opens browser...
spec.pred <- configureMCMC(model_pred, thin = 1, print = TRUE)
spec.pred$addMonitors(c('b')) 
Rmcmc.pred <- buildMCMC(spec.pred)
cm <- compileNimble(model_pred)
Cmcmc.pred <- compileNimble(Rmcmc.pred, project = model_pred) #Error in cModel$.nodeValPointers_byGID : $ operator not defined for this S4 class

ptm <- proc.time()
Cmcmc.pred$run(5000)
samples.pred <- as.matrix(Cmcmc.pred$mvSamples)
proc.time() - ptm

# # #samples.pred1<-samples.pred
# load("samples.pred1.Rdata")
# save(samples.pred,file="samples.pred.Rdata")

# plot(colMeans(samples.pred),colMeans(samples.pred1))
# plot(biomass,colMeans(samples.pred),xlim=c(0,200),ylim=c(0,200),pch=19,xlab="True Biomass",ylab="Predicted Mean Biomass")
# abline(a=0,b=1)

# pdf("pred_validation.pdf")
# plot(biomass,colMeans(samples.pred),xlim=c(0,400),ylim=c(0,400))
# abline(a=0,b=1)
# dev.off()
