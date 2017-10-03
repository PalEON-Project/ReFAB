library(nimble)
library(splines)
load("babySTEPPS/min.list.june24.Rdata")

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
Z.knots = bs(biomass,degree = 3,knots=quantile(biomass,c(.5)),intercept=TRUE)

beta = matrix(NA,ncol(Z.knots),ncol(Y))
p = matrix(NA,nrow(Y),ncol(Y)) ; phi = p
J = nrow(Y)

data = list(Y = counts ,  Z =  Z.knots)

constants = list(n = rowSums(counts), R = ncol(Z.knots), I = ncol(Y), J = nrow(Y))

inits = list(beta = matrix(1,ncol(Z.knots),ncol(Y)),p = matrix(1/20,nrow(Y),ncol(Y)))

dimensions = list(exp.phi = dim(phi), phi.first = dim(phi), Z = dim(Z.knots), beta = dim(beta), p = dim(p), Y = dim(counts), n = nrow(Y))

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
samples <- as.matrix(Cmcmc$mvSamples)

save(samples,file = "nimble.betas.Rdata")
load("nimble.betas.Rdata")

# samples has rows as iterations and columns as variables, you'll have
# to manually remove a burn-in period

pred_code <- nimbleCode({
  
  for(j in 1:J){
    b[j] ~ dunif(1, 400)
    b_trunc <- b[j]
 
    Zb[j,1:5] <- bs_nimble(b_trunc, u[1:8], N0[1:8], N1[1:7], N2[1:6], N3[1:5])
    }

  for(i in 1:I){
  	for(j in 1:J){
  		phi.first[j,i] <- sum(Zb[j,1:5] %*% beta[1:5,i])
  	}
  }
  #phi.first[,] <- Zb[,]%*%beta[,]
  
  for(j in 1:J){
    for(i in 1:I){
      exp.phi[j,i] <- exp(phi.first[j,i])
    }
  }
  
  for(j in 1:J){
   # for(i in 1:10){
   #  phi[j,i] <- exp.phi[j,i]/sum(exp.phi[j,])
   #} 
    p[j,] ~ ddirch(exp.phi[j,] + .025) #http://sourceforge.net/p/mcmc-jags/discussion/610037/thread/c21ef62a
    Y[j,] ~ dmulti(prob = p[j,], size = n[j])
  }
  
})

# x = pol.cal.count[pol.cal.count$Age>=200,]
# x = x[x$Age<=1000,]

# x = x[,-which(colnames(x)==c("PINUSX"))]
# trees <- c("ACERX","CUPRESSA","FRAXINUX","FAGUS","CYPERACE","LARIXPSEU","TSUGAX","QUERCUS","TILIA","BETULA","PICEAX","OSTRYCAR","ULMUS","ABIES","POPULUS")
# other.trees <- c("TAXUS","NYSSA","JUGLANSX","CASTANEA","PLATANUS","SALIX","LIQUIDAM","ALNUSX")
# ten.count = matrix(0,nrow(x),length(trees)+3)
# prairie <- c("CORYLUS","ARTEMISIA","ASTERX","POACEAE","AMBROSIA","CHENOAMX")
# ten.count[,1] <- rowSums(x[,prairie])
# ten.count[,2] <- rowSums(x[,other.trees])
# ten.count[,3:(length(trees)+2)] <- as.matrix(x[,trees])
# ten.count[,(length(trees)+3)] <- rowSums(x[,7:ncol(x)]) - rowSums(ten.count)
# colnames(ten.count)<-c("prairie","other trees",trees,"other herbs")

# ten.count <- round(ten.count)
# ten.count <- ten.count[1:2,]

J = nrow(counts)#nrow(Z.new)
Zb = matrix(NA,J,ncol(Z.knots))
DFS = ncol(Zb)
p = matrix(NA,J,ncol(counts)); phi.first = p; phi = p
#beta.est.sim = matrix(summary(csamp.sim.cal)$statistics[,1],4,ncol(Y))
#load("beta.samps.Rdata")
beta.est = matrix(colMeans(samples[100:nrow(samples),]),ncol(Z.knots),ncol(Y)) #matrix(summary(csamp.real.cal)$statistics[,1],ncol(Z.knots),ncol(Y))# 
new.biomass = seq(1,156,1)
Z.new = bs(new.biomass,intercept=TRUE,knots = attr(Z.knots,"knots"),Boundary.knots = attr(Z.knots,"Boundary.knots"))
u <- c(rep(0,4), 26, rep(190,4))

data.pred = list(Y = counts)

constants.pred = list(beta = beta.est, I = ncol(counts), DFS = DFS, J = J, n = rowSums(counts),  Z =  Z.new, u = u, N0 = rep(0, 8), N1 = rep(0, 7), N2 = rep(0, 6), N3 = rep(0, 5))

inits.pred = list(b=rep(300,J), p = matrix(1/20,nrow(Y),ncol(Y)))

dimensions.pred = list(exp.phi = dim(phi), phi.first = dim(phi), p = dim(p), Zb = dim(Zb), beta = dim(beta.est), Y = dim(counts))

model_pred <- nimbleModel(pred_code, inits = inits.pred, constants = constants.pred, data = data.pred, dimensions = dimensions.pred)

#Cmodel.pred <- compileNimble(model_pred) #opens browser...
spec.pred <- configureMCMC(model_pred, thin = 10, print = TRUE)
spec.pred$addMonitors(c('b')) 
Rmcmc.pred <- buildMCMC(spec.pred)
cm <- compileNimble(model_pred)
Cmcmc.pred <- compileNimble(Rmcmc.pred, project = model_pred) #Error in cModel$.nodeValPointers_byGID : $ operator not defined for this S4 class

ptm <- proc.time()
Cmcmc.pred$run(100)
samples.pred <- as.matrix(Cmcmc.pred$mvSamples)
proc.time() - ptm

plot(biomass,colMeans(samples.pred))


# # b) Here's how to do the Dirichlet-multinomial:

# # set up the "d" function for the distribution
# ddirchmulti <- nimbleFunction(
# run = function(x = double(1), alpha = double(1), size = double(0),
# log_value = integer(0)) {
# returnType(double(0))
# logProb <- lgamma(sum(alpha)) - sum(lgamma(alpha)) +
# sum(lgamma(alpha + x)) - lgamma(sum(alpha) + size)
# if(log_value) {
# return(logProb)
# } else {
# return(exp(logProb))
# }
# })

# # set up the "r" function
# rdirchmulti <- nimbleFunction(
# run = function(n = integer(0), alpha = double(1), size = double(0)) {
# returnType(double(1))
# if(n != 1) nimPrint("rdirchmulti only allows n = 1; using n = 1.")
# p <- rdirch(1, alpha)
# return(rmulti(1, size = size, prob = p))
# })

# # tell NIMBLE about the newly available distribution
# registerDistributions(list(
# ddirchmulti = list(
# BUGSdist = "ddirchmulti(alpha, size)"
# )
# ))

# # now modify your BUGS code and redo the stuff from part (a) above
# # e.g.,  y ~ ddirchmulti(alpha, size)
# # instead of y ~ dmulti(p, size)
# # p ~ ddirch(alpha)


