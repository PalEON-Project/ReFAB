library(nimble)
load("min.list.samples.Rdata")

code <- nimbleCode({
  
  #delta ~ dunif(0,1000)
  
  for(r in 1:R){ 
    for(i in 1:I){
      beta[r,i] ~ dnorm(0,.04)
    }
  }
  
  phi.first <- Z%*%beta
  
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
    p[j,] ~ ddirch(exp.phi[j,]+.025) #http://sourceforge.net/p/mcmc-jags/discussion/610037/thread/c21ef62a
    Y[j,] ~ dmulti(p[j,],n[j])
  }
  
})

z.seq1 = seq(1,400, length=5)
Z.knots = bs(biomass,intercept=TRUE,df=4)

beta = matrix(NA,ncol(Z.knots),ncol(Y))
p = matrix(NA,nrow(Y),ncol(Y)) ; phi = p
J = nrow(Y)

data = list(Y = counts , n = rowSums(counts), Z =  Z.knots)
constants = list(R = ncol(Z.knots), I = ncol(Y), J = nrow(Y))
inits.cal = list(list(beta = matrix(0,ncol(Z.knots),ncol(Y))),list(beta = matrix(3,ncol(Z.knots),ncol(Y))))
dimensions = list(exp.phi = dim(phi), Z = dim(Z.knots), beta = dim(beta), p = dim(p), Y = dim(counts))

model <- nimbleModel(code, inits = inits.cal, constants = constants, data = data.real.cal, dimensions = dimensions)

# # compiled version of the model
# Cmodel <- compileNimble(model)
# 
# # set up MCMC
# spec <- configureMCMC(model, thin = 10)
# spec$addMonitors(c('psi', 'theta'))  # set up monitoring of whatever
# model variables you want posterior samples for - by default, top level
# parameters are already included, so 'mu' in the above example would by
# default be monitored. 'psi' and 'theta' are just for illustration -
#   obviously they are not part of my toy model above
# 
# # create MCMC algorithm for the model
# Rmcmc <- buildMCMC(spec)
# # compiled version of the MCMC
# Cmcmc <- compileNimble(Rmcmc, project = model)
# 
# # run MCMC for 5000 iterations
# Cmcmc$run(5000)
# samples <- as.matrix(Cmcmc$mvSamples)
# 
# # samples has rows as iterations and columns as variables, you'll have
# to manually remove a burn-in period
# 
# b) Here's how to do the Dirichlet-multinomial:
# 
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
# 
# # set up the "r" function
# rdirchmulti <- nimbleFunction(
# run = function(n = integer(0), alpha = double(1), size = double(0)) {
# returnType(double(1))
# if(n != 1) nimPrint("rdirchmulti only allows n = 1; using n = 1.")
# p <- rdirch(1, alpha)
# return(rmulti(1, size = size, prob = p))
# })
# 
# # tell NIMBLE about the newly available distribution
# registerDistributions(list(
# ddirchmulti = list(
# BUGSdist = "ddirchmulti(alpha, size)"
# )
# ))
# 
# # now modify your BUGS code and redo the stuff from part (a) above
# # e.g.,  y ~ ddirchmulti(alpha, size)
# # instead of y ~ dmulti(p, size)
# # p ~ ddirch(alpha)