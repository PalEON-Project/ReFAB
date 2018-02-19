#### Y                     dataset with # of sites rows and # of taxa columns
#### biomass               dataset with length = # of sites
#### Niters                number of iterations
#### Niters.save           number of iterations to output to Rdata file
#### group_rm              number 1-10 for which group to remove

calibration_model <- function(Y, biomass, Z.knots, u, Niters = 5000,
                              Niter.save = 200, group_rm = NA){
  
library(nimble)
source(file.path('genPareto','betabin.R')) # code for user-defined beta-binomial distribution
  
calib_code <- nimbleCode({
    
    for(r in 1:R){ 
      for(i in 1:I){
        beta1[r,i] ~ dnorm(0,.04)
        beta2[r,i] ~ dnorm(0,.04)
      }  
    }

    shape1 <- exp(Z[,] %*% beta1[,])
    shape2 <- exp(Z[,] %*% beta2[,])

    for(j in 1:J){
      Y[j, 1] ~ dbetabin(shape1[j, 1], shape2[j, 1], n[j])
      for(i in 2:(I-1)){
        Y[j, i] ~ dbetabin(shape1[j, i], shape2[j, i], n[j] - sum(Y[j,1:(i-1)]))
      }
    }

  })

#load("Data/calibration.data.Rdata")
model.dir <- c('Workflow Code/')
fig.dir <- c("Figures/")
source("Workflow_Code/utils/bs_nimble.R")

J = nrow(Y)

data = list(Y = as.matrix(Y), Z =  Z.knots)

constants = list(n = rowSums(Y), R = ncol(Z.knots), I = ncol(Y),
                 J = nrow(Y))

inits = list(beta1 = matrix(1, ncol(Z.knots), ncol(Y)),
             beta2 = matrix(1, ncol(Z.knots), ncol(Y)))

dimensions = list(shape1 = c(nrow(Y),ncol(Y)),
                  shape2 = c(nrow(Y),ncol(Y)), 
                  Z = dim(Z.knots), 
                  beta1 = c(ncol(Z.knots), ncol(Y)), 
                  beta2 = c(ncol(Z.knots), ncol(Y)),
                  Y = dim(Y), n = nrow(Y))

# in BUGS code, to calculate the vector of basis matrix values for a given biomass, 
#pass that biomass in as 'u_given', pass in the vector of u values for the knots and 
#pass in N0,N1,N2,N3 of correct length - you can do this simply by providing N0,N1,N2,N3 
#as part of the 'constants' argument given to the 'nimbleModel' function
model <- nimbleModel(calib_code, inits = inits, constants = constants,
                     data = data, dimensions = dimensions)

# compiled version of the model
Cmodel <- compileNimble(model)

# set up MCMC
spec <- configureMCMC(model, print = TRUE, useConjugacy = FALSE)
spec$addMonitors(c('beta1','beta2')) 

# set up monitoring of whatever
# model variables you want posterior samples for - by default, top level
# parameters are already included, so 'mu' in the above example would by
# default be monitored. 'psi' and 'theta' are just for illustration -
# obviously they are not part of my toy model above

# create MCMC algorithm for the model
Rmcmc <- buildMCMC(spec)

# compiled version of the MCMC
Cmcmc <- compileNimble(Rmcmc, project = model)

# run MCMC for 2000 iterations
set.seed(0)
Cmcmc$run(Niters)#50000
samples.mixed <- as.matrix(Cmcmc$mvSamples)
save(samples.mixed, Y, biomass, file = paste0("beta.est.group", group_rm, ".Rdata"))

return(samples.mixed)

}

