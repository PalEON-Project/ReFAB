#####
##### This code runs the validation model for 3/3s of the calibration data
##### using betas from the two thirds calibration fit
##### without beta uncertainty included 
##### no .sh just source
##### 

library(nimble)
library(splines)
library(maps)
library(methods)

## loading threethirds calibration dataset and constants
load("threethirds_v2.0.Rdata")
source(file.path('Workflow_Code','utils','validation_args.R')) #file with constants that should be constant between validation exercises

group_rm <- 11

#### Setting up 3/3 prediction
Y.keep <- Y
biomass.keep <- biomass
Y.pred <- Y.calib <- Y
biomass.pred <- biomass.calib <- biomass

source("Workflow_Code/utils/bs_nimble.R")
Z.test <- matrix(NA,length(biomass.calib),5)
for(i in 1:length(biomass.calib)){
  Z.test[i,] <- bs_nimble(u_given = biomass.calib[i], u = u, N0 = rep(0, (length(u)-1)), N1 = rep(0, (length(u))),
                          N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))
}
Z.knots <- Z.test

##loading 2/3s calibration fit
load(file = paste0("beta.est.group.in", group_rm, ".Rdata"))

burnin <- round(.2 * nrow(samples.mixed))
new.biomass <- 1:bMax
Z.new = matrix(0,nrow=length(new.biomass),ncol=ncol(Z.knots))
for(i in 1:length(new.biomass)){
  u_given <- new.biomass[i]
  Z.new[i,] = bs_nimble(u_given, u=u, N0 = rep(0, (length(u)-1)),
                        N1 = rep(0, (length(u))), 
                        N2 = rep(0, (length(u)+1)), 
                        N3 = rep(0, (length(u)+2)))
}

source(file.path('Workflow_Code','utils','getLik.R'))
outLik <- getLik(Z = Z.new, u = u, beta = colMeans(samples.mixed[burnin:nrow(samples.mixed),]),
                 bMax = bMax, Y = Y.pred)

source('validation.R')
samples.pred <- validation_model(Y = Y.pred, Z.knots = Z.knots, 
                                 samples.mixed = samples.mixed, u = u,
                                 Niters = Niters, bMax = bMax, group_rm = group_rm,
                                 outLik = outLik)
