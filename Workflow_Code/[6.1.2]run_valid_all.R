#####
##### This code runs the calibration and validation for 3/3s dataset
##### creates final betas. 
#####

library(nimble)
library(splines)
library(maps)
library(methods)

load("threethirds_v2.0.Rdata")
source(file.path('Workflow_Code','utils','validation_args.R')) #file with constants that should be constant between validation exercises and full runs
group_rm <- c('FULL')

#### Setting up 3/3 calibration 3/3 prediction
Y.keep <- Y
biomass.keep <- biomass
Y.calib <- Y; Y.pred <- Y
biomass.calib <- biomass; biomass.pred <- biomass

source(file.path("Workflow_Code","utils","bs_nimble.R"))
Z.test <- matrix(NA,length(biomass.calib),5)
for(i in 1:length(biomass.calib)){
  Z.test[i,] <- bs_nimble(u_given = biomass.calib[i], u = u, N0 = rep(0, (length(u)-1)), N1 = rep(0, (length(u))),
                          N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))
}

Z.knots <- Z.test

source(file.path('Workflow_Code','models','calibration.model.R'))
samples.mixed <- calibration_model(Y = Y.calib, biomass = biomass.calib,
                                     Z.knots = Z.knots, u = u, Niters = Niters,
                                     group_rm = group_rm)

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

source(file.path('Workflow_Code','models','validation.R'))
samples.pred <- validation_model(Y = Y.pred, Z.knots = Z.knots, 
                                 samples.mixed = samples.mixed, u = u,
                                 Niters = Niters, bMax = bMax, group_rm = group_rm,
                                 outLik = outLik)

