
arg <- commandArgs(trailingOnly = TRUE)
if (is.na(arg[1])) {
  runnum <- NA
} else {
  group_rm <- as.numeric(arg[1])
}

### after get.data

library(nimble)
library(splines)
library(maps)
library(methods)

ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}

# histogram of betas to see how it looks relative to prior
# the prior might be too tight now. with linexp
# rel to an sd of 5. 
# might want a flat prior like => 1/400 precision

load("twothirds_v1.0.Rdata")

Niters <- 50000
bMax <- 150

#### Setting up taxa removal 
Y.keep <- Y
biomass.keep <- biomass
Y.calib <- Y.pred <- Y[,-group_rm]
biomass.calib <- biomass.pred <- biomass

#### Making sure Z.knots and u are the same between calibration and validation
#Z.knots = bs(biomass.calib, intercept=TRUE, knots = 30, Boundary.knots=c(0,bMax))
u <- c(0,30,bMax) #c(rep(attr(Z.knots,"Boundary.knots")[1],1),attr(Z.knots,"knots"),rep(attr(Z.knots,"Boundary.knots")[2],1))

source("Workflow_Code/utils/bs_nimble.R")
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

load(file = paste0("beta.est.taxa.", group_rm, ".Rdata"))

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
outLik <- getLik(Z = Z.new, u = u, beta = (samples.mixed[nrow(samples.mixed),]),
                 bMax = bMax, Y = Y)

source(file.path('Workflow_Code','models','validation.R'))
samples.pred <- validation_model(Y = Y.pred, Z.knots = Z.knots, 
                                 samples.mixed = samples.mixed, u = u,
                                 Niters = Niters, bMax = bMax, group_rm = group_rm,
                                 outLik = outLik)
