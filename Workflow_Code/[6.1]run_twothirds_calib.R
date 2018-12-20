
#####
##### This code runs the calibration model for 10 groups and the 2/3s fit
##### using the 'job.array.calib.valid.sh' on geo
#####


arg <- commandArgs(trailingOnly = TRUE)
if (is.na(arg[1])) {
  runnum <- NA
} else {
  group_rm <- as.numeric(arg[1])
}

library(nimble)
library(splines)
library(maps)
library(methods)

load('twothirds_v2.0.Rdata')

source(file.path('Workflow_Code','utils','validation_args.R')) #file with constants that should be constant between validation exercises

#### Setting up 10 fold cross validation
Y.keep <- Y
biomass.keep <- biomass

if(is.na(group_rm) | group_rm > 10){ #if true does full 2/3 fit
  Y.calib <- Y; Y.pred <- Y
  biomass.calib <- biomass; biomass.pred <- biomass
}else{
  Y.calib <- Y[-sets10[,group_rm],]; Y.pred <- Y[sets10[,group_rm],]
  biomass.calib <- biomass[-sets10[,group_rm]]; biomass.pred <- biomass[sets10[,group_rm]]
}

u <- c(0,median_use,bMax)

source("Workflow_Code/utils/bs_nimble.R")
Z.test <- matrix(NA,length(biomass.calib),length(u)+2)
for(i in 1:length(biomass.calib)){
  Z.test[i,] <- bs_nimble(u_given = biomass.calib[i], u = u, N0 = rep(0, (length(u)-1)), N1 = rep(0, (length(u))),
                          N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))
}

Z.knots <- Z.test

source(file.path('Workflow_Code','models','calibration.model.R'))
samples.mixed <- calibration_model(Y = Y.calib, biomass = biomass.calib,
                                     Z.knots = Z.knots, u = u, Niters = Niters,
                                     group_rm = group_rm)
  



