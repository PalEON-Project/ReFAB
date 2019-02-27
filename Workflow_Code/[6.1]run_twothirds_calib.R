
#####
##### This code runs the calibration model for 10 groups and the 2/3s fit
##### using the 'job.array.calib.valid.sh' on geo
#####


arg <- commandArgs(trailingOnly = TRUE)
if (is.na(arg[1])) {
  runnum <- NA
} else {
  group_rm <- as.numeric(arg[1])
  runnum <- group_rm
}

library(nimble)
library(splines)
library(maps)
library(methods)

load('twothirds_v2.0.Rdata')

source(file.path('Workflow_Code','utils','validation_args.R')) #file with constants that should be constant between validation exercises
#source(file.path('Workflow_Code','utils','7knot_args.R')) #file with constants that should be constant between validation exercises

#### Setting up 10 fold cross validation
Y.keep <- Y

#### For running with only arborel pollen
arboreal = TRUE
if(arboreal == TRUE){
  Y <- Y.keep[,-which(colnames(Y)%in%c('prairie','other_herbs','CYPERACE'))]
  Niters <- 10000
  group_rm <- paste0(runnum,'nograss')
}

#biomass.keep <- biomass

FULL = FALSE
#### Adds biomass data product uncertainty
if(FULL==TRUE){
  load('biomass_draws.Rdata')
  biomass.keep <- unlist(lapply(biomass_draws,function(x) x[runnum]))
  load("~/ReFAB/new_sites_rm.Rdata") # new 1/3
  load('~/ReFAB/sites_rm.Rdata') # old 1/3
  load("~/ReFAB/TF.Rdata") # old full
  
  bfull <- biomass.keep[TF]
  btt <- bfull[-sites_rm]
  b3 <- bfull[sites_rm]
  newbfull <- biomass.keep[-TF]
  newbtt <- newbfull[-sites_rm_new]
  newb3 <- newbfull[sites_rm_new]
  
  biomass1_3 <- c(b3,newb3)
  biomass2_3 <- c(btt,newbtt)
  biomass <- biomass2_3
}

if(is.na(group_rm) | group_rm > 10){ #if true does full 2/3 fit
  Y.calib <- Y; Y.pred <- Y
  biomass.calib <- biomass; biomass.pred <- biomass
}else{
  Y.calib <- Y[-sets10[,group_rm],]; Y.pred <- Y[sets10[,group_rm],]
  biomass.calib <- biomass[-sets10[,group_rm]]; biomass.pred <- biomass[sets10[,group_rm]]
}

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
  



