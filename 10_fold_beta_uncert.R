
arg <- commandArgs(trailingOnly = TRUE)
if (is.na(arg[1])) {
  runnum <- NA
} else {
  runnum <- as.numeric(arg[1])
}

dat.index <- data.frame(group_rm=sort(rep(1:10,20)),beta_row =rep(round(seq(2000,10000,length.out = 20)),10))

group_rm <- dat.index[runnum, 'group_rm']
beta_row <- dat.index[runnum, 'beta_row']

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

Niters <- 10000
bMax <- 150

#### Setting up 10 fold cross validation
set.seed(5)
sets10 <- matrix(sample(x = 1:100,size = 100, replace = F),10,10)
Y.keep <- Y
biomass.keep <- biomass
Y.calib <- Y[-sets10[,group_rm],]; Y.pred <- Y[sets10[,group_rm],]
biomass.calib <- biomass[-sets10[,group_rm]]; biomass.pred <- biomass[sets10[,group_rm]]

if(is.na(group_rm)){
  Y.calib <- Y; Y.pred <- Y
  biomass.calib <- biomass; biomass.pred <- biomass
}

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
outLik <- getLik(Z = Z.new, u = u, beta = colMeans(samples.mixed),
                 bMax = bMax, Y = Y.pred)

source('validation.R')
samples.pred <- validation_model(Y = Y.pred, Z.knots = Z.knots, 
                                 samples.mixed = samples.mixed, u = u,
                                 Niters = Niters, bMax = bMax, group_rm = group_rm,
                                 outLik = outLik, beta_row)


