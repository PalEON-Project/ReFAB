#####
##### This code runs the validation model for 10 groups with beta uncertainty
##### using the 'job.array.sh' on geo
#####

library(nimble)
library(splines)
library(maps)
library(methods)

arg <- commandArgs(trailingOnly = TRUE)
if (is.na(arg[1])) {
  runnum <- NA
} else {
  runnum <- as.numeric(arg[1])
}
group_rm <- sort(rep(1:10,20))[runnum]

#### only need to load betas from the left out groups don't need to estimate betas 200 times
load(file = paste0("beta.est.group.in", group_rm, ".Rdata")) #'ALL_150' or group_rm

#### way to pick the the same betas across groups because the estimates are correlated so can't be random
dat.index <- data.frame(group_rm=sort(rep(1:10,20)),
                        beta_row =rep(round(seq(nrow(samples.mixed)*.2,nrow(samples.mixed),length.out = 20)),10), #picking betas past burnin
                        counter = rep(1:20,10))

beta_row <- dat.index[runnum, 'beta_row']

## loading twothirds calibration dataset
load("twothirds_v2.0.Rdata")
source(file.path('Workflow_Code','utils','validation_args.R')) #file with constants that should be constant between validation exercises

#### Setting up 10 fold cross validation
Y.keep <- Y
biomass.keep <- biomass

Y.calib <- Y[-sets10[,group_rm],]; Y.pred <- Y[sets10[,group_rm],]
biomass.calib <- biomass[-sets10[,group_rm]]; biomass.pred <- biomass[sets10[,group_rm]]

source("Workflow_Code/utils/bs_nimble.R")
Z.test <- matrix(NA,length(biomass.calib),5)
for(i in 1:length(biomass.calib)){
  Z.test[i,] <- bs_nimble(u_given = biomass.calib[i], u = u, N0 = rep(0, (length(u)-1)), N1 = rep(0, (length(u))),
                          N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))
}

Z.knots <- Z.test

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

### Calculate Likelihoods
source(file.path('Workflow_Code','utils','getLik.R'))
outLik <- getLik(Z = Z.new, u = u, beta = colMeans(samples.mixed),
                 bMax = bMax, Y = Y.pred,knots=length(u)+2)

save(outLik,file=paste0('outLik_group',group_rm,'_beta_',beta_row,'.Rdata'))

### Fit validation model
source(file.path('Workflow_Code','models','validation.R'))
samples.pred <- validation_model(Y = Y.pred, Z.knots = Z.knots, 
                                 samples.mixed = samples.mixed, u = u,
                                 Niters = Niters, bMax = bMax, group_rm = group_rm,
                                 outLik = outLik, beta_row = beta_row)
