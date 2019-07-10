

######## run on geo with job.array.taxa.sh

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

ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}

load("twothirds_v2.0.Rdata") #for calibration
load('threethirds_v2.0.Rdata') #for validation

source(file.path('Workflow_Code','utils','validation_args.R')) #file with constants that should be constant between validation exercises

#### Setting up taxa removal 
Y.keep <- Y
biomass.keep <- biomass
Y.calib <- Y.pred <- Y[,-group_rm]
biomass.calib <- biomass.pred <- biomass

source("Workflow_Code/utils/bs_nimble.R")
Z.test <- matrix(NA,length(biomass.calib),5)
for(i in 1:length(biomass.calib)){
  Z.test[i,] <- bs_nimble(u_given = biomass.calib[i], u = u, N0 = rep(0, (length(u)-1)), N1 = rep(0, (length(u))),
                          N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))
}

Z.knots <- Z.test

source(file.path('Workflow_Code','models','calibration.model.R'))
if(FALSE){
  samples.mixed <- calibration_model(Y = Y.calib, biomass = biomass.calib,
                                     Z.knots = Z.knots, u = u, Niters = Niters,
                                     group_rm = group_rm)
}

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
outLik <- getLik(Z = Z.new, u = u, beta = (samples.mixed[nrow(samples.mixed),]),
                 bMax = bMax, Y = Y.calib, knots = length(u) + 2)

source(file.path('Workflow_Code','models','validation.R'))
samples.pred <- validation_model(Y = Y.pred, Z.knots = Z.knots, 
                                 samples.mixed = samples.mixed, u = u,
                                 Niters = Niters, bMax = bMax, group_rm = group_rm,
                                 outLik = outLik)

if(FALSE){ ### To plot use the following after running on the server
  samps.mat <- array(NA, dim = c(5000,232,22))
  r.saved <- r.saved.23 <- numeric(22)
  
  load('threethirds_v2.0.Rdata') #for validation
  
  load('one.third.idx.Rdata')
  
  for(i in 1:22){
    load(paste0('~/Downloads/iter.samps/samples.pred.group',i,'betaNA.Rdata'))
    samps.mat[,,i] <- samples.pred
  }
  
  pdf(paste0('iterative.taxa.sens',Sys.Date(),'.pdf'))
  par(mfrow=c(2,2))
  for(i in order(r.saved)){ #order(r.saved)
    
    bio.median <- apply(samps.mat[,,i],2,FUN = quantile,.5)
    plot(biomass,bio.median,
         xlim=c(0,bMax), ylim=c(0,bMax), pch=19,
         xlab="True Biomass", ylab="Predicted Mean Biomass",
         main=colnames(Y)[i],col='black',cex = .5)
    points(biomass[one.third.idx], bio.median[one.third.idx],pch=19,col='red',cex = .5)
    abline(a=0,b=1)
    
    lm.mod <- lm(bio.median[-one.third.idx] ~ biomass[-one.third.idx] +0)
    lm.mod.13 <- lm(bio.median[one.third.idx] ~ biomass[one.third.idx] +0)
    
    abline(lm.mod,lty=2)
    abline(lm.mod.13,lty=2,col='red')
    
    r.saved.23[i] <- give_me_R2(preds = bio.median[-one.third.idx],
                        actual = biomass[-one.third.idx])#signif(summary(lm.mod)$r.squared,digits = 3)
    r.saved[i] <- give_me_R2(preds = bio.median[one.third.idx],
                        actual = biomass[one.third.idx])#signif(summary(lm.mod.13)$r.squared,digits = 3)
    
    r2_23 <- signif(r.saved.23[i],digits = 3)
    r2_13 <- signif(r.saved[i],digits = 3)
    
    arrows(x0 = biomass, y0 = apply(samps.mat[,,i],2,FUN = quantile,.05),
           x1 = biomass, y1 = apply(samps.mat[,,i],2,FUN = quantile,.975),
           code = 0, lwd=.5)
    arrows(x0 = biomass[one.third.idx], y0 = apply(samps.mat[,,i],2,FUN = quantile,.05)[one.third.idx],
           x1 = biomass[one.third.idx], y1 = apply(samps.mat[,,i],2,FUN = quantile,.975)[one.third.idx],
           code = 0, lwd=.5,col='red')
    
    legend('bottomright',c(paste('R2 =',r2_23),paste('R2 =',r2_13)),text.col = c('black','red'))
    
  }
  dev.off()  
  
  
  out.table <- data.frame(TwoThirds=r.saved.23,OneThird=r.saved)
  rownames(out.table) <- colnames(Y)
  out.table <- out.table[order(out.table[,2]),]
  xtable(out.table)
  
}




