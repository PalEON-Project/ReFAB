
arg <- commandArgs(trailingOnly = TRUE)
if (is.na(arg[1])) {
  runnum <- NA
} else {
  group_rm <- as.numeric(arg[1])
}

if(group_rm == 11){
  group_rm <- NA
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

load("2018-02-28all.calibration.data.Rdata") 
load("sites_rm.Rdata")

Y <- Y[-sites_rm,]
biomass <- biomass[-sites_rm]

Niters <- 5000
bMax <- 150

#### Setting up 10 fold cross validation
set.seed(5)
sets10 <- matrix(sample(x = 1:100,size = 100, replace = F),10,10)
Y.keep <- Y
biomass.keep <- biomass
Y.calib <- Y[-sets10[,group_rm],]; Y.pred <- Y[sets10[,group_rm],]
biomass.calib <- biomass[-sets10[,group_rm]]; biomass.pred <- biomass[sets10[,group_rm]]

#### Making sure Z.knots and u are the same between calibration and validation
Z.knots = bs(biomass.calib, intercept=TRUE, knots = 30, Boundary.knots=c(0,bMax))
u <- c(0,30,bMax) #c(rep(attr(Z.knots,"Boundary.knots")[1],1),attr(Z.knots,"knots"),rep(attr(Z.knots,"Boundary.knots")[2],1))

source(file.path('Workflow_Code','calibration.model.R'))
samples.mixed <- calibration_model(Y = Y.calib, biomass = biomass.calib,
                                   Z.knots = Z.knots, u = u, Niters = Niters,
                                   group_rm = group_rm)

source('validation.R')
samples.pred <- validation_model(Y = Y.pred, Z.knots = Z.knots, 
                 samples.mixed = samples.mixed, u = u,
                 Niters = Niters, bMax = bMax, group_rm = group_rm)

load(file=paste0('outLik.group.',group_rm,'.Rdata'))
outlier <- which.max(abs(colMeans(samples.pred[,grep('b',colnames(samples.pred))]) - biomass.pred))

source('calibration.figs.R')
calibration.figs(bMax = bMax, Z.knots = Z.knots, Y = Y.keep,
                 samples.mixed = samples.mixed, outLik = outLik,
                 biomass = biomass.keep, samples.pred = samples.pred,
                 group_rm = group_rm,Y.pred = Y.pred,
                 biomass.pred = biomass.pred, outlier = outlier,
                 sets10 = sets10)

if(group_rm == 10){
  samples.pred.mat <- matrix(NA,nrow(samples.pred),length(biomass.keep))
  for(i in 1:10){
    load(file = file.path('~/Downloads','samps.22',paste0('samples.pred.group',i,'.Rdata')))
    #load(file = file.path(paste0('samples.pred.group',i,'.Rdata')))
    samples.pred.mat[,sets10[,i]] <- samples.pred[,grep('b',colnames(samples.pred))]
  }
  pdf(paste0('10.fold.r2.validation.pdf'))
  par(mfrow=c(1,1))
  plot(biomass, colMeans(samples.pred.mat, na.rm = T),
       xlim=c(0,bMax), ylim=c(0,bMax), pch=19,
       xlab="True Biomass", ylab="Predicted Mean Biomass")
  abline(a=0,b=1)
  lm.mod <- lm(biomass~colMeans(samples.pred.mat)+0)
  abline(lm.mod,lty=2)
  mtext(paste("r-squared",summary(lm.mod)$r.squared))
  
  arrows(x0 = biomass, y0 = apply(samples.pred.mat,2,FUN = quantile,.05),
         x1 = biomass, y1 = apply(samples.pred.mat,2,FUN = quantile,.975),
         code = 0, lwd=2)
  
  library(calibrate)
  textxy(biomass, colMeans(samples.pred.mat, na.rm = T),1:100)
  dev.off()
  #10, 1
  #9,10
  #7,5
  
}



