
arg <- commandArgs(trailingOnly = TRUE)
if (is.na(arg[1])) {
  runnum <- NA
} else {
  runnum <- as.numeric(arg[1])
}
load(file = paste0("beta.est.group.in", group_rm, ".Rdata")) #'ALL_150' or group_rm

dat.index <- data.frame(group_rm=sort(rep(1:10,20)),
                        beta_row =rep(round(seq(nrow(samples.mixed)*.2,nrow(samples.mixed),length.out = 20)),10),
                        counter = rep(1:20,10))

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

#load("twothirds_v1.0.Rdata")
load('2018-11-16twothirds.calibration.data.Rdata')

Niters <- 10000
bMax <- 228

#### Setting up 10 fold cross validation
set.seed(5)
sets10 <- matrix(sample(x = 1:100,size = 100, replace = F),10,10)
Y.keep <- Y
biomass.keep <- biomass

if(is.na(group_rm)|group_rm > 10){
  Y.calib <- Y; Y.pred <- Y
  biomass.calib <- biomass; biomass.pred <- biomass
}else{
  Y.calib <- Y[-sets10[,group_rm],]; Y.pred <- Y[sets10[,group_rm],]
  biomass.calib <- biomass[-sets10[,group_rm]]; biomass.pred <- biomass[sets10[,group_rm]]
}

#### Making sure Z.knots and u are the same between calibration and validation
#Z.knots = bs(biomass.calib, intercept=TRUE, knots = 30, Boundary.knots=c(0,bMax))
u <- c(1,43,bMax) #c(rep(attr(Z.knots,"Boundary.knots")[1],1),attr(Z.knots,"knots"),rep(attr(Z.knots,"Boundary.knots")[2],1))

source("Workflow_Code/utils/bs_nimble.R")
Z.test <- matrix(NA,length(biomass.calib),5)
for(i in 1:length(biomass.calib)){
  Z.test[i,] <- bs_nimble(u_given = biomass.calib[i], u = u, N0 = rep(0, (length(u)-1)), N1 = rep(0, (length(u))),
                          N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))
}

Z.knots <- Z.test

#### only need to load betas from the left out groups don't need to estimate betas 200 times
load(file = paste0("beta.est.group.in", group_rm, ".Rdata")) #'ALL_150' or group_rm

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
                 bMax = bMax, Y = Y.pred,knots=length(u)+2)

save(outLik,file=paste0('outLik_group',group_rm,'_beta_',beta_row,'.Rdata'))

source(file.path('Workflow_Code','models','validation.R'))
samples.pred <- validation_model(Y = Y.pred, Z.knots = Z.knots, 
                                 samples.mixed = samples.mixed, u = u,
                                 Niters = Niters, bMax = bMax, group_rm = group_rm,
                                 outLik = outLik, beta_row = beta_row)

####
#### Plotting ####
####

if(FALSE){
dir_to_samples_pred <- c('~/Downloads/archive 4/')

samples.pred.mat <- array(NA,dim=c(5000,max(sets10),20))
for(i in 1:200){
  load(file = file.path(dir_to_samples_pred,
                        paste0('samples.pred.group',dat.index[i,'group_rm'],'beta',dat.index[i,'beta_row'],'.Rdata')))
  #load(file = file.path(paste0('samples.pred.group',i,'.Rdata')))
  if(any(is.na(samples.pred))) print(i)
  samples.pred.mat[,sets10[,dat.index[i,'group_rm']],dat.index[i,'counter']] <- samples.pred[,grep('b',colnames(samples.pred))]
}

pdf(paste0('ALL.calib.r2.validation.beta_uncert',Sys.Date(),'.pdf'))
par(mfrow=c(1,1))
plot(biomass.keep, apply(samples.pred.mat,2,FUN = quantile,.5,na.rm=T),
     xlim=c(0,bMax), ylim=c(0,bMax), pch=19,
     xlab="True Biomass", ylab="Predicted Mean Biomass",main='10 Fold CV')
abline(a=0,b=1)
lm.mod <- lm(biomass.keep ~ apply(samples.pred.mat,2,FUN = quantile,.5,na.rm=T)+0)
abline(lm.mod,lty=2)
mtext(paste("r-squared",summary(lm.mod)$r.squared))

arrows(x0 = biomass.keep, y0 = apply(samples.pred.mat,2,FUN = quantile,.05,na.rm=T),
       x1 = biomass.keep, y1 = apply(samples.pred.mat,2,FUN = quantile,.975,na.rm=T),
       code = 0, lwd=2)

library(calibrate)
textxy(biomass.keep,  apply(samples.pred.mat,2,FUN = quantile,.5),1:100)
dev.off()

pdf('trace.hists.beta.uncert.10fold.pdf')
par(mfrow=c(4,4))
for(i in 1:100){
  hist(samples.pred.mat[,i,],col='gray',freq = F,main=i,xlim=c(0,150))
  plot(samples.pred.mat[,i,1],typ='l',ylim=c(0,150),main = i)
  for(n in 2:20){
    points(samples.pred.mat[,i,n],typ='l')
  }
}
dev.off()

}
