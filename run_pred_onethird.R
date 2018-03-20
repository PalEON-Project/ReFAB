library(nimble)
library(splines)
library(maps)
library(methods)

ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}

load("threethirds_v1.0.Rdata") 
#load("cast.x.Rdata")
load("sites_rm.Rdata")

Niters <- 10000
bMax <- 150
group_rm <- c('TWOTHIRDS_150')

#### Setting up 2/3 calibration 3/3 prediction
Y.keep <- Y
biomass.keep <- biomass
Y.calib <- Y[-sites_rm,]; Y.pred <- Y
biomass.calib <- biomass[-sites_rm]; biomass.pred <- biomass

#### Making sure Z.knots and u are the same between calibration and validation
u <- c(0,30,bMax) #c(rep(attr(Z.knots,"Boundary.knots")[1],1),attr(Z.knots,"knots"),rep(attr(Z.knots,"Boundary.knots")[2],1))

source("Workflow_Code/utils/bs_nimble.R")
Z.test <- matrix(NA,length(biomass.calib),5)
for(i in 1:length(biomass.calib)){
  Z.test[i,] <- bs_nimble(u_given = biomass.calib[i], u = u, N0 = rep(0, (length(u)-1)), N1 = rep(0, (length(u))),
                          N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))
}

Z.knots <- Z.test
source(file.path('Workflow_Code','calibration.model.R'))
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
                 bMax = bMax, Y = Y)


source('validation.R')
samples.pred <- validation_model(Y = Y.pred, Z.knots = Z.knots, 
                                 samples.mixed = samples.mixed, u = u,
                                 Niters = Niters, bMax = bMax, group_rm = group_rm,
                                 outLik = outLik)


pdf(paste0('gold.r2.validation.pdf'))
par(mfrow=c(1,1))
plot(biomass, colMeans(samples.pred, na.rm = T),
     xlim=c(0,bMax), ylim=c(0,bMax), pch=19,
     xlab="True Biomass", ylab="Predicted Mean Biomass")
abline(a=0,b=1)
lm.mod <- lm(biomass~colMeans(samples.pred)+0)
abline(lm.mod,lty=2)

lm.mod.out <- lm(biomass[sites_rm]~colMeans(samples.pred[,sites_rm])+0)
abline(lm.mod.out,lty=2,col='red')

points(biomass[sites_rm],colMeans(samples.pred[,sites_rm], na.rm = T),
       col='red',pch=19)
mtext(paste("r2-all",summary(lm.mod)$r.squared,'r2-out',summary(lm.mod.out)$r.squared))

arrows(x0 = biomass, y0 = apply(samples.pred,2,FUN = quantile,.05),
       x1 = biomass, y1 = apply(samples.pred,2,FUN = quantile,.975),
       code = 0, lwd=2)
arrows(x0 = biomass[sites_rm], y0 = apply(samples.pred[,sites_rm],2,FUN = quantile,.05),
       x1 = biomass[sites_rm], y1 = apply(samples.pred[,sites_rm],2,FUN = quantile,.975),
       code = 0, lwd=2, col = 'red')

dev.off()


