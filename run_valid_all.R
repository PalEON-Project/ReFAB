library(nimble)
library(splines)
library(maps)
library(methods)

ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}

load("2018-02-28all.calibration.data.Rdata") 
#load("cast.x.Rdata")
#load("sites_rm.Rdata")

Niters <- 10000
bMax <- 143
group_rm <- c('ALL_143')

#### Setting up 2/3 calibration 3/3 prediction
Y.keep <- Y
biomass.keep <- biomass
Y.calib <- Y; Y.pred <- Y
biomass.calib <- biomass; biomass.pred <- biomass

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

pdf(paste0('all4calib.r2.validation.pdf'))
par(mfrow=c(1,1))
plot(biomass, colMeans(samples.pred, na.rm = T),
     xlim=c(0,bMax), ylim=c(0,bMax), pch=19,
     xlab="True Biomass", ylab="Predicted Mean Biomass")
abline(a=0,b=1)
lm.mod <- lm(biomass~colMeans(samples.pred)+0)
abline(lm.mod,lty=2)

points(biomass[sites_rm],colMeans(samples.pred[,sites_rm], na.rm = T),
       col='red',pch=19)
mtext(paste("r2-all",summary(lm.mod)$r.squared))

arrows(x0 = biomass, y0 = apply(samples.pred,2,FUN = quantile,.05),
       x1 = biomass, y1 = apply(samples.pred,2,FUN = quantile,.975),
       code = 0, lwd=2)
dev.off()


