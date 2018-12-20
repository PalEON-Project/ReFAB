


ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}


load(file = paste0("beta.est.group.in", group_rm, ".Rdata"))

burnin <- round(.2 * nrow(samples.mixed))
new.biomass <- 1:bMax
Z.new = matrix(0,nrow=length(new.biomass),ncol=5)
for(i in 1:length(new.biomass)){
  u_given <- new.biomass[i]
  Z.new[i,] = bs_nimble(u_given, u=u, N0 = rep(0, (length(u)-1)),
                        N1 = rep(0, (length(u))), 
                        N2 = rep(0, (length(u)+1)), 
                        N3 = rep(0, (length(u)+2)))
}
source(file.path('Workflow_Code','utils','bs_nimble.R'))
outLik <- getLik(Z = Z.new, u = u, beta = (samples.mixed[nrow(samples.mixed),]),
                 bMax = bMax, Y = Y,knots=length(u)+2)

save(outLik, file=paste0('outLik.group.',group_rm,'.Rdata'))

if(FALSE){
  pdf('beta.hists.reorder.pdf')
  par(mfrow=c(4,4))
  for(i in sample(x = 1:ncol(samples.mixed),size = 4)){
    plot(samples.mixed[,i],typ='l')
    hist(samples.mixed[,i],col='gray',freq=F,main=colnames(samples.mixed)[i])
    lines(density(rnorm(nrow(samples.mixed), 0, sd = 5)), lwd =2)
  }
  dev.off()
  
  #u <- u #should have defined knots in calibration
  #u<-c(rep(attr(Z.knots,"Boundary.knots")[1],1),attr(Z.knots,"knots"),rep(attr(Z.knots,"Boundary.knots")[2],1))
  
  pdf(paste0('liks_linexp',group_rm,'.pdf'))
  par(mfrow=c(2,4))
  for(i in 1:nrow(Y)){#1:nrow(Y),
    plot(1:bMax, outLik[i,],main=i,typ='l')
  }
  dev.off()
  par(mfrow=c(1,1))
  plot(biomass,apply(outLik,1,which.max))
  abline(b=1,a=0)
  calibrate::textxy(biomass,apply(outLik,1,which.max),1:length(biomass))
  points(biomass[bimodal_sites],apply(outLik,1,which.max)[bimodal_sites],
         col='red',
         lwd=2)
}

source(file.path('Workflow_Code','models','validation.R'))
samples.pred <- validation_model(Y = Y.pred, Z.knots = Z.knots, 
                                 samples.mixed = samples.mixed, u = u,
                                 Niters = Niters, bMax = bMax, group_rm = group_rm,
                                 outLik = outLik)

#outlier <- which.max(abs(colMeans(samples.pred[,grep('b',colnames(samples.pred))]) - biomass.pred))

if(group_rm == 12){
  
  source('Workflow_Code/older_code/calibration.figs.R')
  calibration.figs(bMax = bMax, Z.knots = Z.knots, Y = Y.keep,
                   samples.mixed = samples.mixed, outLik = outLik,
                   biomass = biomass.keep, samples.pred = samples.pred,
                   group_rm = group_rm,Y.pred = Y.pred,
                   biomass.pred = biomass.pred, outlier = outlier,
                   sets10 = sets10)
  
  samples.pred.mat <- matrix(NA,nrow(samples.pred),length(biomass.keep))
  for(i in 1:10){
    load(file = file.path('~/Downloads','samps.linexp',paste0('samples.pred.group',i,'.Rdata')))
    #load(file = file.path(paste0('samples.pred.group',i,'.Rdata')))
    samples.pred.mat[,sets10[,i]] <- samples.pred[,grep('b',colnames(samples.pred))]
  }
  pdf(paste0('10.fold.r2.validation.linexp.pdf'))
  par(mfrow=c(1,1))
  plot(biomass.keep, colMeans(samples.pred.mat, na.rm = T),
       xlim=c(0,bMax), ylim=c(0,bMax), pch=19,
       xlab="True Biomass", ylab="Predicted Mean Biomass",main='10 Fold CV')
  abline(a=0,b=1)
  lm.mod <- lm(biomass.keep~ apply(samples.pred.mat,2,FUN = quantile,.5)+0)
  abline(lm.mod,lty=2)
  mtext(paste("r-squared",summary(lm.mod)$r.squared))
  
  arrows(x0 = biomass.keep, y0 = apply(samples.pred.mat,2,FUN = quantile,.05),
         x1 = biomass.keep, y1 = apply(samples.pred.mat,2,FUN = quantile,.975),
         code = 0, lwd=2)
  
  library(calibrate)
  textxy(biomass.keep,  apply(samples.pred.mat,2,FUN = quantile,.5),1:100)
  dev.off()
  #10, 1
  #9,10
  #7,5
  
  #calibrated on 2/3 predicted on all
  load('~/Downloads/samps.linexp/samples.pred.groupTWOTHIRDS_150.Rdata')
  load("threethirds_v1.0.Rdata") 
  #load("cast.x.Rdata")
  load("sites_rm.Rdata")
  
  pdf('2_3rds_validation.pdf')
  plot(biomass,colMeans(samples.pred),
       xlim=c(0,bMax), ylim=c(0,bMax), pch=19,
       xlab="True Biomass", ylab="Predicted Mean Biomass",
       main='2/3 calibration')
  lm.mod <- lm(biomass.keep~ apply(samples.pred,2,FUN = quantile,.5)+0)
  abline(lm.mod,lty=2)
  
  mtext(paste("r-squared",summary(lm.mod)$r.squared))
  # points(biomass,colMeans(samples.pred),
  #       pch=19, col='red')
  calibrate::textxy(biomass,colMeans(samples.pred),1:length(biomass))
  abline(b=1,a=0)
  dev.off()
  
}
####
#### Plotting ####
####

if(FALSE){
  dir_to_samples_pred <- c('~/Downloads/preds/')
  
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
  plot(biomass.keep[1:100], apply(samples.pred.mat,2,FUN = quantile,.5,na.rm=T),
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

##threethirds

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

###3/3 fit to 3/3

pdf(paste0('twothirds-calib.r2.validation.pdf'))
par(mfrow=c(1,1))
plot(biomass, colMeans(samples.pred, na.rm = T),
     xlim=c(0,bMax), ylim=c(0,bMax), pch=19,
     xlab="True Biomass", ylab="Predicted Mean Biomass")
abline(a=0,b=1)
lm.mod <- lm(biomass~colMeans(samples.pred)+0)
abline(lm.mod,lty=2)

#points(biomass[sites_rm],colMeans(samples.pred[,sites_rm], na.rm = T),
#       col='red',pch=19)
#points(biomass[bimodal_sites],colMeans(samples.pred)[bimodal_sites],
#      col='red',
#     lwd=2)
mtext(paste("r2-twothirds",summary(lm.mod)$r.squared))

arrows(x0 = biomass, y0 = apply(samples.pred,2,FUN = quantile,.05),
       x1 = biomass, y1 = apply(samples.pred,2,FUN = quantile,.975),
       code = 0, lwd=2)
library(calibrate)
textxy(biomass,colMeans(samples.pred),1:length(biomass))

dev.off()

