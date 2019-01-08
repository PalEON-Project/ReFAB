
source(file.path('Workflow_Code','utils','validation_args.R'))
load('twothirds_v2.0.Rdata')
WANT.NUMS = FALSE

#####
##### Plot 10 Fold CV R2 Validation
#####

dir_to_samples_pred <- c('~/10samps/')
samples.pred.mat <- array(NA,dim=c(5000,max(sets10),20))
dat.index <- data.frame(group_rm=sort(rep(1:10,20)),
                        beta_row =rep(round(seq(Niters*.2,Niters,length.out = 20)),10), #picking betas past burnin
                        counter = rep(1:20,10))
for(i in 1:200){
  load(file = file.path(dir_to_samples_pred,
                        paste0('samples.pred.group',dat.index[i,'group_rm'],'beta',dat.index[i,'beta_row'],'.Rdata')))
  #load(file = file.path(paste0('samples.pred.group',i,'.Rdata')))
  if(any(is.na(samples.pred))) print(i)
  samples.pred.mat[,sets10[,dat.index[i,'group_rm']],dat.index[i,'counter']] <- samples.pred[,grep('b',colnames(samples.pred))]
}

pdf(paste0('10.fold.R2',Sys.Date(),'.pdf'))
par(mfrow=c(1,1))
biomass.keep <- biomass[1:max(sets10)]
plot(biomass.keep, apply(samples.pred.mat,2,FUN = quantile,.5,na.rm=T),
     xlim=c(0,bMax), ylim=c(0,bMax), pch=19,
     xlab="True Biomass (Mg/ha)", ylab="Predicted Mean Biomass (Mg/ha)",main='10 Fold Cross Validation')
abline(a=0,b=1)
lm.mod <- lm(biomass.keep ~ apply(samples.pred.mat,2,FUN = quantile,.5,na.rm=T)+0)
abline(lm.mod,lty=2)
mtext(paste("r-squared = ",signif(summary(lm.mod)$r.squared,digits = 2)))

arrows(x0 = biomass.keep, y0 = apply(samples.pred.mat,2,FUN = quantile,.05,na.rm=T),
       x1 = biomass.keep, y1 = apply(samples.pred.mat,2,FUN = quantile,.975,na.rm=T),
       code = 0, lwd=2)

if(WANT.NUMS ==TRUE) calibrate::textxy(biomass.keep,  apply(samples.pred.mat,2,FUN = quantile,.5),1:max(sets10))
dev.off()

##### Blue Splines Plot for 10 FOLD CV
source(file.path('Workflow_Code','utils','splines_plot.R'))
samples.mixed.all <- rep(0,220)
for(i in 1:10){
  load(paste0('~/Downloads/betas/beta.est.group.in',i,'.Rdata'))
  ppp <- nrow(samples.mixed)
  samples.mixed.all <- rbind(samples.mixed.all,
                           samples.mixed[seq(ppp*.2,ppp,length.out = 250),]) 
  print(i)
}

pdf('splines.10CV.pdf')
splines_plot(samples.mixed = samples.mixed.all,Y = Y,biomass = biomass,
             bMax = bMax)
dev.off()


#####
##### 2/3s fit validation #####
#####

load('~/Downloads/samples.pred.group11betaNA.Rdata')
load("threethirds_v2.0.Rdata") 
load("~/ReFAB/new_sites_rm.Rdata") # new 1/3
load('~/ReFAB/sites_rm.Rdata') # old 1/3
load('~/ReFAB/TF.Rata') # old full

origfull <- samples.pred[,TF]
origtt <- origfull[,-sites_rm]
orig3 <- origfull[,sites_rm]
newfull <- samples.pred[,-TF]
newtt <- newfull[,-sites_rm_new]
new3 <- newfull[,sites_rm_new]

samps1_3 <- cbind(orig3,new3)
samps2_3 <- cbind(origtt,newtt)

bfull <- biomass[TF]
btt <- bfull[-sites_rm]
b3 <- bfull[sites_rm]
newbfull <- biomass[-TF]
newbtt <- newbfull[-sites_rm_new]
newb3 <- newbfull[sites_rm_new]

biomass1_3 <- c(b3,newb3)
biomass2_3 <- c(btt,newbtt)


pdf('2_3rds_validation_withnums.pdf')
par(mfrow=c(1,1))
plot(biomass2_3,colMeans(samps2_3),
     xlim=c(0,bMax), ylim=c(0,bMax), pch=19,
     xlab="True Biomass", ylab="Predicted Mean Biomass",
     main='2/3 validation')
lm.mod <- lm(biomass2_3~ apply(samps2_3,2,FUN = quantile,.5)+0)
abline(lm.mod,lty=2)

mtext(paste("r-squared",summary(lm.mod)$r.squared))

arrows(x0 = biomass2_3, y0 = apply(samps2_3,2,FUN = quantile,.05,na.rm=T),
       x1 = biomass2_3, y1 = apply(samps2_3,2,FUN = quantile,.975,na.rm=T),
       code = 0, lwd=2)

if(WANT.NUMS ==TRUE) calibrate::textxy(biomass2_3,colMeans(samps2_3),1:length(biomass2_3))
abline(b=1,a=0)
dev.off()

#### Likelihood Plots for 10 CV versus 2/3s
dir_to_outs <- c('~/10outs/')
out.save <- array(NA,c(228,150,20))
for(i in 1:200){
  load(file = file.path(dir_to_outs,
                        paste0('outLik_group',dat.index[i,'group_rm'],'_beta_',dat.index[i,'beta_row'],'.Rdata')))
  print(paste0('outLik_group',dat.index[i,'group_rm'],'_beta_',dat.index[i,'beta_row'],'.Rdata'))
  #load(file = file.path(paste0('samples.pred.group',i,'.Rdata')))
  if(any(is.na(outLik))) print(i)
  out.save[,sets10[,dat.index[i,'group_rm']],dat.index[i,'counter']] <- t(outLik)
  outLik <-NULL
}

burnin <- round(.2 * nrow(samples.mixed))
new.biomass <- 1:bMax
source("Workflow_Code/utils/bs_nimble.R")
biomass.calib <- biomass[-sites_rm]

Z.test <- matrix(NA,length(biomass.calib),length(u)+2)
for(i in 1:length(biomass.calib)){
  Z.test[i,] <- bs_nimble(u_given = biomass.calib[i], u = u, N0 = rep(0, (length(u)-1)), N1 = rep(0, (length(u))),
                          N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))
}
Z.knots <- Z.test
Z.new = matrix(0,nrow=length(new.biomass),ncol=ncol(Z.knots))
for(i in 1:length(new.biomass)){
  u_given <- new.biomass[i]
  Z.new[i,] = bs_nimble(u_given, u=u, N0 = rep(0, (length(u)-1)),
                        N1 = rep(0, (length(u))), 
                        N2 = rep(0, (length(u)+1)), 
                        N3 = rep(0, (length(u)+2)))
}

load('~/Downloads/betas/beta.est.group.in11.Rdata')
source(file.path('Workflow_Code','utils','getLik.R'))
outLik <- getLik(Z = Z.new, u = u, beta = colMeans(samples.mixed[burnin:nrow(samples.mixed),]),
                 bMax = bMax, Y = Y, knots=length(u)+2)
save(outLik,'outLiktwothirds.Rdata')

load('twothirds_v2.0.Rdata')
plot(biomass,colMeans(samps2_3))

pdf('max.liks.23.calib.pdf')
par(mfrow=c(3,3))
for(s in 1:150){
  out <- out.save[,s,]
  plot(exp(out[,1]-max(out[,1]))/-sum(out[,1]),typ='l',main=paste('site',s),
       xlab='Biomass',ylab='Likelihood')
  points(exp(outLik[s,]-max(outLik[s,]))/-sum(outLik[s,]),typ='l',col='red')
  abline(v=biomass2_3[s],col='blue')
  abline(v=colMeans(samps2_3)[s],col='red')
  rug(samps2_3[,s],col='red')
  if(length(which(s==seq(8,150,by=8)))==1){
    plot.new()
    legend('center',c('10 fold CV','2/3s','rug is MCMC','settlement'),col=c('black','red','red','blue'),pch=1)
  }
}
dev.off()



#####
##### 2/3s to 1/3s fit R2
#####

pdf(paste0('gold.r2.validation.pdf'))
par(mfrow=c(1,1))
plot(biomass2_3, colMeans(samps2_3, na.rm = T),
     xlim=c(0,bMax), ylim=c(0,bMax), pch=19,
     xlab="True Biomass", ylab="Predicted Mean Biomass",
     main = 'Gold Validation')
abline(a=0,b=1)
lm.mod <- lm(biomass2_3~colMeans(samps2_3)+0)
abline(lm.mod,lty=2)

lm.mod.out <- lm(biomass2_3~colMeans(samps2_3)+0)
abline(lm.mod.out,lty=2,col='red')

points(biomass1_3,colMeans(samps1_3, na.rm = T),
       col='red',pch=19)
mtext(paste("r2-twothirds = ",signif(summary(lm.mod)$r.squared,digits=2),'r2-onethird = ',signif(summary(lm.mod.out)$r.squared,digits = 2)))

arrows(x0 = biomass2_3, y0 = apply(samps2_3,2,FUN = quantile,.05),
       x1 = biomass2_3, y1 = apply(samps2_3,2,FUN = quantile,.975),
       code = 0, lwd=2)
arrows(x0 = biomass1_3, y0 = apply(samps1_3,2,FUN = quantile,.05),
       x1 = biomass1_3, y1 = apply(samps1_3,2,FUN = quantile,.975),
       code = 0, lwd=2, col = 'red')

dev.off()

#####
##### 3/3s fit R2
#####

load('~/Downloads/FULL.preds/samples.pred.group100FULLbetaNA.Rdata')
samples.pred.keep <- samples.pred
samples.pred.ar <- array(NA,c(5000,232,120))
samples.pred.ar[,,100] <- samples.pred
plot(colMeans(samples.pred),unlist(lapply(biomass_draws,function(x) (x[100]))))
for(i in 101:120){
  load(paste0('~/Downloads/FULL.preds/samples.pred.group',i,'FULLbetaNA.Rdata'))
  samples.pred.keep <- rbind(samples.pred.keep,samples.pred)
  samples.pred.ar[,,i] <- samples.pred
  points(colMeans(samples.pred),unlist(lapply(biomass_draws,function(x) (x[i]))))
  }

load('biomass_draws.Rdata')
biomass <- unlist(lapply(biomass_draws,function(x) mean(x[100:120])))
biomass.05 <- unlist(lapply(biomass_draws,function(x) quantile(x[100:120],.05)))
biomass.95 <- unlist(lapply(biomass_draws,function(x) quantile(x[100:120],.95)))

pdf('3_3rds_validation.pdf')
par(mfrow=c(1,1))
plot(biomass,colMeans(samples.pred.keep),
     xlim=c(0,bMax), ylim=c(0,bMax), pch=19,
     xlab="True Biomass", ylab="Predicted Mean Biomass",
     main='3/3 validation')
lm.mod <- lm(biomass~ apply(samples.pred.keep,2,FUN = quantile,.5)+0)
abline(a=0,b=1)
abline(lm.mod,lty=2)

mtext(paste("r-squared",summary(lm.mod)$r.squared))

arrows(x0 = biomass, y0 = apply(samples.pred.keep,2,FUN = quantile,.05,na.rm=T),
       x1 = biomass, y1 = apply(samples.pred.keep,2,FUN = quantile,.975,na.rm=T),
       code = 0, lwd=2)

arrows(x0 = biomass.05, y0 = colMeans(samples.pred.keep),
       x1 = biomass.95, y1 = colMeans(samples.pred.keep),
       code = 0, lwd=2)


if(WANT.NUMS ==TRUE) calibrate::textxy(biomass,colMeans(samples.pred),1:length(biomass))
dev.off()


pdf('max.liks.33.calib.pdf')
par(mfrow=c(3,3))
for(s in 1:150){
  out <- out.save[,s,]
  plot(exp(out[,1]-max(out[,1]))/-sum(out[,1]),typ='l',main=paste('site',s),
       xlab='Biomass',ylab='Likelihood')
  points(exp(outLik[s,]-max(outLik[s,]))/-sum(outLik[s,]),typ='l',col='red')
  abline(v=biomass2_3[s],col='blue')
  abline(v=colMeans(samps2_3)[s],col='red')
  rug(samps2_3[,s],col='red')
  if(length(which(s==seq(8,150,by=8)))==1){
    plot.new()
    legend('center',c('10 fold CV','2/3s','rug is MCMC','settlement'),col=c('black','red','red','blue'),pch=1)
  }
}
dev.off()

#####
####################################### Notes below
#####

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

