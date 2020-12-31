
mse_func <- function(preds,actual) {
  mean((preds - actual)^2)
}

source(file.path('Workflow_Code','utils','give_me_R2.R'))
source(file.path('Workflow_Code','utils','validation_args.R'))

#source(file.path('Workflow_Code','utils','7knot_args.R'))
load('twothirds_v3.0.Rdata')
WANT.NUMS = FALSE

#####
##### Plot 10 Fold CV R2 Validation - 5 Knots
#####

dir_to_samples_pred <- c('~/Dropbox/ReFAB_outputs/samps_5knots_agb/')
samples.pred.mat <- array(NA,dim=c(1000,max(sets10)))
dat.index <- data.frame(group_rm = sort(rep(1:10,1)),
                        beta_row = 10000, #picking betas past burnin
                        counter = rep(1:10,1))
for(i in 1:10){
  load(file = file.path(dir_to_samples_pred,
                        paste0('samples.pred.group',dat.index[i,'group_rm'],'beta',dat.index[i,'beta_row'],'.Rdata')))
  #load(file = file.path(paste0('samples.pred.group',i,'.Rdata')))
  if(any(is.na(samples.pred))) print(i)
  samples.pred.mat[,sets10[,dat.index[i,'group_rm']]] <- samples.pred[,grep('b',colnames(samples.pred))]
}

pdf(paste0('[6-B]10.fold.R2.5knots',Sys.Date(),'.pdf'))
par(mfrow=c(1,1))
biomass.keep <- biomass[1:max(sets10)]
plot(biomass.keep, apply(samples.pred.mat,2,FUN = quantile,.5,na.rm=T),
     xlim=c(0,bMax), ylim=c(0,bMax), pch=19,
     xlab="True Biomass (Mg/ha)", ylab="Predicted Mean Biomass (Mg/ha)",
     main='Five Knots')
abline(a=0,b=1)
lm.mod <- lm(biomass.keep ~ apply(samples.pred.mat,2,FUN = quantile,.5,na.rm=T)+0)
abline(lm.mod,lty=2)
R2 <- give_me_R2( preds = apply(samples.pred.mat,2,FUN = quantile,.5,na.rm=T),
                  actual = biomass.keep)
MSE <- mse_func(preds = apply(samples.pred.mat,2,FUN = quantile,.5,na.rm=T),
                actual = biomass.keep)
legend('bottomright',c(paste("R2 = ",signif(R2,digits = 3)),
                       paste('MSE =', signif(MSE,digits = 3))))

arrows(x0 = biomass.keep, y0 = apply(samples.pred.mat,2,FUN = quantile,.05,na.rm=T),
       x1 = biomass.keep, y1 = apply(samples.pred.mat,2,FUN = quantile,.975,na.rm=T),
       code = 0, lwd=2)

if(WANT.NUMS ==TRUE) calibrate::textxy(biomass.keep,  apply(samples.pred.mat,2,FUN = quantile,.5),1:max(sets10))
dev.off()

#####
##### Plot 10 Fold CV R2 Validation - NO GRASS
#####

dir_to_samples_pred <- c('~/Downloads/samps_nograss/')
Niters <- 10000
samples.pred.mat <- array(NA,dim=c(1000,max(sets10)))
dat.index <- data.frame(group_rm = sort(rep(1:10,1)),
                        beta_row = 10000, #picking betas past burnin
                        counter = rep(1:10,1))
for(i in 1:10){
  load(file = file.path(dir_to_samples_pred,
                        paste0('samples.pred.group',dat.index[i,'group_rm'],'nograssbeta',dat.index[i,'beta_row'],'.Rdata')))
  #load(file = file.path(paste0('samples.pred.group',i,'.Rdata')))
  if(any(is.na(samples.pred))) print(i)
  samples.pred.mat[,sets10[,dat.index[i,'group_rm']]] <- samples.pred[,grep('b',colnames(samples.pred))]
}

pdf(paste0('[6-A]10.fold.R2.nograss',Sys.Date(),'.pdf'))
par(mfrow=c(1,1))
biomass.keep <- biomass[1:max(sets10)]
plot(biomass.keep, apply(samples.pred.mat,2,FUN = quantile,.5,na.rm=T),
     xlim=c(0,bMax), ylim=c(0,bMax), pch=19,
     xlab="True Biomass (Mg/ha)", ylab="Predicted Mean Biomass (Mg/ha)",
     main='Arboreal Only')
abline(a=0,b=1)
lm.mod <- lm(biomass.keep ~ apply(samples.pred.mat,2,FUN = quantile,.5,na.rm=T)+0)
abline(lm.mod,lty=2)
R2 <- give_me_R2( preds = apply(samples.pred.mat,2,FUN = quantile,.5,na.rm=T),
                  actual = biomass.keep)
MSE <- mse_func(preds = apply(samples.pred.mat,2,FUN = quantile,.5,na.rm=T),
                actual = biomass.keep)
legend('bottomright',c(paste("R2 = ",signif(R2,digits = 3)),
                       paste('MSE =', signif(MSE,digits = 3))))

arrows(x0 = biomass.keep, y0 = apply(samples.pred.mat,2,FUN = quantile,.05,na.rm=T),
       x1 = biomass.keep, y1 = apply(samples.pred.mat,2,FUN = quantile,.975,na.rm=T),
       code = 0, lwd=2)

if(WANT.NUMS ==TRUE) calibrate::textxy(biomass.keep,  apply(samples.pred.mat,2,FUN = quantile,.5),1:max(sets10))
dev.off()


#####
##### Plot 10 Fold CV R2 Validation
#####

dir_to_samples_pred <- c('~/Dropbox/ReFAB_outputs/samps_10cv/')#c('~/10samps/')
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

pdf(paste0('[5]new_10.fold.R2',Sys.Date(),'.pdf'))
par(mfrow=c(1,1))
biomass.keep <- biomass[1:max(sets10)]
plot(biomass.keep, apply(samples.pred.mat,2,FUN = quantile,.5,na.rm=T),
     xlim=c(0,bMax), ylim=c(0,bMax), pch=19,
     xlab="True Biomass (Mg/ha)", ylab="Predicted Mean Biomass (Mg/ha)")
abline(a=0,b=1)
lm.mod <- lm(biomass.keep ~ apply(samples.pred.mat,2,FUN = quantile,.5,na.rm=T)+0)
abline(lm.mod,lty=2)
#mtext(paste("r-squared = ",signif(summary(lm.mod)$r.squared,digits = 3)))

R2 <- give_me_R2(preds = apply(samples.pred.mat,2,FUN = quantile,.5,na.rm=T),
                  actual =  biomass.keep)
MSE <- mse_func(preds = apply(samples.pred.mat,2,FUN = quantile,.5,na.rm=T),
                actual =  biomass.keep)
legend('bottomright',c(paste("R2 = ",signif(R2,digits = 3)),
                       paste('MSE =', signif(MSE,digits = 3))))


arrows(x0 = biomass.keep, y0 = apply(samples.pred.mat,2,FUN = quantile,.05,na.rm=T),
       x1 = biomass.keep, y1 = apply(samples.pred.mat,2,FUN = quantile,.975,na.rm=T),
       code = 0, lwd=2)

if(WANT.NUMS ==TRUE) calibrate::textxy(biomass.keep,  apply(samples.pred.mat,2,FUN = quantile,.5),1:max(sets10))
dev.off()

##### Blue Splines Plot for 10 FOLD CV
source(file.path('Workflow_Code','utils','splines_plot.R'))
samples.mixed.all <- rep(0,308)#rep(0,220)#
for(i in 1:10){
  #load(paste0('~/Downloads/betas_10cv/3beta.est.group.in',i,'.Rdata'))
  load(paste0('~/Downloads/5knot_betas/5beta.est.group.in',i,'.Rdata'))
  ppp <- nrow(samples.mixed)
  samples.mixed.all <- rbind(samples.mixed.all,
                           samples.mixed[seq(ppp*.2,ppp,length.out = 250),]) 
  print(i)
}

source('~/ReFAB/Workflow_Code/utils/validation_args.R', echo=TRUE)
source('~/ReFAB/Workflow_Code/utils/7knot_args.R', echo=TRUE)

load('twothirds_v3.0.Rdata')

pdf('5knot_splines.10CV.pdf',height = 5,width = 8)
splines_plot(samples.mixed = samples.mixed.all,Y = Y,biomass = biomass,
             bMax = bMax,order_plot = 1:22)#c(1:3,10,12,18)
dev.off()


#####
##### 2/3s fit validation #####
#####

load('~/Downloads/samps_FULL/samples.pred.group11betaNA.Rdata')
load("threethirds_v3.0.Rdata") 
load("~/ReFAB/new_sites_rm.Rdata") # new 1/3
load('~/ReFAB/sites_rm.Rdata') # old 1/3
load('~/ReFAB/TF.Rdata') # old full

origfull <- samples.pred[,TF]
origtt <- origfull[,-sites_rm]
orig3 <- origfull[,sites_rm]
newfull <- samples.pred[,-TF]
newtt <- newfull[,-sites_rm_new]
new3 <- newfull[,sites_rm_new]

samps1_3 <- cbind(orig3,new3)
samps2_3 <- cbind(origtt,newtt)
#save(samps2_3,file='samps2_3.Rdata')

bfull <- biomass[TF]
btt <- bfull[-sites_rm]
b3 <- bfull[sites_rm]
newbfull <- biomass[-TF]
newbtt <- newbfull[-sites_rm_new]
newb3 <- newbfull[sites_rm_new]

biomass1_3 <- c(b3,newb3)
biomass2_3 <- c(btt,newbtt)

save(biomass1_3,biomass2_3,samps1_3,samps2_3,file='split_calib_dat_v3.0.Rdata')

load('~/Dropbox/ReFAB_outputs/split_calib_dat_v3.0.Rdata')

pdf('new_training_set_2_3rds_validation_agb.pdf')
par(mfrow=c(1,1))
plot(biomass2_3,colMeans(samps2_3),
     xlim=c(0,bMax), ylim=c(0,bMax), pch=19,
     xlab="True Biomass", ylab="Predicted Mean Biomass",
     main='Training Set')
lm.mod <- lm(apply(samps2_3,2,FUN = quantile,.5)~biomass2_3+0)
abline(lm.mod,lty=2)

#mtext(paste("r-squared",signif(summary(lm.mod)$r.squared,digits=3)))

R2 <- give_me_R2(preds = apply(samps2_3,2,FUN = quantile,.5),
                 actual =  biomass2_3)
mtext(paste("r-squared = ",signif(R2,digits = 3)))

arrows(x0 = biomass2_3, y0 = apply(samps2_3,2,FUN = quantile,.05,na.rm=T),
       x1 = biomass2_3, y1 = apply(samps2_3,2,FUN = quantile,.975,na.rm=T),
       code = 0, lwd=2)

if(WANT.NUMS ==TRUE) calibrate::textxy(biomass2_3,colMeans(samps2_3),1:length(biomass2_3))
abline(b=1,a=0)
dev.off()

#####
##### 2/3s to 1/3s fit R2
#####

pdf(paste0(Sys.Date(),'[8]new_gold.r2.validation.pdf'))
par(mfrow=c(1,1))
plot(biomass2_3, apply(samps2_3,2,FUN = quantile,.5),
     xlim=c(0,bMax), ylim=c(0,bMax), pch=19,
     xlab="True Biomass (Mg/ha)", ylab="Predicted Mean Biomass (Mg/ha)")
abline(a=0,b=1)
lm.mod <- lm(apply(samps2_3,2,FUN = quantile,.5)~biomass2_3+0)
abline(lm.mod,lty=2)

lm.mod.out <- lm(apply(samps1_3,2,FUN = quantile,.5)~biomass1_3+0)
abline(lm.mod.out,lty=2,col='red')

points(biomass1_3,colMeans(samps1_3, na.rm = T),
       col='red',pch=19)
#mtext(paste("r2-twothirds = ",signif(summary(lm.mod)$r.squared,digits=2),'r2-onethird = ',signif(summary(lm.mod.out)$r.squared,digits = 2)))

#R2training = signif(summary(lm.mod)$r.squared,digits=3)
#R2testing = signif(summary(lm.mod.out)$r.squared,digits = 3)

R2training <- signif(give_me_R2(preds = apply(samps2_3,2,FUN = quantile,.5),
                 actual =  biomass2_3),digits = 3)
R2testing <- signif(give_me_R2(preds = apply(samps1_3,2,FUN = quantile,.5),
                         actual =  biomass1_3),digits = 3)

mse_train <- signif(mse_func(preds = apply(samps2_3,2,FUN = quantile,.5),
                               actual =  biomass2_3),digits = 3)
mse_test <- signif(mse_func(preds = apply(samps1_3,2,FUN = quantile,.5),
                             actual =  biomass1_3),digits = 3)

arrows(x0 = biomass2_3, y0 = apply(samps2_3,2,FUN = quantile,.05),
       x1 = biomass2_3, y1 = apply(samps2_3,2,FUN = quantile,.975),
       code = 0, lwd=2)
arrows(x0 = biomass1_3, y0 = apply(samps1_3,2,FUN = quantile,.05),
       x1 = biomass1_3, y1 = apply(samps1_3,2,FUN = quantile,.975),
       code = 0, lwd=2, col = 'red')
# legend(
#   'bottomright',
#   c(paste0('Training (R2 = ', R2training,', MSE = ',mse_train,')'),
#     paste0('Testing (R2 = ', R2testing,', MSE = ',mse_test,')')),
#   pch=19,col=c('black','red'))

dev.off()

R2training <- signif(give_me_R2(preds = apply(samps2_3,2,FUN = quantile,.5)[biomass2_3<150],
                                actual =  biomass2_3[biomass2_3<150]),digits = 3)
R2testing <- signif(give_me_R2(preds = apply(samps1_3,2,FUN = quantile,.5)[biomass1_3<150],
                               actual =  biomass1_3[biomass1_3<150]),digits = 3)
print(paste(R2training,R2testing))

save(samps1_3,biomass1_3,samps2_3,biomass2_3,file='valid_out.Rdata')


#####
##### 3/3s fit R2
#####

load('biomass_draws_v3.0.Rdata')
biomass <- unlist(lapply(biomass_draws,function(x) mean(x[100:120])))
biomass.05 <- unlist(lapply(biomass_draws,function(x) quantile(x[100:120],.05)))
biomass.95 <- unlist(lapply(biomass_draws,function(x) quantile(x[100:120],.95)))


load('~/Downloads/samps_FULL/samples.pred.group100FULLbetaNA.Rdata')
samples.pred.keep <- samples.pred
samples.pred.ar <- array(NA,c(5000,232,150))
samples.pred.ar[,,100] <- samples.pred
#plot(colMeans(samples.pred),unlist(lapply(biomass_draws,function(x) (x[100]))))
for(i in 101:150){
  load(paste0('~/Downloads/samps_FULL/samples.pred.group',i,'FULLbetaNA.Rdata'))
  samples.pred.keep <- rbind(samples.pred.keep,samples.pred)
  samples.pred.ar[,,i] <- samples.pred
 # points(colMeans(samples.pred),unlist(lapply(biomass_draws,function(x) (x[i]))))
  }


pdf('[9]new_3_3rds_validation.pdf')
par(mfrow=c(1,1))
plot(biomass,colMeans(samples.pred.keep),
     xlim=c(0,bMax), ylim=c(0,bMax), pch=19,
     xlab="True Biomass (Mg/ha)", ylab="Predicted Mean Biomass (Mg/ha)")
lm.mod <- lm(biomass~ apply(samples.pred.keep,2,FUN = quantile,.5)+0)
abline(a=0,b=1)
abline(lm.mod,lty=2)

#mtext(paste("r-squared",signif(summary(lm.mod)$r.squared,digits=3)))

R2 <- give_me_R2(preds = apply(samples.pred.keep,2,FUN = quantile,.5),
                 actual =  biomass)
MSE <- mse_func(preds = apply(samples.pred.keep,2,FUN = quantile,.5,na.rm=T),
                actual =  biomass)
legend('bottomright',c(paste("R2 = ",signif(R2,digits = 3)),
                       paste('MSE =', signif(MSE,digits = 3))))

arrows(x0 = biomass, y0 = apply(samples.pred.keep,2,FUN = quantile,.05,na.rm=T),
       x1 = biomass, y1 = apply(samples.pred.keep,2,FUN = quantile,.975,na.rm=T),
       code = 0, lwd=2,col='gray')

arrows(x0 = biomass.05, y0 = colMeans(samples.pred.keep),
       x1 = biomass.95, y1 = colMeans(samples.pred.keep),
       code = 0, lwd=2,col='gray')
points(biomass,colMeans(samples.pred.keep),pch=19)


if(WANT.NUMS ==TRUE) calibrate::textxy(biomass,colMeans(samples.pred),1:length(biomass))
dev.off()

#########
######### Likelihood Plots for Unstable Sites #########
#########

load('threethirds_v3.0.Rdata')
dir_to_outs <- c('~/Downloads/outs_FULL/')
out.save <- array(NA,c(232,bMax,50))
for(i in 1:50){
  load(file = file.path(dir_to_outs,
                        paste0('outLik',seq(100,150,1)[i],'FULL.Rdata')))
  if(any(is.na(outLik))) print(i)
  out.save[,,i] <- outLik
}

pdf('[11]unstable_liks_calib_FULL.pdf',height=12,width = 14)
plot_sites <- c(9,13,16,27,28,32,35,36,59,69,70,71,93,102,107,129,187,210)
par(mfrow=c(4,5),oma=rep(0,4),mar=c(4,4,2,2))
for(s in plot_sites){
  out <- out.save[s,,]
  plot(exp(out[,1]-max(out[,1]))/-sum(out[,1]),typ='l',main=paste('site',s),
       xlab='Biomass',ylab='Likelihood', lwd=1.5)
  abline(v=biomass[s],col='blue',lwd=2,lty=2)
  
  rug_do <- apply(out,2,which.max)
  rug(rug_do,col='red',lwd=2)
  
}
#if(length(which(s==seq(15,232,by=16)))==1){
  plot.new()
  legend('center',c('Training Likelihood','Max. Liks. Ests. from \nother betas','Settlement'),col=c('black','red','blue'),lty = c(1,1,2),lwd=2)
#}
dev.off()



###################### Other stuff

dir_to_outs <- c('~/Downloads/outs_10cv')
out.save <- array(NA,c(bMax,150,20))
for(i in 1:200){
  load(file = file.path(dir_to_outs,
                        paste0('outLik_group',dat.index[i,'group_rm'],
                               '_beta_',dat.index[i,'beta_row'],'.Rdata')))
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

load('~/Downloads/betas_10cv/3beta.est.group.in11.Rdata')
source(file.path('Workflow_Code','utils','getLik.R'))
outLik <- getLik(Z = Z.new, u = u, beta = colMeans(samples.mixed[burnin:nrow(samples.mixed),]),
                 bMax = bMax, Y = Y, knots=length(u)+2)
save(outLik,file='outLiktwothirds_v3.0.Rdata')

load('twothirds_v3.0.Rdata')
plot(biomass,colMeans(samps2_3))

pdf('max.liks.23.calib_agb.pdf')
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



pdf('max.liks.33.calib.pdf')
par(mfrow=c(3,3))
for(s in 1:150){
  out <- out.save[,s,]
  plot(exp(outLik[s,]-max(outLik[s,]))/-sum(outLik[s,]),typ='l',main=paste('site',s),
       xlab='Biomass',ylab='Likelihood')
  #points(exp(outLik[s,]-max(outLik[s,]))/-sum(outLik[s,]),typ='l',col='red')
  abline(v=biomass3_3[s],col='blue')
  #abline(v=colMeans(samps2_3)[s],col='red')
  
  
  
  rug(samps3_3[,s],col='red')
  if(length(which(s==seq(8,150,by=8)))==1){
    plot.new()
    legend('center',c('Training likelihood','Max. liks. ests. from \n betas draws','settlement'),col=c('black','red','blue'),pch=1)
  }
}
dev.off()

#####
##### Calculating SD of Max Liks
#####

load('threethirds_v3.0.Rdata')

source(file.path('Workflow_Code','utils','getLik.R'))
out_calc <- list()
for(i in 1:50){
  load(paste0('~/Downloads/betas_FULL/3beta.est.group.in',seq(100,150,1)[i],'FULL.Rdata'))
  out_calc[[i]] <- getLik(Z = Z.new, u = u, beta = (samples.mixed[nrow(samples.mixed),]),
                          bMax = bMax, Y = Y, knots = length(u) + 2)
  
}

par(mfrow=c(3,3))
for(i in 1:20){
  plot(exp(out_calc[[i]][s,]-max(out_calc[[i]][s,]))/-sum(out_calc[[i]][s,]))
}

save(out_calc,file='out_calc_v3.0.Rdata')

out_max_lik <- out_max_lik2 <- list()
exp_out <- array(NA,c(20,232,228))
for(i in 2:20){
  for(s in 1:232){
    exp_out[i,s,] <- (exp(out_calc[[i]][s,]-max(out_calc[[i]][s,]))/-sum(out_calc[[i]][s,]))
  }
  out_max_lik[[i]] <- apply(out_calc[[i]],1,which.max)
  out_max_lik2[[i]] <- apply(exp_out[i,,],1,which.max)
}

out_mat <- do.call(rbind,out_max_lik)
out_mat2 <- do.call(rbind,out_max_lik2)

out_sd <- apply(out_mat,2,sd)
out_sd2 <- apply(out_mat2,2,sd)



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

