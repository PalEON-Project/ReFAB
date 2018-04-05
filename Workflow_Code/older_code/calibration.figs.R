
calibration.figs <- function(bMax, Z.knots, Y, samples.mixed, outLik,
                             biomass, samples.pred,group_rm, Y.pred,
                             biomass.pred, outlier, sets10){
  
  ciEnvelope <- function(x,ylo,yhi,...){
    polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                        ylo[1])), border = NA,...) 
  }
  
  blue       <- col2rgb("blue")
  alphablue  <- rgb(blue[1], blue[2], blue[3], 75, max = 255)

  orange       <- col2rgb("orange")
  alphaorange  <- rgb(orange[1], orange[2], orange[3], 75, max = 255)
    
  source("Workflow_Code/utils/bs_nimble.R")
  
  i.beta1 <- grep("beta1",colnames(samples.mixed))
  i.beta2 <- grep("beta2",colnames(samples.mixed))
  Z <- Z.knots
  
  ####
  #### Beta2Shape Trace Plots ####
  ####
  
  new.biomass <- 1:bMax
  Z.new = matrix(0,nrow=length(new.biomass),ncol=ncol(Z))
  u <- u #should have defined knots in calibration
  u<-c(rep(attr(Z.knots,"Boundary.knots")[1],1),attr(Z.knots,"knots"),rep(attr(Z.knots,"Boundary.knots")[2],1))
  
  for(i in 1:length(new.biomass)){
    u_given <- new.biomass[i]
    Z.new[i,] = bs_nimble(u_given, u=u, N0 = rep(0, (length(u)-1)),N1 = rep(0, (length(u))), N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))
  }
  
  alpha <- beta <- array(NA,dim=c(length(new.biomass),ncol(Y),
                                  nrow(samples.mixed)))
  
  for(i in 1:nrow(samples.mixed)){
    beta1 = matrix(samples.mixed[i,i.beta1],ncol(Z),ncol(Y))
    beta2 = matrix(samples.mixed[i,i.beta2],ncol(Z),ncol(Y))
    
    alpha[,,i] <- (Z.new%*%beta1)
    beta[,,i] <- (Z.new%*%beta2)
  }
  
  pdf(paste0('shape.trace.lim.taxa',group_rm,'.pdf'))
  par(mfrow=c(4,4))
  for(b in seq(1,length(new.biomass),length.out = 5)){
    for(s in 1:ncol(Y)){
      ## Add label for numbering
      plot(alpha[b,s,], typ = 'l',main=paste('alpha',b,s))
      hist(alpha[b,s,],col='gray',freq=F,main=paste('alpha',b,s))
      plot(beta[b,s,], typ = 'l',main=paste('beta',b,s))
      hist(beta[b,s,],col='gray',freq=F,main=paste('alpha',b,s))
    }
  }
  dev.off()
  
  ####
  #### Validataion Max Likelihood Plots ####
  ####
  vals <- 1:bMax
  J = nrow(Y.pred)
  
  pdf(paste0('validation_max_liks_pfts_',group_rm,'.pdf'))
  par(mfrow=c(4,5))
  for(j in 1:J){
    for(s in 1:(ncol(Y)-1)){
      if(sum(outLik[,j,s])!=0){
        plot(vals,exp(outLik[,j,s] - max(outLik[,j,s]))/-sum(outLik[,j,s])
             ,typ='l',ylab=NA,main=paste('site',sets10[,group_rm][j],'taxa',colnames(Y)[s]),xlab = 'biomass')
        abline(v=biomass.pred[j],col='blue',lwd=2)
        abline(v=colMeans(samples.pred)[j],col='orange',lwd=2)
      }else{
        plot.new()
      }
    }
  }
  
  for(j in 1:J){
    plot(apply(outLik[,j,],1,sum)[1:(bMax-1)],typ='l',ylab=NA,
         main=paste('site',sets10[,group_rm][j]),
         xlab = 'biomass')
    abline(v=biomass.pred[j],col='blue',lwd=2)
    abline(v=colMeans(samples.pred)[j],col='orange',lwd=2)
  }
  
  dev.off()
  
  ####
  #### Validataion R2 Plot ####
  ####
  
  pdf(paste0('validation.r2.group',group_rm,'.pdf'))
  par(mfrow=c(1,1))
  plot(biomass.pred,
       colMeans(samples.pred[,grep('b',colnames(samples.pred))]),
       xlim=c(0,bMax),ylim=c(0,bMax),pch=19,xlab="True Biomass",
       ylab="Predicted Mean Biomass")
  abline(a=0,b=1)
  abline(lm(biomass.pred~colMeans(samples.pred[,grep('b',colnames(samples.pred))])+0),lty=2)
  mtext(paste("r-squared",summary(lm(biomass.pred ~ colMeans(samples.pred[,grep('b',colnames(samples.pred))])+0))$r.squared))
  
  arrows(x0 = biomass.pred, y0 = apply(samples.pred,2,FUN = quantile,.05),
         x1 = biomass.pred, y1 = apply(samples.pred,2,FUN = quantile,.975),
         code = 0, lwd=2)
  dev.off()
  
  ####
  #### Splines ####
  ####
  
  if(FALSE){
    K <- 143 # n biomass values
    M <- 200 # n posterior samples
    I # number of taxa
    p_rel <- array(0, c(K, I, M))
    
    for(m in 1:M) {
      devs <- runif(I)
      devs <- rep(devs, each = K)
      ## use same random number for all K biomass values
      p_rel[1:K, 1:I, m] <- qbeta(devs, shape1[1:K , 1:I, m], shape2[1:K , 1:I, m])
    }
  }
  
  #i.beta1 <- i.beta[-i.beta.pine]
  #iter.look <- 250 #4474 low15 #4368 high14
  
  #beta1 = matrix(samples.mixed[iter.look,i.beta1],ncol(Z),ncol(Y))
  #beta2 = matrix(samples.mixed[iter.look,i.beta2],ncol(Z),ncol(Y))
  
  beta1 = matrix(colMeans(samples.mixed[,i.beta1]),ncol(Z),ncol(Y))
  beta2 = matrix(colMeans(samples.mixed[,i.beta2]),ncol(Z),ncol(Y))
  
  alpha <- exp(Z.new%*%beta1)
  beta <- exp(Z.new%*%beta2)
  
  mean.betas.ex <- alpha/(alpha+beta)
  sd.betas.ex <- sqrt(alpha*beta/((alpha+beta)^2*(alpha+beta+1)))
  
  sd.betas.ex[,2] <- sd.betas.ex[,2] * (1 - (sd.betas.ex[,1]))
  mean.betas.ex[,2] <- mean.betas.ex[,2] * (1 - (mean.betas.ex[,1]))
  for(i in 3:(ncol(Y)-1)){
    sd.betas.ex[,i] <- sd.betas.ex[,i] * (1 - rowSums(sd.betas.ex[,1:i-1]))
    mean.betas.ex[,i] <- mean.betas.ex[,i] * (1 - rowSums(mean.betas.ex[,1:i-1]))
  }
  sd.betas.ex[,ncol(Y)] <- 1 - rowSums(sd.betas.ex[,1:(ncol(Y)-1)])
  mean.betas.ex[,ncol(Y)] <- 1 - rowSums(mean.betas.ex[,1:(ncol(Y)-1)])
  
  props <- prop.table(as.matrix(Y),1)
  props.pred <- prop.table(as.matrix(Y.pred),1)
  
  #outlier <- as.numeric(names(sort(abs(biomass-colMeans(samples.pred[,grep('b',colnames(samples.pred))])))))[90]
  
  pdf(paste0('outlier.splines.group.',group_rm,'.pdf'))
  pdf('splines.together.pdf')
  par(mfrow=c(2,2))
  for(i in 1:ncol(Y)){
    plot(1:bMax, mean.betas[,i],
         ylim=c(0,max(mean.betas[,i] + sd.betas[,i]*2,
                      mean.betas[,i] - sd.betas[,i]*2, 
                      props[,i])),
         main = colnames(Y)[i],pch=19)
    points(1:bMax, mean.betas.ex[,i],col='red')
    ciEnvelope(x = 1:bMax, yhi = mean.betas[,i] + sd.betas[,i]*2,
               ylo = mean.betas[,i] - sd.betas[,i]*2,col = alphablue)
    
    ciEnvelope(x = 1:bMax, yhi = mean.betas.ex[,i] + sd.betas.ex[,i]*2,
               ylo = mean.betas.ex[,i] - sd.betas.ex[,i]*2,col = alphaorange)
    points(biomass,props[,i])
    #points(biomass.pred,props.pred[,i],col='red',pch=19)
    legend('topright',c('bMax = 150','bMax = 143'),pch=c(19,19),col=c(alphablue,alphaorange))
    ###Validation outlier
    if(FALSE){
    points(biomass.pred[outlier],
           props.pred[outlier,i],
           col='darkblue',pch=19,cex=2)
    points(colMeans(samples.pred[,grep('b',colnames(samples.pred))])[outlier],
           props.pred[outlier,i],
           col='orange',pch=19,cex=2)
    }
    if(FALSE){
    points(biomass[c(outlier)],
           props[c(outlier),i],
           col='blue',pch=19,cex=2.5)
    points(mean(samples.pred[,c(outlier)]),
           props[c(outlier),i],col='red',pch=19,cex=2.5)
    }
    ###top ten validation outliers
    if(FALSE){
      points(biomass[c(94,75,70,69,100)],
             props[c(94,75,70,69,100),i],
             col='blue',pch=19)
      points(colMeans(samples.pred[,c(94,75,70,69,100)]),
             props[c(94,75,70,69,100),i],col='red',pch=19)
      #points(qually.biomass,qually.props[,i],col='purple',pch=19,cex=2)
    }
  }
  dev.off()
  
  
  ####
  #### Old Spline Plots ####
  ####
  
if(FALSE){  
prop.quants <- matrix(NA,ncol(samples.mixed),3)
for(i in 1:ncol(samples.mixed)){
  prop.quants[i,]<-quantile(samples.mixed[,i],c(.025,.5,.975),na.rm=TRUE)
}
rownames(prop.quants)<-colnames(samples.mixed)

plot.help<- seq(0,3290,145)
total_counts <- rowSums(counts)

if(DRAW == TRUE) pdf(file.path(fig.dir,paste0(Sys.Date(),'tele.betas.calib_horiz_plus.pdf')))
par(mfrow=c(2,2))
for(i in 1:22){
  plot(1:145,prop.quants[grep('p.true1',
                              colnames(samples.mixed))[(plot.help[i]+1):plot.help[i+1]],2],
       pch=21,bg='gray',main=colnames(counts)[i],ylim=range(prop.quants[grep('p.true1',colnames(samples.mixed))[(plot.help[i]+1):plot.help[i+1]],]),ylab='pollen prop',xlab='biomass')
  
  ciEnvelope(x = 1:145,
             ylo = prop.quants[grep('p.true1',colnames(samples.mixed))[(plot.help[i]+1):plot.help[i+1]],1],
             yhi = prop.quants[grep('p.true1',colnames(samples.mixed))[(plot.help[i]+1):plot.help[i+1]],3],
             col = 'lightblue') 
  points(1:145,prop.quants[grep('p.true1',
                                colnames(samples.mixed))[(plot.help[i]+1):plot.help[i+1]],2],
         pch=21,bg='gray',cex=.25)
  points(biomass,counts[,i]/total_counts,cex=.8,pch=19,col='blue')
  
  #points(1:145,samples.mixed[1158,grep('p.true1',
  #                               colnames(samples.mixed))[(plot.help[i]+1):plot.help[i+1]]])
  #points(biomass[95],counts[95,i]/total_counts[95],cex=.8,pch=19,col='red')
}
if(DRAW == TRUE) dev.off()

par(mfrow=c(4,4))
for(i in 1:ncol(samples.mixed)){
  plot(samples.mixed[,i],typ='l',main=colnames(samples.mixed)[i])
}
}
  
if(FALSE){
beta.names <- rep(rep(colnames(Y),each=5),2)
beta.i <- grep('beta',colnames(samples.mixed))

pdf(paste0('beta.posteriors.group',group_rm,'.pdf'))
par(mfrow=c(4,4))
for(i in seq_along(beta.i)){
  plot(samples.mixed[,beta.i[i]],typ='l',main=paste('beta',beta.i[i],'pollen',beta.names[i]))
  hist(samples.mixed[,beta.i[i]],main=paste('beta',beta.i[i]),freq=F)
}
dev.off()
}

}

