splines_plot <- function(samples.mixed,Y,biomass,bMax){
  
  new.biomass <- 1:bMax
  Z.new = matrix(0,nrow=length(new.biomass),ncol=5)
  for(i in 1:length(new.biomass)){
    u_given <- new.biomass[i]
    Z.new[i,] = bs_nimble(u_given, u=u, N0 = rep(0, (length(u)-1)),
                          N1 = rep(0, (length(u))), 
                          N2 = rep(0, (length(u)+1)), 
                          N3 = rep(0, (length(u)+2)))
  }
  
  i.beta1 <- grep(pattern = 'beta1',colnames(samples.mixed))
  i.beta2 <- grep(pattern = 'beta2',colnames(samples.mixed))
  
  beta1 = matrix(colMeans(samples.mixed[,i.beta1]),ncol(Z.new),ncol(Y))
  beta2 = matrix(colMeans(samples.mixed[,i.beta2]),ncol(Z.new),ncol(Y))
  
  alpha <- linexp(Z.new%*%beta1,J=bMax,I=ncol(Y))
  beta <- linexp(Z.new%*%beta2,J=bMax,I=ncol(Y))
  
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
  props.pred <- prop.table(as.matrix(Y),1)
  
  #outlier <- as.numeric(names(sort(abs(biomass-colMeans(samples.pred[,grep('b',colnames(samples.pred))])))))[90]
  
  #pdf(paste0('outlier.splines.group.',group_rm,'.pdf'))
  #pdf('splines.together.pdf')
  #pdf('splines.new_data.pdf')
  par(mfrow=c(2,2))
  for(i in 1:ncol(Y)){
    if(FALSE){#i %in% seq(1,22,3)
      plot.new()
      legend('center',c('PLS Biomass','ReFAB Biomass','PLS Biomass, non outlier'),pch=c(19,19,1),
             col=c('blue','red','black'),title='R2 Outliers')
    }
    plot(1:bMax, mean.betas.ex[,i],
         ylim=c(0,max(mean.betas.ex[,i] + sd.betas.ex[,i]*2,
                      mean.betas.ex[,i] - sd.betas.ex[,i]*2, 
                      props[,i])),
         main = colnames(Y)[i],pch=19,
         ylab = 'Pollen Proporation',
         xlab = 'Biomass (Mg/ha)')
    abline(v=u)
    # points(1:bMax, mean.betas.ex[,i],col='red')
    ciEnvelope(x = 1:bMax, yhi = mean.betas.ex[,i] + sd.betas.ex[,i]*2,
               ylo = mean.betas.ex[,i] - sd.betas.ex[,i]*2,col = adjustcolor('blue',alpha=.5))
    
    #ciEnvelope(x = 1:bMax, yhi = mean.betas.ex[,i] + sd.betas.ex[,i]*2,
    #           ylo = mean.betas.ex[,i] - sd.betas.ex[,i]*2,col = alphaorange)
    points(biomass,props[,i])
    #points(biomass[bimodal_sites],props[,i][bimodal_sites],pch=19,col='blue')
    #points(colMeans(samples.pred)[bimodal_sites],props[,i][bimodal_sites],pch=19,col='red')
    
    textxy(biomass,props[,i],1:length(biomass),offset = 0)
    #textxy(colMeans(samples.pred)[bimodal_sites],props[,i][bimodal_sites],bimodal_sites,offset = 0)
    
    #points(biomass.pred,props.pred[,i],col='red',pch=19)
    #legend('topright',c('bMax = 150','bMax = 143'),pch=c(19,19),col=c(alphablue,alphaorange))
    ###Validation outlier
  }
  #dev.off()
}

splines_plot(samples.mixed,Y,biomass,bMax)
