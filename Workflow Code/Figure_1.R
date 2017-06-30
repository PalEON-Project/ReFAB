biomassCI <- list()
nItsSave = 10000
for(i in 1:43){
  locn <- site.names[i]
  site_number = unique(x.meta[x.meta$site.name == locn,1])
  ten_count_use = ten.count[which(x.meta$site.id == site_number), ]
  
  Y = as.matrix(ten_count_use)
  if(length(Y)>21 & nrow(Y) > 15 &
     max(x.meta[x.meta$site.name == locn,'age_bacon'])>8000 & 
     min(x.meta[x.meta$site.name == locn,'age_bacon'])<2000){
    if(file.exists(paste0('~/ReFAB/MCMC Samples/samplesList_',locn,'.Rda'))){
      load(file = paste0('~/ReFAB/MCMC Samples/samplesList_',locn,'.Rda'))
      biomassCI[[i]] <- apply(samplesList[round(nItsSave/5):nItsSave,1:99],2,quantile,c(0.025,0.5,0.975))
      names(biomassCI)[i] <- x.meta[x.meta$site.name == locn,'site.id'][1]
    }
    
  }
  
}

breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)))

bio.quant <- apply(matrix(unlist(lapply(biomassCI,function(x) x[2,])),19,99,byrow = TRUE),2,quantile,c(0.025,0.5,0.975))
data_binned <-  cut(bio.quant[2,], c(breaks), include.lowest = FALSE, labels = FALSE)
diff.list <- diff.min <- diff.max <- list()
for(i in 1:43){
  if(length(biomassCI[[i]])>1){
    diff.list[[i]] <- diff(rev(biomassCI[[i]][2,]),lag=1)
    diff.min[[i]] <- diff(rev(biomassCI[[i]][1,]),lag=1)
    diff.max[[i]] <- diff(rev(biomassCI[[i]][3,]),lag=1)
  }else{
    diff.list[[i]] <- rep(0,98)
    diff.min[[i]] <- rep(0,98)
    diff.max[[i]] <- rep(0,98)
  }
}

diff.mat <- matrix(unlist(diff.list),43,98,byrow=TRUE)
diff.mat1 <- matrix(unlist(diff.min),43,98,byrow=TRUE)
diff.mat2 <- matrix(unlist(diff.max),43,98,byrow=TRUE)
diff.prop <- diff.prop1 <- diff.prop2 <- numeric(ncol(diff.mat))
for(i in 1:ncol(diff.mat)){
  diff.prop[i] <- length(which(diff.mat[,i]>0))/19
  diff.prop1[i] <- length(which(diff.mat1[,i]>0))/19
  diff.prop2[i] <- length(which(diff.mat2[,i]>0))/19
}
length(which(diff.mat[,1]!=0))

pdf(paste0('average.biomass',Sys.Date(),'.pdf'))
#quartz()
zones=matrix(c(2,1), ncol=1, byrow=TRUE)
layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
par(mar=c(5,4,1,1))
plot(seq(100,9900,100),bio.quant[2,],
     col=colors[data_binned],pch=19,xlim=c(10000,-10),
     ylim=c(0,150),xlab = 'Year Before Present',ylab='Biomass (Mg/ha)')
for(i in 1:43){  
  if(length(biomassCI[[i]])>1){
    lines(seq(100,9900,100),biomassCI[[i]][2,],col='gray')
  }
}
points(seq(100,9900,100),bio.quant[2,],
       col=colors[data_binned],pch=19)
par(mar=c(0,4,1,1))
barplot(height=diff.prop, ylim=c(0, 1), 
        space=0,xlim=c(0,100),ylab = 'Proportion Increasing',
        cex.lab=.6, cex.axis = .5)
abline(h=.5)
dev.off()




barplot(height=diff.prop, ylim=c(0, 1), 
        space=0,xlim=c(0,100),ylab = 'Proportion Increasing',
        cex.lab=.6, cex.axis = .5)

segments(seq(.5,97.5,1),diff.prop1,seq(.5,97.5,1),diff.prop2,col='red',lwd=2)

