


biomassCI <- diff.median <-  list()
prob.of.inc <- matrix(NA,62,99)
all.samps<- numeric(100)
nItsSave = 10000
hem.is <- list()
lat <- long <- name.keep <- list()
for(i in 1:62){
  locn <- names(how.many)[i]
  site_number = unique(x.meta[x.meta$site.name == locn,1])
  ten_count_use = ten.count[which(x.meta$site.id == site_number), ]
  hem.is[[i]]<-sum(ten_count_use[,11])
  Y = as.matrix(ten_count_use)
  if(length(Y)>21 & nrow(Y) > 15 &
     max(x.meta[x.meta$site.name == locn,'age_bacon'])>8000 & 
     min(x.meta[x.meta$site.name == locn,'age_bacon'])<2000){
    if(file.exists(paste0('~/ReFAB/MCMC Samples/samplesList_',locn,'.Rda'))){
      load(file = paste0('~/ReFAB/MCMC Samples/samplesList_',locn,'.Rda'))
      biomassCI[[i]] <- apply(samplesList[round(nItsSave/5):nItsSave,1:100],2,quantile,c(0.025,0.5,0.975))
      names(biomassCI)[i] <- x.meta[x.meta$site.name == locn,'site.id'][1]
      all.samps <- rbind(all.samps,samplesList[round(nItsSave/5):nItsSave,1:100])
      diff.median[[i]] <- apply(diff(t(samplesList[round(nItsSave/5):nItsSave,1:100])),1,quantile,c(0.5))
      test <- diff(t(samplesList[round(nItsSave/5):nItsSave,1:100]))
      lat[[i]] <- x.meta[x.meta$site.name == locn,'lat'][1]
      long[[i]] <- x.meta[x.meta$site.name == locn,'long'][1]
      name.keep[[i]] <- locn
      for(t in 1:99){
        prob.of.inc[i,t] <- length(which(test[t,]<0))/8001
      }
      
      }
  }
}

map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(unlist(long),unlist(lat), pch=19,
       cex=1.1,lwd=.2)
text(unlist(long),unlist(lat),unlist(name.keep),cex=.5)

pdf('boxplots.prob.of.inc.pdf')
boxplot(prob.of.inc,xlim=c(100,1))
abline(h=.5)
title('Probability of Increase')
dev.off()

pdf('barplot.of.prop.inc.pdf')
diff.median.mat <- matrix(unlist(diff.median),99,62)
bars <- numeric(99)
for(t in 1:99){
  bars[t] <- length(which(diff.median.mat[t,]<0))/62
}

barplot(rev(bars))
abline(h=.5)
title('Proportion Median Increasing')
dev.off()



breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)))

#bio.quant <- apply(matrix(unlist(lapply(biomassCI,function(x) x[2,])),61,100,byrow = TRUE),2,quantile,c(0.025,0.5,0.975))
bio.quant <- apply(all.samps,2,quantile,c(0.025,0.5,0.975))

diff.mat.all<- t(diff(t(all.samps[,100:1])))

how.much <- matrix(NA,8001,99)
median.bucket <- matrix(NA,8001,100)
for(i in 1:99){
  for(r in 1:8001){
    how.much[r,i] <- length(which(diff.mat.all[seq(1+r,488062,8001),i]>0))/61
  }
}

#calculate the probability of increase for each site then make a boxplot of all of the sites probability of increase at each time.
#proportion of medians that are greater than zero
save.how.much <- apply(how.much,2,quantile,c(0.025,0.5,0.975))

for(i in 1:100){
  for(r in 1:8001){
     median.bucket[r,i] <- quantile(all.samps[seq(1+r,488062,8001),i],c(.5))
  }
}
save.median.bucket <- apply(median.bucket,2,quantile,c(0.025,0.5,0.975))

data_binned <-  cut(save.median.bucket[2,], c(breaks), include.lowest = FALSE, labels = FALSE)

pdf(paste0('average.biomass',Sys.Date(),'.pdf'))
#quartz()
zones=matrix(c(2,1), ncol=1, byrow=TRUE)
layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
par(mar=c(5,4,1,1))
plot(seq(100,10000,100),save.median.bucket[2,],
     col=colors[data_binned],pch=19,xlim=c(10000,-10),
     ylim=c(0,150),xlab = 'Year Before Present',ylab='Biomass (Mg/ha)')
for(i in 1:length(biomassCI)){  
  if(length(biomassCI[[i]])>1){
    if(hem.is[[i]]>10){
      lines(seq(100,10000,100),biomassCI[[i]][2,],col='darkgreen')
    }else{
      lines(seq(100,10000,100),biomassCI[[i]][2,],col='red')
    }
    
  }
}
ciEnvelope(x = seq(1,9900,100),ylo = save.median.bucket[1,],yhi =save.median.bucket[3,],col='gray')
points(seq(100,10000,100),save.median.bucket[2,],
       col=colors[data_binned],pch=19)
#points(seq(100,10000,100),bio.quant[1,], pch = 19)
#points(seq(100,10000,100),bio.quant[3,], pch=19)
par(mar=c(0,4,1,1))
barplot(height=save.how.much[2,], ylim=c(0, 1), 
        space=0,xlim=c(0,100),ylab = 'Proportion Increasing',
        cex.lab=.6, cex.axis = .5)

segments(seq(.5,98.5,1),save.how.much[1,],seq(.5,98.5,1),save.how.much[3,],col='red',lwd=1)

abline(h=.5)
dev.off()




site_number = unique(x.meta[x.meta$site.name == "Wood Lake",1])
ten_count_use = ten.count[which(x.meta$site.id == site_number), ]
#capitola 38
#kelly's 34
#wood 33

lat.long.reg <- cbind(unlist(long),unlist(lat))
lat.long.reg.df = data.frame(lat.long.reg)
colnames(lat.long.reg.df) = c('x', 'y')

coordinates(lat.long.reg.df) <- ~ x + y
proj4string(lat.long.reg.df) <- CRS('+proj=longlat +ellps=WGS84')

albers <- spTransform(lat.long.reg.df, CRS('+init=epsg:3175'))
albers <- as.matrix(data.frame(albers))




