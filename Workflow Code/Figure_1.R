


biomassCI <- list()
all.samps<- numeric(100)
nItsSave = 10000
for(i in 1:62){
  locn <- names(how.many)[i]
  site_number = unique(x.meta[x.meta$site.name == locn,1])
  ten_count_use = ten.count[which(x.meta$site.id == site_number), ]
  
  Y = as.matrix(ten_count_use)
  if(length(Y)>21 & nrow(Y) > 15 &
     max(x.meta[x.meta$site.name == locn,'age_bacon'])>8000 & 
     min(x.meta[x.meta$site.name == locn,'age_bacon'])<2000){
    if(file.exists(paste0('~/ReFAB/MCMC Samples/samplesList_',locn,'.Rda'))){
      load(file = paste0('~/ReFAB/MCMC Samples/samplesList_',locn,'.Rda'))
      biomassCI[[i]] <- apply(samplesList[round(nItsSave/5):nItsSave,1:100],2,quantile,c(0.025,0.5,0.975))
      names(biomassCI)[i] <- x.meta[x.meta$site.name == locn,'site.id'][1]
      all.samps <- rbind(all.samps,samplesList[round(nItsSave/5):nItsSave,1:100])
      }
    
  }
  
}

breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)))

#bio.quant <- apply(matrix(unlist(lapply(biomassCI,function(x) x[2,])),61,100,byrow = TRUE),2,quantile,c(0.025,0.5,0.975))
bio.quant <- apply(all.samps,2,quantile,c(0.025,0.5,0.975))
data_binned <-  cut(bio.quant[2,], c(breaks), include.lowest = FALSE, labels = FALSE)


diff.mat.all<- t(diff(t(all.samps[,100:1])))

how.much <- matrix(NA,8001,99)
for(i in 1:99){
  for(r in 1:8001){
    how.much[r,i] <- length(which(diff.mat.all[seq(1+r,488062,8001),i]>0))/61
  }
}

save.how.much <- apply(how.much,2,quantile,c(0.025,0.5,0.975))


pdf(paste0('average.biomass',Sys.Date(),'.pdf'))
#quartz()
zones=matrix(c(2,1), ncol=1, byrow=TRUE)
layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
par(mar=c(5,4,1,1))
plot(seq(100,10000,100),bio.quant[2,],
     col=colors[data_binned],pch=19,xlim=c(10000,-10),
     ylim=c(0,150),xlab = 'Year Before Present',ylab='Biomass (Mg/ha)')
for(i in 1:length(biomassCI)){  
  if(length(biomassCI[[i]])>1){
    lines(seq(100,10000,100),biomassCI[[i]][2,],col='gray')
  }
}
points(seq(100,10000,100),bio.quant[2,],
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



