#####
##### Second Derivative #####
#####

calc.second.deriv <- function(biomassCI,h,second.deriv){
  for(i in 1:length(biomassCI)){
    if(length(biomassCI[[i]])>10){
      T <- dim(biomassCI[[i]][,1:20])[2] - h
      t <- h + 1 
      biomassCI[[i]]<-(biomassCI[[i]][,1:20]) #log or no log here?
      second.deriv[[i]]<-sum(((biomassCI[[i]][2,(t:T)+h]-2*biomassCI[[i]][2,(t:T)]+biomassCI[[i]][2,(t:T)-h])/((h*100)^2))^2)
      
    }else{
      second.deriv[[i]]<-NA
    }
  }
  return(second.deriv)
}

pdf(paste('second.deriv.map',Sys.Date(),'.pdf'))
for(i in c(2)){
  h=1
  second.deriv<-calc.second.deriv(biomassCI=biomassCI,h=h,second.deriv=list())
  
  #hist(unlist(second.deriv))
  names(second.deriv) <- datID[,'name']#names(how.many)
  second.deriv.unlist <- unlist(second.deriv)
  
  second.deriv.unlist <- na.omit(second.deriv.unlist)
  
  
  long.keep <- list()
  lat.keep <- list()
  settlebiomass <- list()
  site.name.list <- list()
  site.id.list <- list()
  for(i in 1:length(second.deriv.unlist)){
    long.keep[[i]] <- x.meta[x.meta[,'site.name']==names(second.deriv.unlist)[i],'long'][1]
    lat.keep[[i]] <- x.meta[x.meta[,'site.name']==names(second.deriv.unlist)[i],'lat'][1]
    #settlebiomass[[i]] <- cast.x[cast.x$site.id == names(second.deriv.unlist)[i],ncol(cast.x)][1]
    site.name.list[[i]] <- x.meta[x.meta[,'site.name']==names(second.deriv.unlist)[i],'site.name'][1]
    site.id.list[[i]] <- x.meta[x.meta[,'site.name']==names(second.deriv.unlist)[i],'site.id'][1]
  }
  
  tot.change.mat <- data.frame(unlist(site.id.list),unlist(long.keep),unlist(lat.keep),second.deriv.unlist)
  colnames(tot.change.mat)<-c('site.id','long','lat','second_deriv1000')
  #write.csv(tot.change.mat,file='Second_Deriv_with_meta.csv')
  
  breaks <-  c(0,quantile(probs=c(.05,.25,.5,.75,.975),second.deriv.unlist,na.rm=TRUE),max(second.deriv.unlist)+1)
  
  colors <- colorRampPalette(c("white","yellow",'green','blue','purple'))(length(breaks)-1)
  
  data_binned <-  cut(second.deriv.unlist, c(breaks), include.lowest = FALSE, labels = FALSE)
  
  map('state', xlim=c(-98,-81), ylim=c(41.5,50),bg='grey')
  points(long.keep, lat.keep, pch=19,
         cex=1.3, col=colors[data_binned],lwd=.2)
  title(paste("Second Derivative Sum Point Estimates h = ",h*100))
  
  plotInset(-90,47,-82.5,50,
            expr={
              hist(data_binned,col=colors,xaxt="n",xlab=NA,
                   ylab=NA,main=NA,cex.lab=.5,cex.axis=.5,breaks=0:length(unique(data_binned)))
              axis(side=1,signif(breaks,digits=2),at = 0:length(unique(data_binned)),cex.axis = .5,las=2,line=0)
              mtext(side = 1, "Squared Sum", line = 1.5,cex=.5)
              mtext(side = 2, "Frequency", line = 1.7,cex=.5)
            })
  
}
dev.off()

data.dir = c("/Users/paleolab/babySTEPPS/Data/")
biomass_dat_est <- read.csv(paste0(data.dir,"biomass_prediction_v0.9-7_bam.csv"))

lat.long.reg.df = data.frame(biomass_dat_est$x,biomass_dat_est$y)
colnames(lat.long.reg.df) = c('x', 'y')

coordinates(lat.long.reg.df) <- ~ x + y
proj4string(lat.long.reg.df) <-  CRS('+init=epsg:3175')

albers <- spTransform(lat.long.reg.df,CRS('+proj=longlat +ellps=WGS84'))
albers <- as.matrix(data.frame(albers))

breaks <-  c(0,quantile(probs=c(.05,.25,.5,.75,.975),second.deriv.unlist,na.rm=TRUE),max(second.deriv.unlist)+1)
colors <- colorRampPalette(c("white","yellow",'green','blue','purple'))(length(breaks)-1)
data_binned <-  cut(second.deriv.unlist, c(breaks), include.lowest = FALSE, labels = FALSE)


pdf(paste0('wiggle.map',Sys.Date(),'.pdf'))
map('state', xlim=c(-98,-81), ylim=c(41.5,50),col='black',bg='lightgrey')
#points(albers[which(biomass_dat_est$Hemlock>1),],pch=11,cex=.1,col='darkgrey')
#points(albers[which(biomass_dat_est$Beech>1),],pch=19,cex=.26)
#points(unlist(long.keep), unlist(lat.keep), pch=19,
#       cex=1, col='darkgrey',lwd=.2)
#title(paste("Second Derivative Sum Point Estimates h = ",h*100))
title('Time Series Map')

for(i in 1:97){
  if(length(biomassCI[[i]])>1){
    long.use <- lon[[i]]#x.meta[x.meta[,'site.name']==names(how.many)[i],'long'][1]
    lat.use<- lat[[i]]#x.meta[x.meta[,'site.name']==names(how.many)[i],'lat'][1]
    color.use <- data_binned[i]#which(names(second.deriv.unlist)==names(how.many)[i])
    
    plotInset(long.use-3,lat.use-1,long.use,lat.use+1,
              expr={  
                plot(seq(100,2000,100),
                     biomassCI[[i]][2,],
                     typ='l',xaxt = 'n',
                     yaxt = 'n',xlab = NA,ylab=NA,xlim=c(10000,0),
                     bty = 'n',ylim=c(0,150),col=colors[color.use],lwd=3)
              })
  }
}
plotInset(-90,47.5,-82.5,50.5,
          expr={
            hist(data_binned,col=colors,xaxt="n",xlab=NA,
                 ylab=NA,main=NA,cex.lab=.5,cex.axis=.5,breaks=0:length(unique(data_binned)))
            axis(side=1,c(0,1),
                 at = c(0,length(unique(data_binned))),
                 cex.axis = .5,las=2,line=0)
            mtext(side = 1, "Squared Sum", line = 1,cex=1)
            mtext(side = 2, "Frequency", line = 1.7,cex=.5)
          })
dev.off()

pdf(paste0('wiggle.map.biomass.colors',Sys.Date(),'.pdf'))

#quartz()
map('state', xlim=c(-98,-81), ylim=c(41.5,50),col='black',bg='lightgrey')
points(albers[which(biomass_dat_est$Hemlock>1),],pch=11,cex=.1,col='darkgrey')
points(albers[which(biomass_dat_est$Beech>1),],pch=19,cex=.26)
#points(unlist(long.keep), unlist(lat.keep), pch=19,
#       cex=1, col='darkgrey',lwd=.2)
#title(paste("Second Derivative Sum Point Estimates h = ",h*100))
title('Time Series Map')

for(i in 1:62){
  if(length(biomassCI[[i]])>1){
    long.use <- x.meta[x.meta[,'site.name']==names(how.many)[i],'long'][1]
    lat.use<- x.meta[x.meta[,'site.name']==names(how.many)[i],'lat'][1]
    #color.use <- data_binned[which(names(second.deriv.unlist)==names(how.many)[i])]
    
    breaks <-  c(seq(0,50,10),seq(75,200,25))
    colors <- rev(terrain.colors(length(breaks)))
    data_binned <-  cut(rev(biomassCI[[i]][2,]), c(breaks), include.lowest = FALSE, labels = FALSE)
    
    
    plotInset(long.use-3,lat.use-1,long.use,lat.use+1,
              expr={  
                plot(seq(100,10000,100),
                     biomassCI[[i]][2,],
                     typ='b',xaxt = 'n',
                     yaxt = 'n',xlab = NA,ylab=NA,xlim=c(10000,0),
                     bty = 'n',ylim=c(0,150),col=rev(colors[data_binned]),cex=.3,pch=19)
              })
  }
}

dev.off()
