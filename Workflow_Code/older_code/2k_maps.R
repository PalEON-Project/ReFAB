head(x.meta)

ids <- (unique(x.meta$site.id))

mip.meta <- list()

for(i in 1:length(ids)){
  mip.meta[[i]] <- x.meta[x.meta$site.id==ids[i],]
  if(max(mip.meta[[i]]$age_bacon) < 600 | nrow(mip.meta[[i]]) < 4) mip.meta[[i]] <- NULL
}

mip.meta.use <- do.call(rbind,mip.meta)

length(unique(mip.meta.use$site.id))

#"","name","ID","sigma","group"
datID <- data.frame(unique(mip.meta.use$site.name),1:length(unique(mip.meta.use$site.name)),.12,NA)
colnames(datID) <- c('name','ID','sigma','group')

write.csv(datID,file.path('Cross_Validation','paleon.dataID.csv'))

biomassCI <- lat <- lon <- list()
counter <- 0


for(i in 1:97){
  locn <- datID[i,'name']
  locnClean <- gsub(' ', '-', locn)
  
  lat[[i]] <- x.meta[x.meta$site.name==locn,'lat'][1]
  lon[[i]] <- x.meta[x.meta$site.name==locn,'long'][1]
  
  file.name <- file.path('~/ReFAB/samples.paleon',
                         paste0('samplesList_workInfo_',
                                locnClean,
                                'Sigma0.12GroupNA.Rda.Rda'))
  if(file.exists(file.name)){
    load(file.name)
    biomassCI[[i]] <- apply(samplesList[40:200,1:20],2,quantile,c(0.025,0.5,0.975),na.rm=TRUE)
  }else{
    biomassCI[[i]] <- NULL
    counter <- c(counter,i)
  }
  
}

par(mfrow=c(1,1))
for(i in 2:97) lines(biomassCI[[i]][2,])


pdf(paste0('SiteDiagnositcs.mip.pdf'))

#quartz()
for(i in order(second.deriv.unlist)){


locn <- datID[i,'name']
site_number = unique(x.meta[x.meta$site.name == locn,1])
keep.dataset.id <- unique(x.meta[x.meta$site.id==site_number,4])

par(mfrow=c(1,1))

map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(lon[[i]], lat[[i]],pch=19,cex=1.5)
title(locn)

breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)))

#browser()

bio.quants <- biomassCI[[i]][,1:12]#apply(samplesList[round(nItsSave/5):nItsSave,1:(maxAge/100)],2,quantile,c(0.025,0.5,0.975))

data_binned <-  cut(rev(bio.quants[2,]), c(breaks), include.lowest = FALSE, labels = FALSE)

fig.mat <- matrix(1,28,1)
fig.mat[1:6,]<-1
fig.mat[7:28,]<-seq(2,23,1)
#control.pts<-read.csv(text=getURL('https://raw.githubusercontent.com/PalEON-Project/stepps-baconizing/master/data/bacon_age_markers_v1.csv'))

def.par <- par(no.readonly = TRUE)
layout(fig.mat)
par(mar = c(0, 4, 0, 3), oma = c(5, 0, 4, 0),tcl=.5)

plot(bio.quants[2,], xlim = c(maxAge/100,0), ylim = c(0,150), xaxt='n',
     xlab = 'Years BP', ylab = 'Biomass Mg/ha', col = 'white')
axis(side = 3, at = rev(seq(0,maxAge/100, round((maxAge/100)/6))), labels = rev(seq(0,maxAge,round(maxAge/6))))
ciEnvelope(x=1:(maxAge/100), ylo = bio.quants[1,],yhi = bio.quants[3,],col = 'grey')
points(bio.quants[2,],cex=1.1,pch=16,col = rev(colors[data_binned]))
rug(x.meta[x.meta[,1]==site_number,]$age_bacon/100,lwd=2)
rug(control.pts[which(control.pts[,2]%in%keep.dataset.id),]$geo_age/100,lwd=3,col="red")
#bMax=145
#points(age_index,seq(5, bMax-5, by = 2)[apply(out,2,which.max)])
#points(0,unique(x.meta[x.meta$site.name == locn,'SettleBiomass']),pch=19,col='purple',cex=2)
#legend('topleft','Mx.Lik.',pch=1)

ten_count_use = ten.count[which(x.meta$site.id == site_number), ]
Y = as.matrix(ten_count_use)
Y[which(is.na(Y))] <-0
prop.use <- prop.table(as.matrix(Y),margin=1)    

for(p in rev(match(names(sort(colMeans(prop.use))),colnames(prop.use)))){
  prop.plot<- cbind(as.vector(x.meta[which(x.meta$site.id==site_number),]$age_bacon),as.matrix(prop.use[,p]))      	
  prop.plot<-prop.plot[order(prop.plot[,1]),]
  plot(x=prop.plot[,1],y=prop.plot[,2],type="l",xlim=c(maxAge,-10),
       ylim=c(0,max(prop.use[,p])),ylab=NA,yaxt='n', xaxt='n')
  #axis(2,at=signif(max(prop.use)/2),labels=signif(max(prop.use[,p])/2,digits=1))
  ciEnvelope(prop.plot[,1],rep(0,nrow(prop.use)),prop.plot[,2],col="darkblue")
  legend('topleft',colnames(prop.use)[p])
  #legend('topleft',paste(signif(mean(prop.use[,p]),digits=2)),bty="n")
}

}

dev.off()


####
#### Biomass Maps ####
####

breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)-1))

pdf('biomass.maps.mip.pdf')
for(t in 20:1){
  biomass.use <- unlist(lapply(biomassCI,function(x){x[2,t]}))
  data_binned <-  cut(biomass.use, c(breaks), include.lowest = FALSE, labels = FALSE)
  
  map('state', xlim=c(-98,-81), ylim=c(41.5,50))
  points(unlist(lon),unlist(lat), pch=21,
         cex=1.1, bg=colors[data_binned],lwd=.2)
  title(paste0("Biomass Point Estimates t = ",t*100))
  plotInset(-90,47,-82.5,50,
            expr={
              hist(data_binned,col=colors,xaxt="n",xlab=NA,
                   ylab=NA,main=NA,cex.lab=.5,cex.axis=.5)
              axis(side=1,breaks,at = seq(1,12,1),cex.axis = .5,las=2,line=0)
              mtext(side = 1, "Biomass (Mg/ha)", line = 1.5,cex=.5)
              mtext(side = 2, "Frequency", line = 1.7,cex=.5)
            })
}
dev.off()

####
#### Difference Maps ####
####

breaks <- seq(-10,10,2)
colors <- colorRampPalette(c('red',"white",'blue'))(length(breaks)-1)

biom.med <- lapply(biomassCI,function(x){x[2,]})
biomass.mat <- do.call(rbind,biom.med)

#quartz()
#pdf(paste0('difference.maps.',Sys.Date(),'.pdf'))
pdf('mip.diff.maps.100.pdf')
par(mfrow=c(3,4), oma = rep(.1,4))
b <- 1
for(r in rev(seq(2,12,b))){ #just thousand year time bins
  #if(r-b==0) b = b-5
  x <- biomass.mat[,r] - biomass.mat[,(r-b)]
  data_binned <-  cut(x, c(breaks), include.lowest = FALSE, labels = FALSE)
  map('state', xlim=c(-97,-83), ylim=c(41.75,49),mar=rep(0,4))
  points(unlist(lon),unlist(lat),pch=21,bg=colors[data_binned],col='lightgray')
  title(paste(r*100, '-', (r-b)*100))
}
plot.new()
legend('center',legend=as.character(breaks[1:length(breaks)-1]),pch=rep(19,length(breaks)-1),col=colors)
dev.off()



corrplot(cor(t(biomass.mat[,])),order = 'FPC')

ser <- which(rowSums(cor(t(biomass.mat)))<0)

map('state', xlim=c(-97,-83), ylim=c(41.75,49),mar=rep(0,4))
points(unlist(lon)[ser],unlist(lat)[ser],pch=21,bg='red',col='lightgray')
points(unlist(lon)[-ser],unlist(lat)[-ser],pch=21,bg='blue',col='lightgray')
#might just be picking the places with not that much data



