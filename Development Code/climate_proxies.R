blue = col2rgb("blue")
alphablue = rgb(blue[1],blue[2],blue[3],75,max=255)

#plots a confidence interval around an x-y plot (e.g. a timeseries)
ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}

co2<-read.table('~/Downloads/antarctica2015co2composite.txt',stringsAsFactors = FALSE)
head(co2)

inso<-read.table('~/Downloads/j_45north.txt',stringsAsFactors = FALSE)
head(inso)

temp<-read.csv('Marcott_fig1_tempA.csv',header = FALSE,
               stringsAsFactors = FALSE,na.strings = 'NaN')
temp<-na.omit(temp)
temp<-temp[-1,]

breaks <-  c(seq(0,40,10),seq(80,160,40))
colors <- rev(terrain.colors(length(breaks)-1))
data_binned1 <-  cut(calc.all1[,27], c(breaks), include.lowest = FALSE, labels = FALSE)

dbs<-numeric(max(plot.which.keep))

dbs[plot.which.keep] <- data_binned1[order(calc.all1[,27])]

names(biomassCI)<-unique(x.meta[,1])

pdf('climate_proxies_biomass.pdf')
par(mfrow=c(4,1),mar=c(2,4,2,2))
     
plot(co2$V1[-1],co2$V2[-1],typ='l',lwd=2,main='CO2 (PPM)',
     xlab='Age BP',xlim=c(0,10000),ylab='CO2 Concentration (ppm)',
     ylim=c(250,340))
#ciEnvelope(co2$V1[-1],as.numeric(co2$V2[-1])-10,
#           as.numeric(co2$V2[-1])+10,col=alphablue)

plot(inso$V1[-1]*1000,inso$V2[-1],typ='l',lwd=2,
     main = 'Summer Insolation',
     xlab='Age BP',xlim=c(0,10000),
     ylab='Summer Insolation', ylim=c(9.7,9.715))

plot(temp$V1[-1],temp$V2[-1],typ='l',lwd=2,
     main = 'Global Temperature (C)',
     xlab='Age BP',xlim=c(0,10000),
     ylab='Global Temperature (C)')
ciEnvelope(as.numeric(temp$V1),na.omit(as.numeric(temp$V2)-as.numeric(temp$V2)*as.numeric(temp$V3)),
           na.omit(as.numeric(temp$V2)+as.numeric(temp$V2)*as.numeric(temp$V3)),col=alphablue)

#run stuff in pick.sites.R

	  plot(seq(100,9900,100),biomassCI[[i]][2,],cex=.1,ylim=c(0,150),xlim=c(-10,10000),ylab="Biomass (Mg / Ha)",xlab="Years BP",main='Biomass')
      #title(x.meta[x.meta[,1]== site.id.list[i],5][1],outer=TRUE)
      #axis(3,at=seq(0,10000,1000),labels=seq(0,10000,1000),padj=1)
for(i in plot.which.keep){
      #ciEnvelope(seq(100,9900,100),biomassCI[[i]][1,],biomassCI[[i]][3,],col=colors[dbs[i]])
      if(length(biomassCI[[i]])>2){
	      plot.col <- cut(calc.all1[calc.all1$SiteID==names(biomassCI)[i],27], c(breaks), include.lowest = FALSE, labels = FALSE)
	      points(seq(100,9900,100),biomassCI[[i]][2,],cex=.8,pch=16,col=colors[plot.col])
      }
}
dev.off()
