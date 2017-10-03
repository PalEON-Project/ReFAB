library(nimble)
library(splines)
library(maps)
library(plyr)
library(oce)
library(RCurl)
library(corrplot)

data.dir = "/Users/paleolab/babySTEPPS/Data/"
fig.dir = "/Users/paleolab/babySTEPPS/Figures/"
model.dir = "/Users/paleolab/babySTEPPS/Code/"

setwd("/Users/paleolab/babySTEPPS/")
load("add.bacon2.Rdata")
load("2016-05-31nimble.betas.Rdata")
source(paste0(model.dir,"bs_nimble.R"))

#plots a confidence interval around an x-y plot (e.g. a timeseries)
ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}

x = new.pol1[new.pol1$age_bacon>=200,]
x = x[x$age_bacon<=10000,]
ten.count = x[,7:26]
ten.count.save = ten.count
ten.count = round(ten.count.save)

x.meta = x[,c('SiteID','LatitudeNorth',"LongitudeWest","dataset.id","sitename","age_bacon")]

load(file="biomass.CI.Rdata")

list.buddy <- seq(1,length(unique(x.meta[,1])),1)
lat.save <- list()
lon.save <- list()
for(i in 1:length(unique(x.meta[,1]))){
	lat.save[[i]]<-unique(x.meta[x.meta[,1]==unique(x.meta[,1])[i],2])
	lon.save[[i]]<-unique(x.meta[x.meta[,1]==unique(x.meta[,1])[i],3])
}
#lat.save[[1]]<-45.78650
lat.save <- unlist(lat.save)
lon.save <- unlist(lon.save)

site.id.list <- unique(x.meta[,1])
dataset.id.list <- unique(x.meta[,4])
name.list <-unique(x.meta[,5])

control.pts<-read.csv(text=getURL('https://raw.githubusercontent.com/PalEON-Project/stepps-baconizing/master/data/bacon_age_markers_v1.csv'))

plot.which=numeric(500)

############# Conditions for picking sites

for(i in list.buddy[order(lat.save)]){#
	if(
	min(x.meta[x.meta[,1]== site.id.list[i],]$age_bacon)<1000 & 
	max(x.meta[x.meta[,1]== site.id.list[i],]$age_bacon)>9000)
	#min(control.pts[which(control.pts[,1]%in%unique(x.meta[x.meta[,1]==site.id.list[i],4])),]$geo_age)<2000 &
	#max(control.pts[which(control.pts[,1]%in%unique(x.meta[x.meta[,1]==site.id.list[i],4])),]$geo_age)>8000)
	{
    plot.which[i] <- site.id.list[i]
    }
}

plot.which.keep<-plot.which[order(lat.save)]
plot.which.keep<-plot.which.keep[plot.which.keep!=0]
plot.which.keep<- plot.which.keep[-c(1,6,11)]

site.pick <- site.id.list[plot.which.keep]

calc.all <- cast.x
calc.all1 <- calc.all[which(calc.all$site.id%in%site.pick),]


breaks <-  c(seq(0,40,10),seq(80,160,40))
colors <- rev(terrain.colors(length(breaks)-1))
data_binned <-  cut(calc.all1[,ncol(calc.all1)], c(breaks), include.lowest = FALSE, labels = FALSE)

pdf("pick_sites1.pdf")

map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(calc.all1$lon, calc.all1$lat, pch=21,
		cex=1.1, bg=colors[data_binned],lwd=.2)
title("Biomass Point Estimates at Settlement")
plotInset(-90,47,-82.5,50,
          expr={
          	hist(calc.all1[,27],col=colors,xaxt="n",xlab=NA,
          	ylab=NA,main=NA,cex.lab=.5,cex.axis=.5,breaks=breaks)
          	axis(side=1,breaks,at = breaks,cex.axis = .5,las=2,line=0)
          	mtext(side = 1, "Biomass (Mg/ha)", line = 1.5,cex=.5)
         mtext(side = 2, "Frequency", line = 1.7,cex=.5)
          	})

for(i in which(unique(x.meta[,1])%in%site.pick)){
	  par(mfrow=c(1,1))
	  map('state', xlim=c(-97.3,-83), ylim=c(41.5,50))
      points(x.meta[x.meta[,1]==unique(x.meta[,1])[i],3],      
      x.meta[x.meta[,1]==unique(x.meta[,1])[i],2], pch=19, cex=1)  
      #title(unique(x.meta[,1])[i])	
      	title(x.meta[x.meta[,1]== site.id.list[i],5][1])
      	
     fig.mat <- matrix(1,33,1)
     fig.mat[1:6,]<-1
     fig.mat[7:13,]<-2
     fig.mat[14:33,]<-seq(3,22,1)
     
     layout(fig.mat)
     
      par(mar = c(0, 4, 0, 3), oma = c(5, 0, 4, 0),tcl=.5)
      
	  ########## plot biomass time series
      plot(seq(100,9900,100),biomassCI[[i]][2,],cex=.1,ylim=c(0,200),xlim=c(-10,10000),ylab="Biomass",xlab="Years BP",main=NA,, xaxt='n')
      title(x.meta[x.meta[,1]== site.id.list[i],5][1],outer=TRUE)
      axis(3,at=seq(0,10000,1000),labels=seq(0,10000,1000),padj=1)

      ciEnvelope(seq(100,9900,100),biomassCI[[i]][1,],biomassCI[[i]][3,],col="lightblue")

      points(seq(100,9900,100),biomassCI[[i]][2,],cex=.5,pch=16)
      
      ########### plot biomass change time series
      plot(seq(150,9850,100),biomassCI[[i]][2,1:98]-biomassCI[[i]][2,2:99],cex=.1,xlim=c(-10,10000),ylab="Biomass",xlab="Years BP",main=NA,, xaxt='n')
      title(x.meta[x.meta[,1]== site.id.list[i],5][1],outer=TRUE)
     
      ciEnvelope(seq(150,9850,100),biomassCI[[i]][1,1:98]-biomassCI[[i]][1,2:99],biomassCI[[i]][3,1:98]-biomassCI[[i]][3,2:99],col="lightblue")
      abline(h=0,lty=2)

      points(seq(150,9850,100),biomassCI[[i]][2,1:98]-biomassCI[[i]][2,2:99],cex=.5,pch=16)
      
      keep.dataset.id <- unique(x.meta[x.meta[,1]==site.id.list[i],4])
      
      rug(x.meta[x.meta[,1]== site.id.list[i],]$age_bacon,lwd=2)
      rug(control.pts[which(control.pts[,1]%in%keep.dataset.id),]$geo_age,lwd=3,col="red")

      ten.count.use = ten.count[which(x.meta$SiteID==site.id.list[i]),]
      prop.use <- prop.table(as.matrix(ten.count.use),margin=1)    
     
      for(p in 1:ncol(prop.use)){
        prop.plot<- cbind(as.vector(x.meta[which(x.meta$SiteID==site.id.list[i]),]$age_bacon),as.matrix(prop.use[,p]))      	
        prop.plot<-prop.plot[order(prop.plot[,1]),]
        plot(x=prop.plot[,1],y=prop.plot[,2],type="l",xlim=c(0,10000),ylim=c(0,max(prop.use)),ylab=NA,yaxt='n', xaxt='n')
      	 #axis(2,at=signif(max(prop.use)/2),labels=signif(max(prop.use[,p])/2,digits=1))
      	 ciEnvelope(prop.plot[,1],rep(0,nrow(prop.use)),prop.plot[,2],col="darkblue")
      	 legend('topright',colnames(prop.use)[p],bty="n")
      } 

    }

dev.off()

breaks <-  c(seq(0,40,10),seq(80,160,40))
colors <- rev(terrain.colors(length(breaks)-1,alpha=.6))


pdf("time_series_by_settlement_biomass.pdf")  		
	
	  
      #title(x.meta[x.meta[,1]== site.id.list[i],5][1],outer=TRUE)
      #axis(3,at=seq(0,10000,1000),labels=seq(0,10000,1000),padj=1)

par(mfrow=c(4,1))
for(i in q1){
	plot(seq(100,9900,100),biomassCI[[i]][2,],cex=.1,
	  ylim=c(0,150),xlim=c(-10,10000),ylab="Biomass (Mg / Ha)",xlab="Years BP",main=NA)
	  
	  save.col <- cut(biomassCI[[i]][1,2], c(breaks), include.lowest = FALSE, labels = FALSE)
	
      ciEnvelope(seq(100,9900,100),biomassCI[[i]][1,],biomassCI[[i]][3,],col=colors[save.col])
      
      data_binned <-  
      points(seq(100,9900,100),biomassCI[[i]][2,],cex=.2,pch=16)
      abline(h=0)
     }
     
dev.off()



######QUESTION 1


ints <- list()
for(i in 1:7){
	ints[[i]]<-which(data_binned%in%i)[1:2]
	}

ints<-na.omit(unlist(ints))

data_binned[ints]
q1 <- which(names(biomassCI)%in%calc.all1$site.id[ints])

site.pick1 <- calc.all1$SiteID[ints]
calc.all2 <- calc.all[which(calc.all$site.id%in%site.pick1),]
#calc.all2<-calc.all2[-5,]
data_binned1 <-  cut(calc.all2[,ncol(calc.all2)], c(breaks), include.lowest = FALSE, labels = FALSE)

#q1<-q1[order(calc.all2[,ncol(calc.all2)])]
dbs<-numeric(500)
dbs[q1] <- data_binned1[order(calc.all2[,27])]
   
     fig.mat <- matrix(1,60,1)
     fig.mat[1:5,1]<-1
     fig.mat[6:10,1]<-2
     fig.mat[11:15,1]<-3
     fig.mat[16:20,1]<-4
     fig.mat[21:25,1]<-5
     fig.mat[26:30,1]<-6
     fig.mat[31:35,1]<-7
     fig.mat[36:40,1]<-8
     fig.mat[41:45,1]<-9
     
     fig.mat[46:50,1]<-10
     fig.mat[51:55,1]<-11
     fig.mat[56:60,1]<-12
   
pdf("map_sites_draft.pdf")
map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(calc.all2$LongitudeWest, calc.all2$LatitudeNorth, pch=21,
		cex=2, bg=colors[data_binned1],lwd=.2)
text(calc.all2$LongitudeWest, calc.all2$LatitudeNorth,labels=data_binned1)
dev.off()
     
pdf("Differences_by_settlement_biomass.pdf")

map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(calc.all2$LongitudeWest, calc.all2$LatitudeNorth, pch=21,
		cex=2, bg=colors[data_binned1],lwd=.2)
text(calc.all2$LongitudeWest, calc.all2$LatitudeNorth,labels=data_binned1)

     layout(fig.mat)
     
     par(mar = c(0, 4, 0, 3), oma = c(5, 0, 4, 0),tcl=.5)   
    
for(i in q1){	
	
	 # plot(seq(100,9900,100),biomassCI[[i]][2,],cex=.1,ylim=c(0,200),xlim=c(-10,10000),ylab="Biomass",xlab="Years BP",main=NA,, xaxt='n')
      # #title(x.meta[x.meta[,1]== site.id.list[i],5][1],outer=TRUE)
      # axis(3,at=seq(0,10000,1000),labels=seq(0,10000,1000),padj=1)

      # ciEnvelope(seq(100,9900,100),biomassCI[[i]][1,],biomassCI[[i]][3,],col=colors[dbs[i]])

      # points(seq(100,9900,100),biomassCI[[i]][2,],cex=.5,pch=16)
color.pts <- numeric(98)
pts <- numeric(98); pts.025 <- numeric(98) ; pts.975 <- numeric(98)
for(p in 1:98){
	norm.pts <- rnorm(100,as.vector(biomassCI[[i]][2,1:98]-biomassCI[[i]][2,2:99])[p],max(as.vector(biomassCI[[i]][2,1:98]-biomassCI[[i]][2,2:99])/10+rnorm(1,.1,.01))) #HUGE HACK
	pts[p]<-mean(norm.pts)
	pts.025[p]<-quantile(norm.pts,.025)
	pts.975[p]<-quantile(norm.pts,.975)
}

for(r in 1:length(pts)){
	if(pts[r]<=0) {
		color.pts[r]<-"blue"
	}else{
		color.pts[r]<-"red"
	}
}

	plot(seq(150,9850,100),pts,cex=1,pch=19,xlim=c(-10,10000),ylab="Biomass",xlab="Years BP",main=NA,xaxt='n',col=color.pts)
      ciEnvelope(seq(150,9850,100),pts.025,pts.975,col=colors[dbs[i]])
      abline(h=0,lty=2)     
      points(seq(150,9850,100),pts,cex=1,pch=19,col=color.pts)
      #abline(v=7200)
      abline(v=c(7200,3000))
      
}

axis(1,at=seq(0,10000,1000),labels=seq(0,10000,1000),padj=1)
title("Differences by Settlement Biomass",outer=TRUE)  
dev.off()


####Question 3
pdf("relative_differences_2.pdf")
map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(calc.all2$LongitudeWest, calc.all2$LatitudeNorth, pch=21,
		cex=2, bg=colors[data_binned1],lwd=.2)
text(calc.all2$LongitudeWest, calc.all2$LatitudeNorth,labels=data_binned1)


 layout(fig.mat)
     
 par(mar = c(0, 4, 0, 3), oma = c(5, 0, 4, 0),tcl=.5)   
    
for(i in q1){
	pts <- numeric(98); pts.025 <- numeric(98) ; pts.975 <- numeric(98)
	norm.pts <- as.vector(biomassCI[[i]][2,seq(1,98,5)]/biomassCI[[i]][2,seq(2,99,5)])
    hist(norm.pts,col=colors[dbs[i]],xlim=c(.85,1.125),ylim=c(0,100),main=NA,breaks=20,freq=FALSE)
    #text(.875, 50,labels=paste0("Biomass Cat.",dbs[i]))

pts.025<-quantile(norm.pts,.025)
pts.975<-quantile(norm.pts,.975)
abline(v=c(pts.025,pts.975),lty=2,lwd=3)
abline(v=1,lty=1,col="red",lwd=3)


}
axis(1,at=seq(0,10000,1000),labels=seq(0,10000,1000),padj=1)
title("Absolute Biomass Change",outer=TRUE) 
dev.off()

######Question 2

cor.calc.mat=rep(NA,21)

pdf('biomass_corrs_pick_plot.pdf')

for(i in q1){
	  par(mfrow=c(1,1))
	  map('state', xlim=c(-97.3,-83), ylim=c(41.5,50))
      points(x.meta[x.meta[,1]==unique(x.meta[,1])[i],3],      
      x.meta[x.meta[,1]==unique(x.meta[,1])[i],2], pch=19, cex=1)  
      #title(unique(x.meta[,1])[i])	
      	title(x.meta[x.meta[,1]== site.id.list[i],5][1])
      	
     fig.mat <- matrix(1,33,1)
     fig.mat[1:6,]<-1
     fig.mat[7:13,]<-2
     fig.mat[14:33,]<-seq(3,22,1)
     
     layout(fig.mat)
     
      par(mar = c(0, 4, 0, 3), oma = c(5, 0, 4, 0),tcl=.5)
      
	  ########## plot biomass time series
      plot(seq(100,9900,100),biomassCI[[i]][2,],cex=.1,ylim=c(0,200),xlim=c(-10,10000),ylab="Biomass",xlab="Years BP",main=NA,, xaxt='n')
      title(x.meta[x.meta[,1]== site.id.list[i],5][1],outer=TRUE)
      axis(3,at=seq(0,10000,1000),labels=seq(0,10000,1000),padj=1)

      ciEnvelope(seq(100,9900,100),biomassCI[[i]][1,],biomassCI[[i]][3,],col=colors[dbs[i]])

      points(seq(100,9900,100),biomassCI[[i]][2,],cex=.5,pch=16)
      
      ########### plot biomass change time series
      plot(seq(150,9850,100),biomassCI[[i]][2,1:98]-biomassCI[[i]][2,2:99],cex=.1,xlim=c(-10,10000),ylab="Biomass",xlab="Years BP",main=NA,, xaxt='n')
      title(x.meta[x.meta[,1]== site.id.list[i],5][1],outer=TRUE)
     
      ciEnvelope(seq(150,9850,100),biomassCI[[i]][1,1:98]-biomassCI[[i]][1,2:99],biomassCI[[i]][3,1:98]-biomassCI[[i]][3,2:99],col=colors[dbs[i]])
      abline(h=0,lty=2)

      points(seq(150,9850,100),biomassCI[[i]][2,1:98]-biomassCI[[i]][2,2:99],cex=.5,pch=16)
      
      keep.dataset.id <- unique(x.meta[x.meta[,1]==site.id.list[i],4])
      
      rug(x.meta[x.meta[,1]== site.id.list[i],]$age_bacon,lwd=2)
      rug(control.pts[which(control.pts[,1]%in%keep.dataset.id),]$geo_age,lwd=3,col="red")

      ten.count.use = ten.count[which(x.meta$SiteID==site.id.list[i]),]
      prop.use <- prop.table(as.matrix(ten.count.use),margin=1)    
     
      for(p in 1:ncol(prop.use)){
        prop.plot<- cbind(as.vector(x.meta[which(x.meta$SiteID==site.id.list[i]),]$age_bacon),as.matrix(prop.use[,p]))      	
        prop.plot<-prop.plot[order(prop.plot[,1]),]
        plot(x=prop.plot[,1],y=prop.plot[,2],type="l",xlim=c(0,10000),ylim=c(0,max(prop.use)),ylab=NA,yaxt='n', xaxt='n')
      	 #axis(2,at=signif(max(prop.use)/2),labels=signif(max(prop.use[,p])/2,digits=1))
      	 ciEnvelope(prop.plot[,1],rep(0,nrow(prop.use)),prop.plot[,2],col="darkblue")
      	 legend('topright',colnames(prop.use)[p],bty="n")
      } 
      round.keep <- round(prop.plot[,1]/100,digits=0)
      if(max(round.keep)>97){
      	prop.use<-prop.use[-which(round.keep>97),]
        round.keep<-round.keep[-which(round.keep>97)]
      }
      
      
      Biomass <- as.vector(biomassCI[[i]][2,round.keep])
      cor.calc <- cor(x=cbind(prop.use,Biomass),y=cbind(prop.use,Biomass))
      cor.calc.mat<-cbind(cor.calc.mat,cor.calc[,21])
      
      cor.calc[is.na(cor.calc)]<-0
        par(mfrow=c(1,1))
      corrplot(cor.calc,type="upper",order="FPC")

    }
    
cor.calc.mat<-cor.calc.mat[,-1]
colnames(cor.calc.mat)<-name.list[q1]
cor.calc.mat[is.na(cor.calc.mat)]<-0
corrplot(cor.calc.mat)
    
dev.off()


###########
########### Ecosystem Stability ###########
###########

pdf('log_v_reg_biomass.pdf')
par(mfrow=c(3,3))
for(i in 1:length(biomassCI)){
	if(length(biomassCI[[i]])>10){
		map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(x.meta[x.meta[,1]==unique(x.meta[,1])[i],3],      
      x.meta[x.meta[,1]==unique(x.meta[,1])[i],2], pch=19, cex=1) 
		plot(log(biomassCI[[i]][2,]),ylim=c(0,10),main='log biomass')
		plot(biomassCI[[i]][2,],ylim=c(0,200),main = 'biomass')
	}
} 
dev.off()

h=5
t=6
second.deriv<-list()
##different time steps
##all possible steps
calc.second.deriv <- function(biomassCI,h,second.deriv){
	for(i in 1:length(biomassCI)){
	if(length(biomassCI[[i]])>10){
		T <- dim(biomassCI[[i]])[2] - h
		t <- h + 1 
		biomassCI[[i]]<-log(biomassCI[[i]])
		second.deriv[[i]]<-sum(((biomassCI[[i]][2,(t:T)+h]-2*biomassCI[[i]][2,(t:T)]+biomassCI[[i]][2,(t:T)-h])/((h*100)^2))^2)
		
	}else{
		second.deriv[[i]]<-NA
	}
}
return(second.deriv)
}
pdf(paste('second.deriv.map.subset.log',Sys.Date(),'.pdf'))
for(i in c(1,2,5,10,20)){
h=i
second.deriv<-calc.second.deriv(biomassCI=biomassCI,h=h,second.deriv=list())

#hist(unlist(second.deriv))
names(second.deriv) <- unique(x.meta[,1])[1:182]
second.deriv.unlist <- unlist(second.deriv)
second.deriv.unlist <- second.deriv.unlist[plot.which.keep]
#rownames(calc.all1)<-seq(1,34,1)
second.deriv.copy <- second.deriv
second.deriv.unlist.names <- unlist(second.deriv.copy)
second.deriv.unlist.names <- second.deriv.unlist.names[plot.which.keep]

for(q in 1:length(second.deriv.unlist.names)){
  names(second.deriv.unlist.names)[q] <- as.character(x.meta[x.meta$site.id==names(second.deriv.unlist[q]),'site.name'][1])
}

breaks <-  c(0,quantile(probs=c(.05,.25,.5,.75,.975),second.deriv.unlist,na.rm=TRUE),max(second.deriv.unlist)+1)
 
  long.keep <- list()
  lat.keep <- list()
  for(i in 1:length(second.deriv.unlist)){
  	long.keep[[i]] <- x.meta[x.meta[,1]==names(second.deriv.unlist)[i],'long'][1]
  	lat.keep[[i]] <- x.meta[x.meta[,1]==names(second.deriv.unlist)[i],'lat'][1]
  }

colors <- colorRampPalette(c("white","yellow",'green','blue'))(length(breaks)-1)

data_binned <-  cut(second.deriv.unlist, c(breaks), include.lowest = FALSE, labels = FALSE)

map('state', xlim=c(-98,-81), ylim=c(41.5,50))
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



first.deriv<-list()
first.deriv.calc <- function(biomassCI,first.deriv,h){
	for(i in 1:183){
	if(length(biomassCI[[i]])>1){
first.deriv[[i]]<-sum((biomassCI[[i]][2,rev(seq(2,99,h))]-biomassCI[[i]][2,rev(seq(1,98,h))])/(rev(seq(200,9900,h*100))-rev(seq(100,9800,h*100))))
	}else{
		first.deriv[[i]]<-NA
	}
    }  
    return(first.deriv)
}

h=10
first.deriv <- first.deriv.calc(biomassCI=biomassCI,first.deriv=list(),h=h)

#hist(unlist(first.deriv))
names(first.deriv) <- unique(x.meta[,1])
first.deriv.unlist <- unlist(first.deriv)
first.deriv.unlist.comp <- first.deriv.unlist[which(names(first.deriv.unlist)%in%calc.all1$SiteID)]
rownames(calc.all1)<-seq(1,34,1)
calc.all3<-cbind(calc.all1[-c(14,17),],first.deriv.unlist.comp)

breaks <-  c(min(first.deriv.unlist.comp)-1,-.5,-.05,0,.05,.5,max(first.deriv.unlist.comp)+1)
colors <- colorRampPalette(c("blue",'white', "red"))(length(breaks)-1)
#rev(rainbow(length(breaks)-1,start=.2,end=0))#want different colors

data_binned <-  cut(calc.all3[,30], c(breaks), include.lowest = FALSE, labels = FALSE)

pdf(paste('first.deriv.map',h*100,'.pdf'))
map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(calc.all3$LongitudeWest, calc.all3$LatitudeNorth, pch=21,
		cex=1.3, bg=colors[data_binned],lwd=.2)
title(paste("First Derivative Sum Point Estimates h = ",h*100))

plotInset(-90,47,-82.5,50,
          expr={
          	hist(data_binned,col=colors,xaxt="n",xlab=NA,
          	ylab=NA,main=NA,cex.lab=.5,cex.axis=.5,breaks=0:7)
          	axis(side=1,signif(breaks[unique(data_binned)],digits=2),at = sort(unique(data_binned)),cex.axis = .5,las=2,line=0)
          	mtext(side = 1, "Squared Sum", line = 1.5,cex=.5)
            mtext(side = 2, "Frequency", line = 1.7,cex=.5)
          	})
dev.off()

par(mfrow=c(4,1))

pdf('derivative.timeseries2.pdf')
for(i in plot.which.keep[c(6,8,10,13)]){
par(mfrow=c(4,1),mar=c(2,4,2,2))
plot(seq(100,9900,100),biomassCI[[i]][2,],pch=1,cex=.1,ylab='Biomass',xlab="Years BP",ylim=c(0,350))
ciEnvelope(seq(100,9900,100),biomassCI[[i]][1,],biomassCI[[i]][3,],col="lightblue")
points(seq(100,9900,100),biomassCI[[i]][2,],pch=19,cex=.5)

plot(seq(100,9800,100),diff(biomassCI[[i]][2,])/diff(seq(100,9900,100)),typ='l',ylab='First Derivative',xlab="Years BP")
abline(h=0,col='red')
mtext(paste("sum = ",signif(sum(diff(biomassCI[[i]][2,])/diff(seq(100,9900,100))),digits=2)),side=3)

second.deriv.plot <- ((biomassCI[[i]][2,(t:T)+h]-2*biomassCI[[i]][2,(t:T)]+biomassCI[[i]][2,(t:T)-h])/(h*100^2))^2
plot(seq(100,9900,length.out=length(second.deriv.plot)), second.deriv.plot
,typ='l',,ylab='Second Derivative',xlab="Years BP")
abline(h=0,col='red')
mtext(paste("sum = ",signif(sum(((biomassCI[[i]][2,(t:T)+h]-2*biomassCI[[i]][2,(t:T)]+biomassCI[[i]][2,(t:T)-h])/(h^2))^2),digits=2)),side=3)

map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(x.meta[x.meta[,1]==unique(x.meta[,1])[i],3],      
      x.meta[x.meta[,1]==unique(x.meta[,1])[i],2], pch=19, cex=1) 

fig.mat <- matrix(1,33,1)
     fig.mat[1:6,]<-1
     fig.mat[7:13,]<-2
     fig.mat[14:33,]<-seq(3,22,1)
     
     layout(fig.mat)
     
      par(mar = c(0, 4, 0, 3), oma = c(5, 0, 4, 0),tcl=.5)
      
	  ########## plot biomass time series
      plot(seq(100,9900,100),biomassCI[[i]][2,],cex=.1,ylim=c(0,350),xlim=c(-10,10000),ylab="Biomass",xlab="Years BP",main=NA,, xaxt='n')
      title(x.meta[x.meta[,1]== site.id.list[i],5][1],outer=TRUE)
      axis(3,at=seq(0,10000,1000),labels=seq(0,10000,1000),padj=1)

      ciEnvelope(seq(100,9900,100),biomassCI[[i]][1,],biomassCI[[i]][3,],col="lightblue")

      points(seq(100,9900,100),biomassCI[[i]][2,],cex=.5,pch=16)
      
      ########### plot biomass change time series
       plot(seq(150,9850,100),biomassCI[[i]][2,1:98]-biomassCI[[i]][2,2:99],cex=.1,xlim=c(-10,10000),ylab="Biomass",xlab="Years BP",main=NA,, xaxt='n')
      title(x.meta[x.meta[,1]== site.id.list[i],5][1],outer=TRUE)
     
      ciEnvelope(seq(150,9850,100),biomassCI[[i]][1,1:98]-biomassCI[[i]][1,2:99],biomassCI[[i]][3,1:98]-biomassCI[[i]][3,2:99],col="lightblue")
      abline(h=0,lty=2)

      points(seq(150,9850,100),biomassCI[[i]][2,1:98]-biomassCI[[i]][2,2:99],cex=.5,pch=16)
      
      keep.dataset.id <- unique(x.meta[x.meta[,1]==site.id.list[i],4])
      
      rug(x.meta[x.meta[,1]== site.id.list[i],]$age_bacon,lwd=2)
      rug(control.pts[which(control.pts[,1]%in%keep.dataset.id),]$geo_age,lwd=3,col="red")

      ten.count.use = ten.count[which(x.meta$SiteID==site.id.list[i]),]
      prop.use <- prop.table(as.matrix(ten.count.use),margin=1)    
     
      for(p in 1:ncol(prop.use)){
        prop.plot<- cbind(as.vector(x.meta[which(x.meta$SiteID==site.id.list[i]),]$age_bacon),as.matrix(prop.use[,p]))      	
        prop.plot<-prop.plot[order(prop.plot[,1]),]
        plot(x=prop.plot[,1],y=prop.plot[,2],type="l",xlim=c(0,10000),ylim=c(0,max(prop.use)),ylab=NA,yaxt='n', xaxt='n')
      	 #axis(2,at=signif(max(prop.use)/2),labels=signif(max(prop.use[,p])/2,digits=1))
      	 ciEnvelope(prop.plot[,1],rep(0,nrow(prop.use)),prop.plot[,2],col="darkblue")
      	 legend('topright',colnames(prop.use)[p],bty="n")
      } 

      
	
}
dev.off()


