setwd("/Users/paleolab/babySTEPPS/")
data.dir = c("/Users/paleolab/babySTEPPS/Data/")
fig.dir = c('/Users/paleolab/babySTEPPS/Figures/')

#####
##### Install Packages #####
#####
library(reshape)
library(ggplot2)
library(sp)
library(rgdal)
library(fields)
library(maptools)
library(neotoma)
require(grid)
require(plyr)
require(maps)
library(gridExtra)
gpclibPermit()
library(mgcv)
library(splines)
library(boot)
library(gtools)
library(rjags)
library(oce)
library(mgcv)
library(nimble)

usShp <- readShapeLines(file.path(data.dir, '/us_alb.shp'), proj4string=CRS('+init=epsg:3175'))
usShp@data$id <- rownames(usShp@data)
usFortified <- fortify(usShp, region='id')
###### MNWI_dat.R #####

#####
##### Download Data for MN and WI and MI #####
#####
gpids <- get_table(table.name='GeoPoliticalUnits')
gpid <- gpids[which(gpids$GeoPoliticalName == 'Minnesota'),1]
gpid1 <- gpids[which(gpids$GeoPoliticalName == 'Wisconsin'),1]
gpid2 <- gpids[which(gpids$GeoPoliticalName == 'Michigan'),1]
meta <- get_dataset(datasettype='pollen', gpid=c(gpid), ageyoung=0) #Pollen data for all of MN&WI
meta1 <- get_dataset(datasettype='pollen', gpid=c(gpid1), ageyoung=0) #Pollen data for all of MN&WI
meta2 <- get_dataset(datasettype='pollen', gpid=c(gpid2), ageyoung=0) #Pollen data for all of MN&WI

site.locs <- ldply(meta, function(x) c(x$site.data$long, x$site.data$lat))
site.locs1 <- ldply(meta1, function(x) c(x$site.data$long, x$site.data$lat))
site.locs2 <- ldply(meta2, function(x) c(x$site.data$long, x$site.data$lat))
site.locs<-rbind(site.locs,site.locs1,site.locs2)

mnwi1 <- rep(0,length(c(meta)))
mnwi2 <- rep(0,length(c(meta1)))
mnwi3 <- rep(0,length(c(meta2)))
for(i in 1:length(meta)) mnwi1[i] <- meta[[i]]$dataset.meta$dataset.id
for(i in 1:length(meta1)) mnwi2[i] <- meta1[[i]]$dataset.meta$dataset.id
for(i in 1:length(meta2)) mnwi3[i] <- meta2[[i]]$dataset.meta$dataset.id
datasets <- as.vector(c(mnwi1,mnwi2,mnwi3))
dat.mnwi <- get_download(x = as.vector(c(mnwi1,mnwi2,mnwi3)))
save(dat.mnwi,file="mnwi4.rdata")

  comp.tax <- compile_taxa(dat.mnwi[[1]], 'WhitmoreSmall')
	
	samp.meta <- comp.tax$sample.meta
	counts <- comp.tax$counts
	dat.meta <- matrix(unlist(comp.tax$dataset$site.data),nrow(counts),8,byrow=TRUE)
	colnames(dat.meta) <- names(comp.tax$dataset$site.data)
	pol.cal.count<-cbind(dat.meta,samp.meta,counts)
		
for(i in 2:length(dat.mnwi)){
	comp.tax <- compile_taxa(dat.mnwi[[i]], 'WhitmoreSmall')
	samp.meta <- comp.tax$sample.meta
	dat.meta <- matrix(unlist(comp.tax$dataset$site.data),nrow(samp.meta),8,byrow=TRUE)
	colnames(dat.meta) <- names(comp.tax$dataset$site.data)
	counts <- comp.tax$counts
	cbound<-cbind(dat.meta,samp.meta,counts)
	pol.cal.count <- smartbind(pol.cal.count,cbound)
}

pol.cal.count <- data.frame(pol.cal.count, stringsAsFactors=FALSE)
pol.cal.count$lat <- as.numeric(as.character(pol.cal.count$lat))
pol.cal.count$long <- as.numeric(as.character(pol.cal.count$long))

##### Fix problems with matrix
rownames(pol.cal.count)<-seq(1,nrow(pol.cal.count),1)
pol.cal.count[is.na(pol.cal.count)]<-0
#save(pol.cal.count,file="pol.cal.count.mnwi2.csv")

#####
##### Add Bacon Pollen Dates #####
#####
bacon<-read.csv(paste0(data.dir,'/pollen_ts_bacon_v6.csv'))

### makes a dataframe with only sites that have bacon ages
new.pol1<-cbind(pol.cal.count[pol.cal.count$dataset.id==unique(bacon$id[1]),],bacon[bacon$id==unique(bacon$id[1]),])

for(i in 2:26){
	if(nrow(pol.cal.count[pol.cal.count$dataset.id==unique(bacon$id)[i],])==nrow(bacon[bacon$id==unique(bacon$id)[i],])){
		new.pol<-cbind(pol.cal.count[pol.cal.count$dataset.id==unique(bacon$id)[i],],
	    bacon[bacon$id==unique(bacon$id)[i],])
	    new.pol1<-rbind(new.pol1,new.pol)
	}else{
		my.pol <- pol.cal.count[pol.cal.count$dataset.id==unique(bacon$id)[i],]
		bacon.pol <- bacon[bacon$id==unique(bacon$id)[i],]
		new.pol2 <- which(my.pol$age%in%bacon.pol$age_default)
		new.pol<-cbind(my.pol[new.pol2,],bacon.pol[which(bacon.pol$age_default%in%my.pol$age),])
	    new.pol1<-rbind(new.pol1,new.pol)
	}
}

for(i in 27:length(unique(bacon$id))){
	if(nrow(pol.cal.count[pol.cal.count$dataset.id==unique(bacon$id)[i],])==nrow(bacon[bacon$id==unique(bacon$id)[i],])){
		new.pol<-cbind(pol.cal.count[pol.cal.count$dataset.id==unique(bacon$id)[i],],
	    bacon[bacon$id==unique(bacon$id)[i],])
	    new.pol1<-rbind(new.pol1,new.pol)
	}else{
		my.pol <- pol.cal.count[pol.cal.count$dataset.id==unique(bacon$id)[i],]
		bacon.pol <- bacon[bacon$id==unique(bacon$id)[i],]
		new.pol2 <- which(my.pol$age%in%bacon.pol$age_default)
		new.pol<-cbind(my.pol[new.pol2,],bacon.pol[which(bacon.pol$age_default%in%my.pol$age),])
	    new.pol1<-rbind(new.pol1,new.pol)
	    #print(unique(bacon$id)[i])
	}
}

###Tests

new.pol1$PINE - new.pol1$PINUSX #should equal zero
new.pol1$HEMLOCK - new.pol1$TSUGA

new.pol2 <- new.pol1[,c('site.id','lat','long','dataset.id','site.name','age_bacon',colnames(new.pol1)[20:99])]

#####
###### Format settlement pollen data #####
#####
x = new.pol2[new.pol2$age_bacon>=50,]
x = x[x$age_bacon<=150,]

melt.x <- melt(x,id.vars=colnames(new.pol2)[1:6])
cast.x <- cast(melt.x,lat + long + site.id + site.name ~ variable,sum)
cast.x=as.data.frame(cast.x)
#cast.x = as.data.frame(cast.x[-c(1:3),])
#row_keep = rep(0,nrow(cast.x))
#plot_biomass_pollen = matrix(0,nrow(cast.x),(ncol(x)-3))

#####
##### Adding settlement biomass data #####
#####

biomass_dat_est <- read.csv(paste0(data.dir,"biomass_prediction_v0.9-7_bam.csv"))

##### Changing pollen coordinates so that we can find the right rows when we find the biomass for each pond.
lat.long.reg <- cbind(as.numeric(as.character(cast.x$long)),as.numeric(as.character(cast.x$lat)))
lat.long.reg.df = data.frame(lat.long.reg)
colnames(lat.long.reg.df) = c('x', 'y')

coordinates(lat.long.reg.df) <- ~ x + y
proj4string(lat.long.reg.df) <- CRS('+proj=longlat +ellps=WGS84')

albers <- spTransform(lat.long.reg.df, CRS('+init=epsg:3175'))
albers <- as.matrix(data.frame(albers))

centers_biomass = cbind(biomass_dat_est$x,biomass_dat_est$y)
idx_cores = vector(length=nrow(cast.x))

for(i in 1:nrow(cast.x)){   
  core_site = albers[i,]
  d = rdist(matrix(core_site, ncol=2), as.matrix(centers_biomass))
  if(min(d)<8000){
  	idx_cores[i] = which.min(d) 
  }else{
  	idx_cores[i] = NA 
  }
  
}

if(DRAW == TRUE) pdf(paste0(dump.dir,"check_points.pdf")) else (quartz())
plot(albers[,1], albers[,2])
points(centers_biomass[idx_cores,1], centers_biomass[idx_cores,2], col='blue', pch=8)
plot(usShp, add=T, lwd=2)
if(DRAW == TRUE) dev.off()

biomass <- list()
for(i in 1:nrow(cast.x)){ 
  biomass[[i]] = biomass_dat_est[idx_cores[i],'Total']
 }
 
cast.x <- cbind(cast.x,unlist(biomass))
cast.x <- cast.x[-which(is.na(cast.x[,ncol(cast.x)])),]

breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)-1))
data_binned <-  cut(cast.x[,ncol(cast.x)], c(breaks), include.lowest = FALSE, labels = FALSE)

pdf(paste0(fig.dir,paste0("biomass.pts.settlement",Sys.Date(),".pdf")))
map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(cast.x$long,cast.x$lat, pch=21,
		cex=1.1, bg=colors[data_binned],lwd=.2)
title("Biomass Point Estimates at Settlement")
plotInset(-90,47,-82.5,50,
          expr={
          	hist(data_binned,col=colors,xaxt="n",xlab=NA,
          	ylab=NA,main=NA,cex.lab=.5,cex.axis=.5)
          	axis(side=1,breaks,at = seq(1,12,1),cex.axis = .5,las=2,line=0)
          	mtext(side = 1, "Biomass (Mg/ha)", line = 1.5,cex=.5)
         mtext(side = 2, "Frequency", line = 1.7,cex=.5)
          	})
dev.off()

#####
##### Creating a calibration dataset with the species we want to use
#####

###with pine

trees <- c('PINUSX',"ALNUSX","JUGLANSX","ACERX","CUPRESSA","FRAXINUX","FAGUS","CYPERACE","LARIXPSEU","TSUGAX","QUERCUS","TILIA","BETULA","PICEAX","OSTRYCAR","ULMUS","ABIES","POPULUS")
other.trees <- c("TAXUS","NYSSA","CASTANEA","PLATANUS","SALIX","LIQUIDAM")
ten.count = matrix(0,nrow(cast.x),length(trees)+3)
prairie <- c("ARTEMISIA","ASTERX","POACEAE","AMBROSIA","CHENOAMX","CORYLUS")
ten.count[,1] <- unlist(rowSums(cast.x[,prairie]))
ten.count[,2] <- unlist(rowSums(cast.x[,other.trees]))
ten.count[,3:(length(trees)+2)] <- as.matrix(cast.x[,trees])
ten.count[,(length(trees)+3)] <- rowSums(cast.x[,5:84]) - rowSums(ten.count)
colnames(ten.count)<-c("prairie","other trees",trees,"other herbs")

#####
##### Remove sites and take subset of points
#####

if(DRAW == TRUE) pdf(paste0(fig.dir,paste0("all_sites",Sys.Date(),".pdf")))
map('state', ylim=range(cast.x$lat)+c(-2, 2), xlim=range(cast.x$long)+c(-1, 1),main=NA)
points(cast.x$long, cast.x$lat, pch=19, cex=1,col="gray")
title(main="all sites")

set.seed(4)
sites_rm = sample(1:nrow(cast.x),round(nrow(cast.x)/3))

points(cast.x$long[-sites_rm], cast.x$lat[-sites_rm], pch=19, cex=1,col="blue")
if(DRAW == TRUE) dev.off()

biomass = cast.x[,ncol(cast.x)]
total_counts = round(rowSums(ten.count))
counts = round(ten.count)
props = counts/rowSums(counts)

total_counts_spp = colSums(counts)

props = props[,order(total_counts_spp,decreasing=TRUE)]

pdf(paste0(fig.dir,"scatter.newdata",Sys.Date(),".pdf"))
#quartz()
par(mfrow=c(4,4))
for(i in 1:ncol(props)){
  if(length(unique(props[-sites_rm,i]))>=9){
    plot(biomass[-sites_rm],props[-sites_rm,i],main=colnames(props)[i],xlab="biomass",ylab="pollen prop",pch = 19, cex = .5)
  }
}  
dev.off()

props = as.data.frame(props)

#plot_biomass_pollen = plot_biomass_pollen[-c(which(props$Other>.5),which(props$POACEAE>.8)),]
#####
##### Drawing Splines (with all data) #####
#####

Z = bs(biomass,intercept=TRUE,df=5)
betas = matrix(0,ncol(Z),ncol(counts)); betas.save = betas

#if(DRAW == TRUE) pdf(paste0(dump.dir,"splines.new.pdf"))
quartz()
par(mfrow=c(3,3))
for(i in 1:ncol(counts)){
  gam_mod = gam(cbind(counts[,i],total_counts-counts[,i]) ~ s(biomass),family=binomial(link="logit"))
  plot(biomass,counts[,i]/total_counts,pch=19,cex=.4,col='grey',ylab="Pollen Prop",main=colnames(counts)[i])
  points(biomass,predict(gam_mod,type="response"),pch=19,col="green")
  
  glm_mod = glm(cbind(counts[,i],total_counts-counts[,i]) ~ Z - 1,family=binomial(link="logit"))   
  points(biomass,counts[,i]/total_counts,pch=19,cex=.4,col='grey')
  new.biomass = seq(1,200,1)
  Z.new = bs(new.biomass,intercept=TRUE,df = ncol(Z))
  lines(new.biomass, predict(glm_mod,newdata=list(Z=Z.new),type="response"),col="blue")  
  
  betas[,i] = glm_mod$coefficients
  betas.save[,i] = glm_mod$coefficients
  
}

Z.knots<- Z

u<-c(rep(attr(Z.knots,"Boundary.knots")[1],1),attr(Z.knots,"knots"),rep(attr(Z.knots,"Boundary.knots")[2],1))

source('~/babystepps/Code/bs_nimble.R')

Z.knots.check = matrix(0,nrow=u[length(u)],ncol=(length(u)+2));

for(i in 1:u[length(u)]){
	u_given <-i
	Z.knots.check[i,] = bs_nimble(u_given, u=u, N0 = rep(0, (length(u)-1)),N1 = rep(0, (length(u))), N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))
}

for(i in 1:length(biomass)){
    u_given <- biomass[i]
	Z.knots[i,] = bs_nimble(u_given, u=u, N0 = rep(0, (length(u)-1)),N1 = rep(0, (length(u))), N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))
}

#if(DRAW == TRUE) dev.off()

#####
##### Create final calibration datasets #####
#####

Y = counts[-sites_rm,] #remove sites "sites_rm" defined above
Y <- Y[,rev(order(colMeans(Y)))]
biomass = biomass[-sites_rm]
counts = counts[-sites_rm,]
counts <- Y[,rev(order(colMeans(Y)))]
total_counts = rowSums(counts)

save(Y,biomass,file='calibration.data.Rdata')

#save.image(file="add.bacon2.Rdata")

u<-c(rep(attr(Z,"Boundary.knots")[1],1),attr(Z,"knots"),rep(attr(Z,"Boundary.knots")[2],1))

x = new.pol1[new.pol1$age_bacon>=200,]
x = x[x$age_bacon<=10000,]

x.meta = x[,c('site.id','lat',"long","dataset.id","site.name","age_bacon")]

trees <- c("PINUSX","ALNUSX","JUGLANSX","ACERX","CUPRESSA","FRAXINUX","FAGUS","CYPERACE","LARIXPSEU","TSUGAX","QUERCUS","TILIA","BETULA","PICEAX","OSTRYCAR","ULMUS","ABIES","POPULUS")
other.trees <- c("TAXUS","NYSSA","CASTANEA","PLATANUS","SALIX","LIQUIDAM")
ten.count = matrix(0,nrow(x),length(trees)+3)
prairie <- c("ARTEMISIA","ASTERX","POACEAE","AMBROSIA","CHENOAMX","CORYLUS")
ten.count[,1] <- unlist(rowSums(x[,prairie]))
ten.count[,2] <- unlist(rowSums(x[,other.trees]))
ten.count[,3:(length(trees)+2)] <- as.matrix(x[,trees])
ten.count[,(length(trees)+3)] <- rowSums(x[,20:99]) - rowSums(ten.count)
colnames(ten.count)<-c("prairie","other trees",trees,"other herbs")

ten.count.save = ten.count
ten.count = round(ten.count.save)

ten.count <- ten.count[,colnames(counts)]

save(x.meta,ten.count,file = 'prediction.data.Rdata')

#####
##### Plots #####
#####

pdf(paste0(fig.dir,"all.sites.neotoma",Sys.Date(),".pdf"))
map('state', xlim=range(site.locs$V1)+c(-2, 2), ylim=range(site.locs$V2)+c(-1, 1))
points(site.locs$V1, site.locs$V2, pch=19, cex=1,col="black")
title("All Pollen Sites")
dev.off()

pdf(paste0(fig.dir,"chrono.10k",Sys.Date(),".pdf"))
par(mfrow=c(1,1))
plot(ponds1$age,ponds1$lat,pch=19,cex=.5,xlim=c(0,10000),main="All Pollen Records",ylab="Latitude",xlab="Age BP")
dev.off()


plot.seq = seq(0,20000,500)

par(mfrow = c(2,2))
for(i in 2:length(plot.seq)){
  map('state', xlim=range(as.numeric(as.character(ponds1$long)))+c(-2, 2), ylim=range(as.numeric(as.character(ponds1$lat)))+c(-1, 1))
  points(ponds1[ponds1$age>plot.seq[i-1]&ponds1$age<plot.seq[i],]$long, 
         ponds1[ponds1$age>plot.seq[i-1]&ponds1$age<plot.seq[i],]$lat, 
         pch=19, cex=.5)
  title(c(plot.seq[i-1],"-",plot.seq[i]))
}

#### Plot basis functions ####
plot(Z.knots.check[,1],xlim=c(0,u[length(u)]),pch=19,ylim=c(0,1),xlab="Biomass")
for(i in 2:ncol(Z.knots.check)){
  points(Z.knots.check[,i],col=i,pch=19)
}
abline(v=u,lwd=2)
title("Basis Functions")

