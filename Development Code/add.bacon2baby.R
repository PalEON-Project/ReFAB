setwd("/Users/paleolab/babySTEPPS/")
data.dir = c("~/babySTEPPS/Data/")
fig.dir = c("~/babySTEPPS/Figures/")


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

usShp <- readShapeLines(file.path(data.dir, '/us_alb.shp'), proj4string=CRS('+init=epsg:3175'))
usShp@data$id <- rownames(usShp@data)
usFortified <- fortify(usShp, region='id')

bacon<-read.csv(paste0(data.dir,'/pollen_ts_bacon_v6.csv'))

load(paste0(data.dir,"pol.cal.count.mnwi3.Rdata"))

new.pol1<-cbind(pol.cal.count[pol.cal.count$dataset.id==unique(bacon$id[1]),],bacon[bacon$id==unique(bacon$id[1]),])

for(i in 2:26){
	if(nrow(pol.cal.count[pol.cal.count$dataset.id==unique(bacon$id)[i],])==nrow(bacon[bacon$id==unique(bacon$id)[i],])){
		new.pol<-cbind(pol.cal.count[pol.cal.count$dataset.id==unique(bacon$id)[i],],
	    bacon[bacon$id==unique(bacon$id)[i],])
	    new.pol1<-rbind(new.pol1,new.pol)
	}else{
		my.pol <- pol.cal.count[pol.cal.count$dataset.id==unique(bacon$id)[i],]
		bacon.pol <- bacon[bacon$id==unique(bacon$id)[i],]
		new.pol2 <- which(my.pol$Age%in%bacon.pol$age_default)
		new.pol<-cbind(my.pol[new.pol2,],bacon.pol[which(bacon.pol$age_default%in%my.pol$Age),])
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
		new.pol2 <- which(my.pol$Age%in%bacon.pol$age_default)
		new.pol<-cbind(my.pol[new.pol2,],bacon.pol[which(bacon.pol$age_default%in%my.pol$Age),])
	    new.pol1<-rbind(new.pol1,new.pol)
	    #print(unique(bacon$id)[i])
	}
}

save(new.pol1,file=paste0(data.dir,"baby.bacon.Rdata"))

#These aren't in pol.cal.count yet
#3081 dataset id
#CLH4   Lily04 45.901 -92.272 wisconsin
#CLH5   Lone02 45.932 -92.2359 wisconsin

#new.pol1[which((new.pol1$Age - new.pol1$age_default)==0&(new.pol1$QUERCUS - new.pol1$OAK)!=0),c("QUERCUS","OAK","Age","age_default","dataset.id")]

#unique(new.pol1[which((new.pol1$Age - new.pol1$age_default)!=0),'dataset.id'])

new.pol.final <- new.pol1


calc.all <- matrix(0,length(unique(new.pol1$dataset.id)),26)
for(i in 1:length(unique(new.pol1 $dataset.id))){
	    ag.keep<- new.pol1[new.pol1 $dataset.id==unique(new.pol1 $dataset.id)[i],1:6]
	    calc.all[i,1:6] <- unlist(ag.keep[1,])
		calc.all[i,7:26] <- colSums(new.pol1[new.pol1$dataset.id==unique(new.pol1$dataset.id)[i],7:26])
}

colnames(calc.all) <- colnames(new.pol1[,1:26])

calc.all <-as.data.frame(calc.all)


biomass_dat_est <- read.csv(paste0(data.dir,"biomass_prediction_v0.9-7_bam.csv"))

lat.long.reg <- cbind(calc.all$LongitudeWest, calc.all$LatitudeNorth)
lat.long.reg.df = data.frame(lat.long.reg)
colnames(lat.long.reg.df) = c('x', 'y')

coordinates(lat.long.reg.df) <- ~ x + y
proj4string(lat.long.reg.df) <- CRS('+proj=longlat +ellps=WGS84')

albers <- spTransform(lat.long.reg.df, CRS('+init=epsg:3175'))
albers <- as.matrix(data.frame(albers))

centers_biomass = cbind(biomass_dat_est$x,biomass_dat_est$y)
idx_cores = vector(length=nrow(calc.all))

for(i in 1:nrow(calc.all)){   
  core_site = albers[i,]
  d = rdist(matrix(core_site, ncol=2), as.matrix(centers_biomass))
  idx_cores[i] = which.min(d) 
}

if(DRAW == TRUE) pdf(paste0(dump.dir,"check_points.pdf")) else (quartz())
plot(albers[,1], albers[,2])
points(centers_biomass[idx_cores,1], centers_biomass[idx_cores,2], col='blue', pch=8)
plot(usShp, add=T, lwd=2)
#sites in michigan
if(DRAW == TRUE) dev.off()

calc.all <- cbind(calc.all,numeric(nrow(calc.all)))

for(i in 1:nrow(calc.all)){ 
  calc.all[i,27] = biomass_dat_est[idx_cores[i],"Total"]  
}

breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)-1))
data_binned <-  cut(calc.all[,27], c(breaks), include.lowest = FALSE, labels = FALSE)

pdf(paste0(fig.dir,paste0("biomass.pts.settlement.all.cores",Sys.Date(),".pdf")))
map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(calc.all$LongitudeWest, calc.all$LatitudeNorth, pch=21,
		cex=1.1, bg=colors[data_binned],lwd=.2)
title("Biomass Point Estiamtes at Settlement")
plotInset(-90,47,-82.5,50,
          expr={
          	hist(calc.all[,27],col=colors,xaxt="n",xlab=NA,
          	ylab=NA,main=NA,cex.lab=.5,cex.axis=.5,breaks=breaks)
          	axis(side=1,breaks,at = breaks,cex.axis = .5,las=2,line=0)
          	mtext(side = 1, "Biomass (Mg/ha)", line = 1.5,cex=.5)
         mtext(side = 2, "Frequency", line = 1.7,cex=.5)
          	})
dev.off()


#### making calibration dataset


settle <- new.pol.final[new.pol.final$age_bacon<=150&new.pol.final$age_bacon>=50,] # think harder about this part # problem because you don't plot every core's biomass at time of settlement



ag.settle <- matrix(0,length(unique(settle$dataset.id)),26)
for(i in 1:length(unique(settle$dataset.id))){
	    ag.keep<- settle[settle$dataset.id==unique(settle$dataset.id)[i],1:6]
	    ag.settle[i,1:6] <- unlist(ag.keep[1,])
		ag.settle[i,7:26] <- colSums(settle[settle$dataset.id==unique(settle$dataset.id)[i],7:26])
}

colnames(ag.settle) <- colnames(settle[,1:26])

ag.settle<-as.data.frame(ag.settle)

##### Changing pollen coordinates so that we can find the right rows when we find the biomass for each pond.

#### Biomass for calibration sites
biomass_dat_est <- read.csv(paste0(data.dir,"biomass_prediction_v0.9-7_bam.csv"))

lat.long.reg <- cbind(ag.settle$LongitudeWest,ag.settle$LatitudeNorth)
lat.long.reg.df = data.frame(lat.long.reg)
colnames(lat.long.reg.df) = c('x', 'y')

coordinates(lat.long.reg.df) <- ~ x + y
proj4string(lat.long.reg.df) <- CRS('+proj=longlat +ellps=WGS84')

albers <- spTransform(lat.long.reg.df, CRS('+init=epsg:3175'))
albers <- as.matrix(data.frame(albers))

centers_biomass = cbind(biomass_dat_est$x,biomass_dat_est$y)
idx_cores = vector(length=nrow(ag.settle))
idx_cores_all = vector(length=length(unique(new.pol.final$dataset.id)))

for(i in 1:nrow(ag.settle)){   
  core_site = albers[i,]
  d = rdist(matrix(core_site, ncol=2), as.matrix(centers_biomass))
  idx_cores[i] = which.min(d) 
}

if(DRAW == TRUE) pdf(paste0(dump.dir,"check_points.pdf")) else (quartz())
plot(albers[,1], albers[,2])
points(centers_biomass[idx_cores,1], centers_biomass[idx_cores,2], col='blue', pch=8)
plot(usShp, add=T, lwd=2)
#sites in michigan
if(DRAW == TRUE) dev.off()

ag.settle <- cbind(ag.settle,numeric(nrow(ag.settle)))

for(i in 1:nrow(ag.settle)){ 
  ag.settle[i,27] = biomass_dat_est[idx_cores[i],"Total"]  
}

breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)-1))
data_binned <-  cut(ag.settle[,27], c(breaks), include.lowest = FALSE, labels = FALSE)

pdf(paste0(fig.dir,paste0("biomass.pts.settlement",Sys.Date(),".pdf")))
map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(ag.settle$LongitudeWest,ag.settle$LatitudeNorth, pch=21,
		cex=1.1, bg=colors[data_binned],lwd=.2)
title("Biomass Point Estiamtes at Settlement")
plotInset(-90,47,-82.5,50,
          expr={
          	hist(ag.settle[,27],col=colors,xaxt="n",xlab=NA,
          	ylab=NA,main=NA,cex.lab=.5,cex.axis=.5,breaks=breaks)
          	axis(side=1,breaks,at = breaks,cex.axis = .5,las=2,line=0)
          	mtext(side = 1, "Biomass (Mg/ha)", line = 1.5,cex=.5)
         mtext(side = 2, "Frequency", line = 1.7,cex=.5)
          	})
dev.off()

#####
##### Remove sites and take subset of points
#####
#plot_biomass_pollen=plot_biomass_pollen[,-which(colSums(plot_biomass_pollen)==0)]

plot_biomass_pollen = cbind(ag.settle, albers)

if(DRAW == TRUE) pdf(paste0(fig.dir,paste0("all_sites",Sys.Date(),".pdf")))
map('state', ylim=range(plot_biomass_pollen[,2])+c(-2, 2), xlim=range(plot_biomass_pollen[,3])+c(-1, 1),main=NA)
points(plot_biomass_pollen[,3], plot_biomass_pollen[,2], pch=19, cex=1,col="gray")
title(main="all sites")

set.seed(4)
sites_rm = sample(1:nrow(plot_biomass_pollen),round(nrow(plot_biomass_pollen)/3))

points(plot_biomass_pollen[-sites_rm,3], plot_biomass_pollen[-sites_rm,2], pch=19, cex=1,col="blue")
if(DRAW == TRUE) dev.off()

library(mgcv)

biomass = plot_biomass_pollen[,27]
total_counts = round(rowSums(plot_biomass_pollen[,7:26]))
counts = round(plot_biomass_pollen[,7:26])
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

Z.knots.check = matrix(0,nrow=u[length(u)],ncol=(length(u)+2));

for(i in 1:u[length(u)]){
	u_given <-i
	Z.knots.check[i,] = bs_nimble(u_given, u=u, N0 = rep(0, (length(u)-1)),N1 = rep(0, (length(u))), N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))
}

for(i in 1:length(biomass)){
    u_given <- biomass[i]
	Z.knots[i,] = bs_nimble(u_given, u=u, N0 = rep(0, (length(u)-1)),N1 = rep(0, (length(u))), N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))
}

#### Plot basis functions ####
plot(Z.knots.check[,1],xlim=c(0,u[length(u)]),pch=19,ylim=c(0,1),xlab="Biomass")
for(i in 2:ncol(Z.knots.check)){
	points(Z.knots.check[,i],col=i,pch=19)
}
abline(v=u,lwd=2)
title("Basis Functions")


#if(DRAW == TRUE) dev.off()

#####
##### Create final datasets #####
#####

Y = counts[-sites_rm,] #remove sites "sites_rm" defined above
biomass = biomass[-sites_rm]
counts = counts[-sites_rm,]
total_counts = rowSums(counts)

final_coors = plot_biomass_pollen[-sites_rm,c(2,3,(ncol(plot_biomass_pollen)-1),ncol(plot_biomass_pollen))]

save.image(file="add.bacon.Rdata")





