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


#####
###### Load biomass and pollen data #####
#####

#### Biomass
biomass_dat_est <- read.csv(paste0(data.dir,"biomass_prediction_v0.9-7_bam.csv"))

#### Pollen
load(file="~/stepps-data/data/pol.cal.count.mnwi3.Rdata") #all MN and WI and MI part that we want follows from MNWI_dat.R

#x[x$SiteID==29,] #there are two siteID 29s

#####
###### Format settlement pollen data #####
#####
x = pol.cal.count[pol.cal.count$Age>=100,]
x = x[x$Age<=200,]

melt.x <- melt(x,id.vars=c("SiteID","LatitudeNorth","LongitudeWest","dataset.id","ContactName"))
cast.x <- cast(melt.x,LatitudeNorth + LongitudeWest + SiteID +ContactName ~ variable,sum)
cast.x=as.data.frame(cast.x)
cast.x = as.data.frame(cast.x[-c(1:3),])
row_keep = rep(0,nrow(cast.x))
plot_biomass_pollen = matrix(0,nrow(cast.x),(ncol(x)-3))


##### Changing pollen coordinates so that we can find the right rows when we find the biomass for each pond.
lat.long.reg <- cbind(cast.x$LongitudeWest,cast.x$LatitudeNorth)
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
  idx_cores[i] = which.min(d) 
}

if(DRAW == TRUE) pdf(paste0(dump.dir,"check_points.pdf")) else (quartz())
plot(albers[,1], albers[,2])
points(centers_biomass[idx_cores,1], centers_biomass[idx_cores,2], col='blue', pch=8)
plot(usShp, add=T, lwd=2)
if(DRAW == TRUE) dev.off()

for(i in 1:nrow(cast.x)){ 
  plot_biomass_pollen[i,1] = sum(biomass_dat_est[idx_cores[i],5:25])
  plot_biomass_pollen[i,2:(ncol(x)-3)] = as.numeric(cast.x[i,c(1,2,6:ncol(cast.x))])
  #plot_biomass_pollen[i,4:78] = plot_biomass_pollen[i,4:78]/sum(plot_biomass_pollen[i,4:78])
}

breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)-1))
data_binned <-  cut(plot_biomass_pollen[,1], c(breaks), include.lowest = FALSE, labels = FALSE)

	 
pdf(paste0(fig.dir,paste0("biomass.pts.settlement",Sys.Date(),".pdf")))
map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(plot_biomass_pollen[,3],plot_biomass_pollen[,2], pch=21,
		cex=1.1, bg=colors[data_binned],lwd=.2)
title("Biomass Point Estiamtes at Settlement")
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
##### Remove sites and take subset of points
#####
colnames(plot_biomass_pollen)<-c("Biomass","LatNorth","LongWest",colnames(cast.x[6:ncol(cast.x)])) ####MUST RUN

#plot_biomass_pollen=plot_biomass_pollen[,-which(colSums(plot_biomass_pollen)==0)]
plot_biomass_pollen = cbind(plot_biomass_pollen, albers)

if(DRAW == TRUE) pdf(paste0(dump.dir,paste0("all_sites",Sys.Date(),".pdf")))
quartz()
par(mfrow=c(1,2))
map('state', ylim=range(plot_biomass_pollen[,2])+c(-2, 2), xlim=range(plot_biomass_pollen[,3])+c(-1, 1))
points(plot_biomass_pollen[,3], plot_biomass_pollen[,2], pch=19, cex=1)
title(main="all sites")

set.seed(4)
sites_rm = sample(1:nrow(plot_biomass_pollen),round(nrow(plot_biomass_pollen)/3))

map('state', xlim=range(plot_biomass_pollen[,3])+c(-2, 2), ylim=range(plot_biomass_pollen[,2])+c(-1, 1))
points(plot_biomass_pollen[-sites_rm,3], plot_biomass_pollen[-sites_rm,2], pch=19, cex=1)
title(main="remaining sites")
if(DRAW == TRUE) dev.off()

library(mgcv)

biomass = plot_biomass_pollen[,1]
total_counts = round(rowSums(plot_biomass_pollen[,4:(ncol(plot_biomass_pollen)-2)]))
counts = round(plot_biomass_pollen[,4:(ncol(plot_biomass_pollen)-2)])
colnames(counts) <- colnames(plot_biomass_pollen[,4:(ncol(plot_biomass_pollen)-2)])

props = counts/rowSums(counts)

total_counts_spp = colSums(counts)

props = props[,order(total_counts_spp,decreasing=TRUE)]

pdf("scatter.newdata.pdf")
#quartz()
par(mfrow=c(4,4))
for(i in 1:ncol(props)){
  if(length(unique(props[-sites_rm,i]))>=9){
    plot(plot_biomass_pollen[-sites_rm,1],props[-sites_rm,i],main=colnames(props)[i],xlab="biomass",ylab="pollen prop",pch = 19, cex = .5)
  }

}  
dev.off()

props = as.data.frame(props)

#plot_biomass_pollen = plot_biomass_pollen[-c(which(props$Other>.5),which(props$POACEAE>.8)),]

#####
##### Creating a dataset with the species we want to use
#####

# counts = counts[,-which(colnames(counts)==c("PINUSX"))]
# trees <- c("ALNUSX","JUGLANSX","ACERX","CUPRESSA","FRAXINUX","FAGUS","CYPERACE","LARIXPSEU","TSUGAX","QUERCUS","TILIA","BETULA","PICEAX","OSTRYCAR","ULMUS","ABIES","POPULUS")
# other.trees <- c("TAXUS","NYSSA","CASTANEA","PLATANUS","SALIX","LIQUIDAM")
# ten.count = matrix(0,nrow(counts),length(trees)+3)
# prairie <- c("ARTEMISIA","ASTERX","POACEAE","AMBROSIA","CHENOAMX","CORYLUS")
# ten.count[,1] <- rowSums(counts[,prairie])
# ten.count[,2] <- rowSums(counts[,other.trees])
# ten.count[,3:(length(trees)+2)] <- counts[,trees]
# ten.count[,(length(trees)+3)] <- rowSums(counts) - rowSums(ten.count)
# colnames(ten.count)<-c("prairie","other trees",trees,"other herbs")

# props = as.data.frame(props)

ten.count = count
total_counts = rowSums(counts)

#####
##### Drawing Splines #####
#####

Z = bs(biomass,intercept=TRUE,df=4)
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
  new.biomass = seq(1,400,1)
  Z.new = bs(new.biomass,intercept=TRUE,df = ncol(Z))
  lines(new.biomass, predict(glm_mod,newdata=list(Z=Z.new),type="response"),col="blue")  
  
  betas[,i] = glm_mod$coefficients
  betas.save[,i] = glm_mod$coefficients
  
}
#if(DRAW == TRUE) dev.off()

#####
##### Simulating Data #####
#####
# 
# rownames(Z)<-NULL
# delta = 50#bigger
# phi.b = matrix(0,nrow(counts),ncol(counts)); p = phi.b
# Y = matrix(0,nrow(counts),ncol(counts))
# phi.b = exp(Z%*%betas)/rowSums(exp(Z%*%betas))
# 
# for(j in 1:nrow(counts)){
#   p[j,] = rdirichlet(1,phi.b[j,]*delta)
#   Y[j,] = rmultinom(1,prob = p[j,], size = rowSums(counts)[j])
# }
# 
# colnames(Y)<-colnames(counts)
# size = rowSums(Y)
# 
# print("Finished formatting data. Saved all data to data_formatted.Rdata")

#####
##### Create final datasets #####
#####

Y = counts[-sites_rm,] #remove sites "sites_rm" defined above
biomass = biomass[-sites_rm]
counts = counts[-sites_rm,]
total_counts = rowSums(counts)

final_coors = plot_biomass_pollen[-sites_rm,c(2,3,(ncol(plot_biomass_pollen)-1),ncol(plot_biomass_pollen))]
