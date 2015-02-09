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
require(maptools)
require(ggplot2)

usShp <- readShapeLines(file.path(data.dir, '/us_alb.shp'), proj4string=CRS('+init=epsg:3175'))
usShp@data$id <- rownames(usShp@data)
usFortified <- fortify(usShp, region='id')


#####
###### Load Biomass data #####
#####
biomass_dat <- read.csv(paste0(data.dir,"plss_biomass_v0.9-1.csv"))
biomass_dat[is.na(biomass_dat)] <-0 #take out NAs
biomass_dat_est <- read.csv(paste0(data.dir,"biomass_estimate_v1.9.csv"))

##### Get rid of weird groups in biomass data
biomass_dat<-biomass_dat[,colnames(biomass_dat)!="No.tree"]
biomass_dat<-biomass_dat[,colnames(biomass_dat)!="NonTree"]
biomass_dat<-biomass_dat[,colnames(biomass_dat)!="Unknown"]
biomass_dat<-biomass_dat[,colnames(biomass_dat)!="Unknown.tree"]


##### organize the dataframe
biomass_dat4<-rowSums(biomass_dat[,5:ncol(biomass_dat)])
centers_biomass = cbind(biomass_dat$x,biomass_dat$y)
biomass_dat5 <- data.frame(x=centers_biomass[,1],y=centers_biomass[,2],
                           as.numeric(rowSums(biomass_dat[,5:ncol(biomass_dat)])))               

colnames(biomass_dat5)<-c("x","y","biomass")
load(file=paste0(data.dir,"pol.cal.count.mnwi.csv")) #all MN and WI follows from MNWI_dat.R

hk_counts = read.csv(paste0(data.dir,"hotchkiss_lynch_calcote_counts_v0.csv")) #have to add someday
hk_meta = read.csv(paste0(data.dir,"hotchkiss_lynch_calcote_meta_v0.csv"))


x = pol.cal.count[pol.cal.count$Age>=100,]
x = x[x$Age<=200,]

##### Making a pollen proportion data frame that includes biomass for the pollen source grid cell.
melt.x <- melt(x,id.vars=c("SiteID","LatitudeNorth","LongitudeWest","dataset.id","ContactName"))
cast.x <- cast(melt.x,SiteID + LatitudeNorth + LongitudeWest + dataset.id +ContactName ~ variable,sum)
cast.x=as.data.frame(cast.x)
row_keep = rep(0,nrow(cast.x))
plot_biomass_pollen = matrix(0,nrow(cast.x),78) #75 for just bigwoods #78 for mnwi


##### Changing pollen coordinates so that we can find the right rows when we find the biomass for each pond.
new.site.locs <- cbind(cast.x$LongitudeWest,cast.x$LatitudeNorth)
centers_pol = data.frame(new.site.locs)
colnames(centers_pol) = c('x', 'y')

coordinates(centers_pol) <- ~ x + y
proj4string(centers_pol) <- CRS('+proj=longlat +ellps=WGS84')

centers_polA <- spTransform(centers_pol, CRS('+init=epsg:3175'))
centers_polA <- as.matrix(data.frame(centers_polA))

cast.x$LongitudeWest <- centers_polA[,1]
cast.x$LatitudeNorth <- centers_polA[,2]

centers_biomass = cbind(biomass_dat5$x,biomass_dat5$y)
idx_cores = vector(length=nrow(cast.x))

for(i in 1:nrow(cast.x)){   
  core_site = centers_polA[i,]
  d = rdist(matrix(core_site, ncol=2), as.matrix(centers_biomass))
  idx_cores[i] = which.min(d) 
}

if(DRAW == TRUE) pdf(paste0(dump.dir,"check_points.pdf"))
plot(centers_polA[,1], centers_polA[,2])
points(centers_biomass[idx_cores,1], centers_biomass[idx_cores,2], col='blue', pch=8)
plot(usShp, add=T, lwd=2) 
if(DRAW == TRUE) dev.off()

for(i in 1:nrow(cast.x)){ 
  plot_biomass_pollen[i,1] = sum(biomass_dat_est[idx_cores[i],])
  plot_biomass_pollen[i,2:78] = as.numeric(cast.x[i,c(2,3,7:ncol(cast.x))])
  #plot_biomass_pollen[i,4:78] = plot_biomass_pollen[i,4:78]/sum(plot_biomass_pollen[i,4:78])
}  

##### Fixing up the data frame
#hist(plot_biomass_pollen[,1],breaks=20)
#plot_biomass_pollen <- plot_biomass_pollen[-c(25,23),]#-c(25,23) #too much poaceae for bigwoods data set
colnames(plot_biomass_pollen)<-c("Biomass","LatNorth","LongWest",colnames(cast.x[7:ncol(cast.x)])) ####MUST RUN

plot_biomass_pollen=plot_biomass_pollen[,-which(colSums(plot_biomass_pollen)==0)]
if(SAVE == TRUE){
  save(plot_biomass_pollen,file=paste0(dump.dir,"plot_biomass_pollen.Rdata"))
}
#load("plot_biomass_pollen.Rdata")

colnames(plot_biomass_pollen)


##### Changing pollen coordinates back so we can plot the pie plots #### when you rerun from the beginning make sure this is nessecary...
new.site.locs <- cbind(plot_biomass_pollen[,3],plot_biomass_pollen[,2])
centers_pol = data.frame(new.site.locs)
colnames(centers_pol) = c('x', 'y')
coordinates(centers_pol) <- ~ x + y
proj4string(centers_pol) <- CRS('+init=epsg:3175')
centers_polA <- spTransform(centers_pol, CRS('+proj=longlat +ellps=WGS84'))
centers_polA <- as.matrix(data.frame(centers_polA))
plot_biomass_pollen[,2] <- centers_polA[,1]
plot_biomass_pollen[,3] <- centers_polA[,2]

if(DRAW == TRUE) pdf(paste0(dump.dir,"all_sites.pdf"))
map('state', xlim=range(plot_biomass_pollen[,2])+c(-2, 2), ylim=range(plot_biomass_pollen[,3])+c(-1, 1))
points(plot_biomass_pollen[,2], plot_biomass_pollen[,3], pch=19, cex=1)
points(cast.x[,3],cast.x[,2],col="white")
if(DRAW == TRUE) dev.off()

#head(plot_biomass_pollen)
library(mgcv)

biomass = plot_biomass_pollen[,1]
total_counts = round(rowSums(plot_biomass_pollen[,4:ncol(plot_biomass_pollen)]))
counts = round(plot_biomass_pollen[,4:ncol(plot_biomass_pollen)])
colnames(counts) <- colnames(plot_biomass_pollen[,4:ncol(plot_biomass_pollen)])

#going down to ncol(counts) spp for first attempt at model. also truncating biomass to 400.
trees <- c("POACEAE","PINUSX","CYPERACE","LARIXPSEU","TSUGAX","QUERCUS","TILIA","BETULA","PICEAX","OSTRYCAR","ULMUS","ABIES","POPULUS")
ten.count = matrix(0,142,length(trees)+1)
prairie <- c("AMBROSIA","ARTEMISIA","ASTERX","CHENOAMX","FABACEAE")
ten.count[,1] <- rowSums(counts[,prairie])
ten.count[,2:(length(trees)+1)] <- counts[,trees]
colnames(ten.count)<-c("PRAIRIE",trees)
ten.count = ten.count[-61,] #getting rid of grass pond
biomass = biomass[-61]
for(i in 1:length(biomass)){
  if(biomass[i]>400) biomass[i] = 400
}

library(splines)
Z = bs(biomass,intercept=TRUE) #add knots here

#plot(biomass,Z[,1])
#points(biomass,Z[,2],col="blue")
#points(biomass,Z[,3],col="red")
#points(biomass,Z[,4],col="green")
#points(biomass)

rownames(Z)<-NULL

library(boot)
library("mgcv")

counts = ten.count
total_counts = rowSums(counts)

Z = bs(biomass,intercept=TRUE,df=4)
betas = matrix(0,ncol(Z),ncol(counts))

if(DRAW == TRUE) pdf(paste0(dump.dir,"splines1.pdf"))
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
  
}
if(DRAW == TRUE) dev.off()

#Z.new%*%betas
library(gtools)
rownames(Z)<-NULL
delta = 50#bigger
phi.b = matrix(0,nrow(counts),ncol(counts)); p = phi.b
Y = matrix(0,nrow(counts),ncol(counts))
phi.b = exp(Z%*%betas)/rowSums(exp(Z%*%betas))

for(j in 1:nrow(counts)){
  p[j,] = rdirichlet(1,phi.b[j,]*delta)
  Y[j,] = rmultinom(1,prob = p[j,], size = rowSums(counts)[j])
}

colnames(Y)<-colnames(counts)
size = rowSums(Y)

save.image(paste0(dump.dir,"data_formatted.Rdata"))

print("Finished formatting data. Saved all data to data_formatted.Rdata")

set.seed(3)
sites_rm = sample(1:141,50)
Y = Y[-sites_rm,]
biomass = biomass[-sites_rm]
counts = counts[-sites_rm,]
total_counts = rowSums(counts)


if (nrow(Y) < 141) print("removed 50 sites for data validation")


