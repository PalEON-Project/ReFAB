library(neotoma)
library(dplyr)
library(maptools)#
library(ggplot2)
library(reshape)#
library(fields)
library(oce)#
library(splines)
library(mgcv)
library(ncdf4)#
library(raster)
library(plyr)

data.dir = c("Data/")
fig.dir = c('Figures/')
DRAW = TRUE

### Shape files for map plots
usShp <- readShapeLines(file.path(data.dir, '/us_alb.shp'), proj4string=CRS('+init=epsg:3175'))
usShp@data$id <- rownames(usShp@data)
usFortified <- fortify(usShp, region='id')

#####
##### Download Data for MN and WI and MI #####
#####


meta <- get_dataset(datasettype='pollen', gpid=c("Wisconsin", 
                                                 "Michigan",
                                                 "Minnesota",
                                                 "Illinois",
                                                 "Indiana"), ageyoung=0)

meta_dl <- get_download(meta)

comp.tax <- compile_taxa(meta_dl, 'WhitmoreSmall')
pol_cal_count <- compile_downloads(comp.tax)
pol_cal_count <- pol_cal_count %>% mutate(.id = as.numeric(.id))
pol_cal_count_save <- pol_cal_count

save(pol_cal_count,file=paste0('nimble_pull',Sys.Date(),'.Rdata'))
load("~/ReFAB/nimble_pull2018-02-06.Rdata")
load("~/ReFAB/nimble_pull2018-07-13.Rdata")
load('~/ReFAB/nimble_pull2018-10-31.Rdata')

if(DRAW==TRUE) pdf(file.path(fig.dir,paste0('all_pollen_cores_map',Sys.Date(),'.pdf')))
map('state', ylim=range(pol_cal_count$lat)+c(4, 2), xlim=range(pol_cal_count$long)+c(-1, 1),main=NA)
points(pol_cal_count$long, pol_cal_count$lat, pch=19, cex=.4,col="black")
#title(paste0('N = ',length(unique(pol_cal_count$dataset))))
if(DRAW== TRUE) dev.off()

#####
##### Add Bacon Pollen Dates #####
#####

cal <- read.csv(paste0(data.dir,'/cal_data_mid_depth_2015-06-10.csv')) #andria calibration
hotch <- read.csv('hk_counts3.csv') #adding hotchkiss sites separately because they aren't in neotoma
pol_cal_count <- pol_cal_count %>% mutate(.id = as.numeric(.id))
cal.new.pol1 <- merge(x = pol_cal_count, y = cal, 
                     by.x = c('dataset','depth'), 
                     by.y = c('id','depth'))
colnames(cal.new.pol1)[(9:10)] <- c('lat','long')
not_found <- cal[-which(cal$id%in%cal.new.pol1$dataset),]
cal.new.pol2 <- merge(x = hotch, y = not_found, 
                      by.x = c('lat','long','BETULA','PICEAX'), 
                      by.y = c('lat','long','BIRCH','SPRUCE'))
cal.new.pol <- rbind.fill(cal.new.pol1,cal.new.pol2)

####
#### Expert Elicitation
####

JWW <- read.csv('Data/elicitation_additional_data_UMW_JWW.csv')
STJ <- read.csv('Data/elicitation_additional_data_UMW_STJ.csv')

experts <- droplevels(cbind(JWW[,1:4],STJ[,4]))
expert_cal <- list()#matrix(NA,nrow(experts),ncol=ncol(pol_cal_count))
for(i in 1:nrow(experts)){
  pol_rows <- which(as.character(pol_cal_count$site.name) == droplevels(experts[i,1]))
  pol_site <- pol_cal_count[pol_rows,]
  
  depth_get <- as.numeric(experts[i,4])
  if(depth_get!=0) expert_cal[[i]] <- pol_site[depth_get,]
  
}

expert.pol <- do.call(rbind,expert_cal)

####
#### Map to look at difference of old versus new site coverage
####

pdf(file.path(fig.dir,'new.sites.pdf'))
map('state', ylim=range(cal.new.pol1$lat)+c(-2, 2),
    xlim=range(cal.new.pol1$long)+c(-2, 2),main=NA)
points(cal.new.pol$long, cal.new.pol$lat, pch=19, cex=1,col="red")
points(cal.new.pol1$long, cal.new.pol1$lat, pch=19, cex=1,col="black")
points(expert.pol$long, expert.pol$lat, pch=19, cex=1,col="blue")
legend('topright',c('Original, Neotoma','Hotchkiss','New from EE'),
       col=c('black','red','blue'),pch=19)
dev.off()

#####
###### Format settlement pollen data #####
#####

x <- plyr::rbind.fill(cal.new.pol,expert.pol)

length(unique(x$dataset))

all.pollen.taxa.names <- colnames(pol_cal_count)[11:length(colnames(pol_cal_count))]
melt.x <- melt(x, id.vars=c('.id','age','lat','long','site.name'),
               measure.vars = all.pollen.taxa.names)
cast.x <- cast(melt.x, lat + long + .id + site.name ~ variable, sum) #summing over pollen taxa for years around settlement
cast.x <- cast.x.full <- as.data.frame(cast.x)
save(cast.x,file='cast.x.Rdata')

#####
##### Remove sites and take subset of points
#####

set.seed(3)
sites_rm = sample(1:nrow(cast.x),round(nrow(cast.x)/3))
save(sites_rm, file = paste0(Sys.Date(),'sites_rm.Rdata'))

if(DRAW == TRUE) pdf(paste0(fig.dir,paste0("all_sites",Sys.Date(),".pdf")))
map('state', ylim=range(pol_cal_count$lat)+c(-2, 2), xlim=range(pol_cal_count$long)+c(-1, 1),main=NA)
points(cast.x$long, cast.x$lat, pch=19, cex=1,col="gray")
points(cast.x$long[-sites_rm], cast.x$lat[-sites_rm], pch=19, cex=1,col="blue")
title(main="all sites")
if(DRAW == TRUE) dev.off()

#plot_biomass_pollen = plot_biomass_pollen[-c(which(props$Other>.5),which(props$POACEAE>.8)),]

#####
##### Adding settlement biomass data #####
#####

###old ###biomass_dat_est <- read.csv('Data/biomass_prediction_v0.9-10_bam.csv')
###stem
#nc <- nc_open(file.path(data.dir,'PLS_biomass_western_point_v0.999.nc'))
#nc_draws <- nc_open(file.path(data.dir,'PLS_biomass_western_v0.999.nc'))

###AGB
nc <- nc_open(file.path('~','Downloads','PLS_biomass_agb_western_point_v1.0rc1.nc'))
nc_draws <- nc_open(file.path('~','Downloads','PLS_biomass_agb_western_v1.0rc1.nc'))


x <- nc$dim$x$vals
y <- nc$dim$y$vals
data <- ncvar_get(nc,varid = c('Total'))
data_draws <- ncvar_get(nc_draws,varid = c('Total'))


rownames(data) <- rownames(data_draws) <- x
colnames(data) <- colnames(data_draws) <- y

r1 <- raster(list(x=x,y=y,z=data))

#can do
#plot(r1)

biomass_dat_est <- as.data.frame(rasterToPoints(r1))
save(biomass_dat_est,file='biomass_dat_est_agb.Rdata')
colnames(biomass_dat_est) <- c('x','y','Total')

##### Changing pollen coordinates so that we can find the right rows when we find the biomass for each pond.
lat.long.reg <- cbind(as.numeric(as.character(cast.x$long)),as.numeric(as.character(cast.x$lat)))
lat.long.reg.df = data.frame(lat.long.reg)
colnames(lat.long.reg.df) = c('x', 'y')

coordinates(lat.long.reg.df) <- ~ x + y
proj4string(lat.long.reg.df) <- CRS('+proj=longlat +ellps=WGS84')

albers <- spTransform(lat.long.reg.df, CRS('+init=epsg:3175'))
albers <- as.matrix(data.frame(albers))

save(albers,file=paste0(Sys.Date(),'calibration.albers.Rdata'))

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

### Map checking for biomass grid cell that pollen core is contained within
if(DRAW == TRUE) pdf(paste0(fig.dir,"check_points.pdf")) else (quartz())
plot(albers[,1], albers[,2])
points(centers_biomass[idx_cores,1], centers_biomass[idx_cores,2],
       col='blue', pch=8)
plot(usShp, add=T, lwd=2)
if(DRAW == TRUE) dev.off()

biomass <- list()
biomass_draws <- list()
for(i in 1:length(idx_cores)){ 
  biomass[[i]] = biomass_dat_est[idx_cores[i],'Total']
  biomass_draws[[i]] <- data_draws[as.character(biomass_dat_est[idx_cores[i],'x']),as.character(biomass_dat_est[idx_cores[i],'y']),]
}

#only need to save for full calibration so no need to take to bottom of script
save(biomass_draws,file='biomass_draws_v3.0.Rdata')

breaks <-  c(seq(0,50,10),seq(75,250,25),435)
colors <- (tim.colors(length(breaks)-1))
data_binned <-  cut(biomass_dat_est$Total, c(breaks), include.lowest = FALSE, labels = FALSE)

data_binned_biomass <- cut(unlist(biomass), c(breaks), include.lowest = FALSE, labels = FALSE)

### checking to make sure the biomass we are getting is from the closest grid cell
if(DRAW == TRUE) pdf('look_biomass.pdf')
plot(centers_biomass,col=colors[data_binned],cex=.4,pch=19)
legend('topright',legend=breaks,col=colors,pch=19)
points(albers,cex=.3,pch=23,bg=colors[data_binned_biomass],col='black',lwd=.2)
if(DRAW == TRUE) dev.off()

if(DRAW == TRUE) pdf('settlement.rainbow.pdf')
plot(r1,col=tim.colors(32))
points(albers[,1],albers[,2],col='black')
if(DRAW == TRUE) dev.off()

cast.x <- cbind(cast.x,unlist(biomass))

#### Remove places where we do not have settlement biomass # should have all of it eventually
#### cast.x <- cast.x[-which(is.na(cast.x[,ncol(cast.x)])),]
#### biomass <- biomass[-which(is.na(biomass))]

biomass <- unlist(biomass)

breaks <-  c(seq(0,40,10),seq(50,300,50))
colors <- rev(terrain.colors(length(breaks)-1))
data_binned <-  cut(cast.x[,ncol(cast.x)], c(breaks), include.lowest = FALSE, labels = FALSE)

if(DRAW==TRUE) pdf(paste0(fig.dir,paste0("biomass.pts.settlement",Sys.Date(),".pdf")))

pdf('[3]only_calibration_map_agb.pdf')
par(mfrow=c(1,1))
map('state', xlim=c(-100,-80), ylim=c(39.5,49.5))
points(cast.x$long,cast.x$lat, pch=21,
       cex=1.1, bg=colors[data_binned],lwd=.2)
#text(cast.x$long,cast.x$lat,labels=1:nrow(cast.x),cex=.3)
#bimodal_sites <- c(13,31,35,36,15,34,61,11,14,18,21,33)
#points(cast.x$long[-sites_rm][bimodal_sites],cast.x$lat[-sites_rm][bimodal_sites],col='red',lwd=2)
#title("Biomass Point Estimates at Settlement")
plotInset(-90,47,-82.5,50,
          expr={
            hist(data_binned,col=colors,xaxt="n",xlab=NA,
                 ylab=NA,main=NA,cex.lab=.5,cex.axis=.5)
            axis(side=1,breaks,at = seq(1,length(breaks),1),cex.axis = .5,las=2,line=0)
            mtext(side = 1, "Biomass (Mg/ha)", line = 1.5,cex=.5)
            mtext(side = 2, "Frequency", line = 1.7,cex=.5)
          })
dev.off()
if(DRAW==TRUE) dev.off()

#####
##### Creating a calibration dataset with the species we want to use
#####

trees <- c("JUGLANSX","FRAXINUX","OSTRYCAR","ULMUS","TILIA","CARYA",
           "FAGUS","TSUGAX","QUERCUS","BETULA",
           'PINUSX',"ACERX","ALNUSX",
           "CYPERACE","PICEAX","ABIES","POPULUS",
           "LARIXPSEU","CUPRESSA") #
other.trees <- c("CASTANEA","PLATANUS","SALIX","LIQUIDAM","TAXUS","NYSSA")#NULL#c()
drop.taxa <- NA#c('other_herbs')

source(file.path('Workflow_Code','utils','taxa_selection.R'))
Y <- taxa_selection(trees = trees, other.trees = other.trees,
                    cast.x = cast.x, sites_rm = 0,
                    all.pollen.taxa.names = all.pollen.taxa.names,
                    prairie.include = T,bigwoods.include=F, other.herbs.include = T,
                    other.trees.include = T, drop.taxa = drop.taxa,
                    PFT.do = F)

Y.all <- Y
total_counts = round(rowSums(Y,na.rm = TRUE))
props = Y/rowSums(Y,na.rm = TRUE)

total_counts_spp = colSums(Y)

props = props[,order(total_counts_spp,decreasing=TRUE)]

if(DRAW==TRUE) pdf(paste0(fig.dir,"scatter.newdata",Sys.Date(),".pdf"))
#quartz()
par(mfrow=c(2,2))
for(i in 1:ncol(props)){
  if(length(unique(props[,i]))>=9){
    plot(biomass,props[,i],main=colnames(props)[i],xlab="biomass",ylab="pollen prop",pch = 19, cex = .5)
    calibrate::textxy(biomass,props[,i],1:length(biomass),cex=.2)
       }
}  
if(DRAW==TRUE) dev.off()

props = as.data.frame(props)

#####
##### Drawing Splines (with all data)
#####

Z = Z.knots = bs(unlist(biomass), intercept=TRUE, knots = median(biomass), Boundary.knots=c(0,max(biomass)))
u <- c(0,median(biomass),max(biomass)) #c(rep(attr(Z.knots,"Boundary.knots")[1],1),attr(Z.knots,"knots"),rep(attr(Z.knots,"Boundary.knots")[2],1))
betas = matrix(0,ncol(Z),ncol(Y)); betas.save = betas

counts <- Y
counts[is.na(counts)] <- 0 ### important. Could put above in merge or ten.count creation.

if(DRAW == TRUE) pdf(paste0(fig.dir,"splines.new.pdf"))
#quartz()
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

if(DRAW == TRUE) dev.off()

library(nimble)
source('Workflow_Code/utils/bs_nimble.R')

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
if(DRAW==TRUE) pdf(file.path(fig.dir,'basis_function_check.pdf'))
par(mfrow=c(1,1))
plot(Z.knots.check[,1],xlim=c(0,u[length(u)]),pch=19,ylim=c(0,1),xlab="Biomass")
for(i in 2:ncol(Z.knots.check)){
  points(Z.knots.check[,i],col=i,pch=19)
}
abline(v=u,lwd=2)
title("Basis Functions")
if(DRAW==TRUE) dev.off()

#####
##### Create final calibration datasets #####
#####

save(Y,biomass,Z,file=paste0(Sys.Date(),'all.calibration.data.Rdata'))


#####
##### Getting the original golden 1/3 back per Chris's demands
##### eventhough the counts are different now because we are using settlement horizon
##### when we originally developed the model we were using age-depth model estimated age at time of settelement
#####
Y.save <- Y
biomass.save <- biomass
albers.save <- albers
sites_rm_save <- sites_rm
load("~/ReFAB/2018-03-01calibration.albers.Rdata")
albers <- albers[-c(1,2,6),] #removing southern MI sites
TF <- list()
for(i in 1:nrow(albers)){
  TF[[i]] <-  which(rdist(x1=matrix(albers[i,],ncol=2),x2=as.matrix(albers.save))==min(rdist(x1=matrix(albers[i,],ncol=2),x2=as.matrix(albers.save)),na.rm = T))
}

TF <- unique(unlist(TF))

plot(albers.save[,1],albers.save[,2]) #new dataset
points(albers[,1],albers[,2],pch=8,col='blue') #old dataset
points(albers.save[TF,1],albers.save[TF,2],pch=1,col='green') #matching new dataset to old dataset
load('sites_rm.Rdata')
TF_keep <- cbind(albers,TF)[-sites_rm,]
points(TF_keep[,1],TF_keep[,2],pch=1,col='red') #old 2/3s datasset
save(TF,file='TF.Rdata')
old_full <- Y.save[TF,]
old_full_biomass <- biomass.save[TF]

old_two_thirds <- old_full[-sites_rm,]
old_two_thirds_biomass <- old_full_biomass[-sites_rm]

new_pool <- Y.save[-TF,]
new_pool_biomass <- biomass.save[-TF]

set.seed(4)
sites_rm_new <- sample(1:nrow(new_pool),round(nrow(new_pool)/3)) #just removing from new_pool that's why it's small
save(sites_rm_new,file='new_sites_rm.Rdata')

new_two_thirds <- new_pool[-sites_rm_new,]
new_two_thirds_biomass <- new_pool_biomass[-sites_rm_new]

ag_two_thirds <- rbind(old_two_thirds,new_two_thirds)
ag_two_thirds_biomass <- c(old_two_thirds_biomass,new_two_thirds_biomass)


old_albers <- albers.save[TF,]
old_two_thirds_albers <- old_albers[-sites_rm,]
new_albers <- albers.save[-TF,]
new_two_thirds_albers <- new_albers[-sites_rm_new,]
ag_two_thirds_albers <- rbind(old_two_thirds_albers,new_two_thirds_albers)

cast.x.save <- cbind(cast.x,1:232)
old.cast.x <- cast.x[TF,]
old.two.thirds.cast.x <- old.cast.x[-sites_rm,]
old.one.third.cast.x <- old.cast.x[sites_rm,]
new.cast.x <- cast.x[-TF,]
new.two.thirds.cast.x <- new.cast.x[-sites_rm_new,]
new.one.third.cast.x <- new.cast.x[sites_rm_new,]
ag.two.thirds.cast.x <- rbind(old.two.thirds.cast.x,new.two.thirds.cast.x)
ag.one.thirds.cast.x <- rbind(old.one.third.cast.x,new.one.third.cast.x)

one.third.idx <- ag.one.thirds.cast.x[,ncol(ag.one.thirds.cast.x)]

save(one.third.idx,file=paste0(Sys.Date(),'one.third.idx.Rdata'))

save(ag.two.thirds.cast.x,file=paste0(Sys.Date(),'two.thirds.cast.x.Rdata'))

#Y = Y[-sites_rm,] #remove sites "sites_rm" defined above
#biomass = biomass[-sites_rm]

Y <- ag_two_thirds
biomass <- ag_two_thirds_biomass
Z = Z.knots = bs(unlist(biomass), intercept=TRUE, knots = median(biomass), Boundary.knots=c(0,max(biomass)))
dim(Y)[1]-length(biomass) # should be zero

save(Y,biomass,file=paste0(Sys.Date(),'twothirds.calibration.data.Rdata'))

### moved prediction dataset creation to [2]create_prediction_dataset.R
