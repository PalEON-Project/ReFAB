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

diff(order((pol_cal_count %>% filter(.id == 1966))$depth))

map('state', ylim=range(pol_cal_count$lat)+c(-2, 2), xlim=range(pol_cal_count$long)+c(-1, 1),main=NA)
points(pol_cal_count$long, pol_cal_count$lat, pch=19, cex=1,col="gray")

save(pol_cal_count,file=paste0('nimble_pull',Sys.Date(),'.Rdata'))
load("~/ReFAB/nimble_pull2018-02-06.Rdata")
load("~/ReFAB/nimble_pull2018-07-13.Rdata")
#####
##### Add Bacon Pollen Dates #####
#####

bacon <- read.csv(paste0(data.dir,'/sediment_ages_v1.0_varves.csv')) #andria bacon
cal <- read.csv(paste0(data.dir,'/cal_data_mid_depth_2015-06-10.csv')) #andria calibration

pol_cal_count <- pol_cal_count %>% mutate(.id = as.numeric(.id))

#current calibration dataset -> subset of andria's dataset
cal.new.pol <- merge(x = pol_cal_count, y = cal, 
                     by.x = c('dataset','depth'), 
                     by.y = c('id','depth'))

pol_all_settle <- pol_cal_count[which(pol_cal_count$age > 50 & pol_cal_count$age < 250),]
#possible new calibration sites
pol_all_settle_nocal <-pol_all_settle[-which(pol_all_settle$dataset%in%cal$id),]

### Plot sites that are currently missing from calibration data set
pdf('calibration.sites.pdf')
map('state', ylim=range(cal$lat)+c(-2, 2),
    xlim=range(cal$long)+c(-1, 1),main=NA)

points(cal$long, cal$lat, pch=19, col="blue",cex=.5)
points(cal.new.pol$long.x, cal.new.pol$lat.x, cex=.5,col="green")
points(pol_all_settle_nocal$long,
       pol_all_settle_nocal$lat,cex=.5,pch=19,col='red')

legend('topright',c('Andria cal','Ann current cal','Neotoma Age b/t 50-250'),
       pch=c(19,1,19), col= c('blue','green','red'))
dev.off()

new.pollen <- merge(x = pol_cal_count, y = bacon, 
                    by.x = c('.id','age','depth','BETULA'), 
                    by.y = c('id','age_default','depth','BIRCH'))

aa <- bacon %>% 
  dplyr::select(id, age_default, depth, starts_with('bacon_draw')) %>% 
  left_join(pol_cal_count,
            by = c('id' = '.id', 'age_default' = 'age', 'depth'))


pol_list <- list()
pol_ids <- unique(pol_cal_count$dataset)
for(i in 1:length(pol_ids)){
  pol_list[[i]] <- pol_cal_count[pol_cal_count$dataset==pol_ids[i],]
}
names(pol_list) <- pol_ids

### Test
sum(aa$BIRCH - aa$BETULA) #should equal zero
new.pollen[which(new.pollen$PINE != new.pollen$PINUSX),c('depth','PINE',"PINUSX")]
#####
###### Format settlement pollen data #####
#####

#original
#x <- new.pollen[new.pollen$age_bacon >= 100, ]
#x <- x[x$age_bacon<=250,]

#just settlement horizon
#colnames(cal.new.pol)[9] <- c('lat')
#colnames(cal.new.pol)[10] <- c('long')
#x <- cal.new.pol

#settlement horizon + previous 100 years
pol_settle <- list()
for(i in 1:length(pol_list)){
  #cal_list <- merge(x = pol_list[[i]], y = cal, 
  #                 by.x = c('dataset','depth'), 
  #                by.y = c('id','depth'))
  cal_stop <- cal[which(cal$id==names(pol_list)[i]),]
  cal_list <- pol_list[[i]][which(pol_list[[i]]$depth==cal_stop$depth),]
  # this is not using bacon ages I think... Ask Andria if cal has bacon ages. 
  #We need the bacon ages to aggregate 100 years back from settlement horizon.
  i_cal <- which(pol_list[[i]]$age >= cal_list$age & pol_list[[i]]$age < (cal_list$age + 100))
  pol_settle[[i]] <- pol_list[[i]][i_cal,]
  
}
## might not be getting all of the sites
x <- pol_settle_mat <- do.call(rbind,pol_settle)
length(unique(x$dataset))

all.pollen.taxa.names <- colnames(pol_cal_count)[11:length(colnames(pol_cal_count))]
melt.x <- melt(x, id.vars=c('.id','age','lat','long','site.name'),
               measure.vars = all.pollen.taxa.names)
cast.x <- cast(melt.x, lat + long + .id + site.name ~ variable, sum) #summing over pollen taxa for years around settlement
cast.x <- cast.x.full <- as.data.frame(cast.x)


#####
##### Remove sites and take subset of points
#####

set.seed(4)
sites_rm = sample(1:nrow(cast.x),round(nrow(cast.x)/3))
#save(sites_rm, file = paste0(Sys.Date(),'sites_rm.Rdata'))

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

biomass_dat_est <- read.csv(paste0(data.dir,"biomass_prediction_v0.9-10_bam.csv"))
nc <- nc_open(file.path(data.dir,'PLS_biomass_western_point_v0.99.nc'))

x <- nc$dim$x$vals
y <- nc$dim$y$vals
data <- ncvar_get(nc,varid = c('total'))

rownames(data) <- x
colnames(data) <- y

r1 <- raster(list(x=x,y=y,z=data))
plot(r1)

biomass_dat_est <- as.data.frame(rasterToPoints(r1))
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
points(centers_biomass[idx_cores,1], centers_biomass[idx_cores,2], col='blue', pch=8)
plot(usShp, add=T, lwd=2)
if(DRAW == TRUE) dev.off()

biomass <- list()
for(i in 1:length(idx_cores)){ 
  biomass[[i]] = biomass_dat_est[idx_cores[i],'Total']
}

cast.x <- cbind(cast.x,unlist(biomass))

#### Remove places where we do not have settlement biomass # should have all of it eventually
#### cast.x <- cast.x[-which(is.na(cast.x[,ncol(cast.x)])),]
#### biomass <- biomass[-which(is.na(biomass))]

biomass <- unlist(biomass)

breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)-1))
data_binned <-  cut(cast.x[,ncol(cast.x)], c(breaks), include.lowest = FALSE, labels = FALSE)

pdf(paste0(fig.dir,paste0("biomass.pts.settlement",Sys.Date(),".pdf")))
map('state', xlim=c(-98,-81), ylim=c(41.5,49))
points(cast.x$long,cast.x$lat, pch=21,
       cex=1.1, bg=colors[data_binned],lwd=.2)
#text(cast.x$long[-sites_rm],cast.x$lat[-sites_rm],labels=1:100,cex=.3)
#bimodal_sites <- c(3,98,62,100,73,92,63,55, 76, 4, 87, 79, 51)
#points(cast.x$long[-sites_rm][bimodal_sites],cast.x$lat[-sites_rm][bimodal_sites],col='red',lwd=2)
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

### I do this after indexing because some of the places where we have settlement pollen we don't have settlement biomass.

trees <- c("FAGUS","TSUGAX","QUERCUS","BETULA",
           'PINUSX',"JUGLANSX","ACERX","FRAXINUX",
           "OSTRYCAR","ULMUS","TILIA","ALNUSX",
           "CYPERACE","PICEAX","ABIES","POPULUS",
           "CARYA","LARIXPSEU","CUPRESSA") #
other.trees <- c("CASTANEA","PLATANUS","SALIX","LIQUIDAM","TAXUS","NYSSA")#NULL#c()
drop.taxa <- NA#c('other_herbs')

source(file.path('Workflow_Code','utils','taxa_selection.R'))
Y <- taxa_selection(trees = trees, other.trees = other.trees,
                    cast.x = cast.x, sites_rm = 0,
                    all.pollen.taxa.names = all.pollen.taxa.names,
                    prairie.include = T, other.herbs.include = T,
                    other.trees.include = T, drop.taxa = drop.taxa,
                    PFT.do = F)

set.seed(5)
sites_pred <- sample(1:nrow(Y),round(nrow(Y)/3))

Y.all <- Y
#Y <- Y[-sites_pred,]


total_counts = round(rowSums(Y,na.rm = TRUE))
props = Y/rowSums(Y,na.rm = TRUE)

total_counts_spp = colSums(Y)

props = props[,order(total_counts_spp,decreasing=TRUE)]

if(DRAW==TRUE) pdf(paste0(fig.dir,"scatter.newdata",Sys.Date(),".pdf"))
#quartz()
par(mfrow=c(4,4))
for(i in 1:ncol(props)){
  if(length(unique(props[-sites_rm,i]))>=9){
    plot(biomass[-sites_rm],props[-sites_rm,i],main=colnames(props)[i],xlab="biomass",ylab="pollen prop",pch = 19, cex = .5)
  }
}  
if(DRAW==TRUE) dev.off()

props = as.data.frame(props)

#####
##### Drawing Splines (with all data) ##### Not super important just for checking
#####

Z = Z.knots = bs(unlist(biomass), intercept=TRUE, knots = 30, Boundary.knots=c(0,232))
u <- c(0,30,232) #c(rep(attr(Z.knots,"Boundary.knots")[1],1),attr(Z.knots,"knots"),rep(attr(Z.knots,"Boundary.knots")[2],1))
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

Y = Y[-sites_rm,] #remove sites "sites_rm" defined above
biomass = biomass[-sites_rm]
Z = Z.knots = bs(unlist(biomass), intercept=TRUE, knots = 30, Boundary.knots=c(0,232))
dim(Y)[1]-length(biomass) # should be zero

save(Y,biomass,Z,file=paste0(Sys.Date(),'twothirds.calibration.data.Rdata'))

#####
##### Create final prediction datasets (BACON) #####
#####

x = new.pollen[new.pollen$age_bacon>=200,]
x = x[x$age_bacon<=10000,]

x.meta = x[,c('id','lat',"long","dataset","site.name","age_bacon")]
x.bacon <- x[,grep(pattern = 'bacon',colnames(x))]
colnames(x.meta)[1] <- c('site.id')

pred.x <- x[,which(colnames(x)%in%colnames(cast.x))]

ten.count <- taxa_selection(trees = trees, other.trees = other.trees,
                    cast.x = pred.x, sites_rm = 0,
                    all.pollen.taxa.names = all.pollen.taxa.names,
                    prairie.include = T, other.herbs.include = T,
                    other.trees.include = T, drop.taxa = drop.taxa,
                    PFT.do = F)
load("~/ReFAB/twothirds_v1.0.Rdata")
ten.count <- ten.count[,colnames(Y)]

#save(x.meta, ten.count, x.bacon , Z, u, file = 'prediction.data_v1.Rdata')


#### Writing data frame to feed to job array
IDs <- unique(x.meta$site.id)
name.list <- list()
site.data <- as.data.frame(matrix(NA,length(IDs),5))

source(file.path('Workflow_Code','utils','test_site.R'))
for(i in 1:length(IDs)){ 
  which_rows <- which(x.meta$site.id == IDs[i])
  x.meta.use <- x.meta[which_rows,]
  test_site(x.meta.use)
  
  ten_count_use = ten.count[which_rows, ]
  
  site.data[i,] <- data.frame(max.age = max(x.meta[which_rows,'age_bacon']),
                              min.age = min(x.meta[which_rows,'age_bacon']),
                              lat = x.meta[which_rows,'lat.x'][1],
                              long = x.meta[which_rows,'long.x'][1],
                              n.samps = length(which_rows))
  name.list[[i]] <-  x.meta[which_rows,'site.name'][1]
  
}

site.data <- cbind(site.data,unlist(name.list),unlist(IDs))
colnames(site.data) <- c('max.age','min.age','lat','long','n.samps','site.name','site.id')

site_keep <- site.data[which(site.data$max.age>8000&site.data$min.age<2000&site.data$n.samps>10),]

map('state', ylim=range(pol_cal_count$lat)+c(-2, 2), xlim=range(pol_cal_count$long)+c(-1, 1),main=NA)
#points(pol_cal_count$long, pol_cal_count$lat, pch=19, cex=1,col="gray")
#points(site.data$long, site.data$lat, pch=19, cex=1,col="blue")
points(site_keep$long,site_keep$lat,col='red',lwd=1,pch=19)

site_keep_bacon <- site_keep

if(which(table(site_keep$site.name)>1)) {
  print('Site Doubled Up.')
  print(table(site_keep$site.name)[which(table(site_keep$site.name)>1)])
  }

n.sites <- nrow(site_keep)
n.betas <- 20

dataID <- data.frame(name = sort(rep(site_keep$site.name,n.betas)), ID = 1:n.sites,
                     sigma = rep(0.12,n.sites*n.betas), beta = rep(1:n.betas,n.sites))

write.csv(dataID, file='dataID_bacon_v3.csv')

#####
##### Create final prediction datasets (NON - BACON) #####
#####

nonbacon <- pol_cal_count[-which(pol_cal_count$.id%in%unique(x$.id)),]

x <- nonbacon

x.meta = x[,c('.id','lat',"long","dataset","site.name","age")]
colnames(x.meta)[1] <- c('site.id')

pred.x <- x[,which(colnames(x)%in%colnames(cast.x))]

ten.count <- taxa_selection(trees = trees, other.trees = other.trees,
                            cast.x = pred.x, sites_rm = 0,
                            all.pollen.taxa.names = all.pollen.taxa.names,
                            prairie.include = T, other.herbs.include = T,
                            other.trees.include = T, drop.taxa = drop.taxa,
                            PFT.do = F)
ten.count <- ten.count[,colnames(Y)]

save(x.meta, ten.count, Z, u, file = 'prediction.data_v3_nonbacon.Rdata')

#### Writing data frame to feed to job array
IDs <- unique(x.meta$site.id)
name.list <- list()
site.data <- as.data.frame(matrix(NA,length(IDs),5))

source(file.path('Workflow_Code','utils','test_site.R'))
for(i in 1:length(IDs)){ 
  which_rows <- which(x.meta$site.id == IDs[i])
  x.meta.use <- x.meta[which_rows,]
  test_site(x.meta.use)
  
  ten_count_use = ten.count[which_rows, ]
  
  site.data[i,] <- data.frame(max.age = max(x.meta[which_rows,'age']),
                              min.age = min(x.meta[which_rows,'age']),
             lat = x.meta[which_rows,'lat'][1],
             long = x.meta[which_rows,'long'][1],
             n.samps = length(which_rows))
  name.list[[i]] <-  x.meta[which_rows,'site.name'][1]
  
}

site.data <- cbind(site.data,unlist(name.list),unlist(IDs))
colnames(site.data) <- c('max.age','min.age','lat','long','n.samps','site.name','site.id')

site_keep <- site.data[which(site.data$max.age>8000&site.data$min.age<2000&site.data$n.samps>10),]

pdf('new.prediction.sites.pdf')
map('state', ylim=range(site_keep$lat)+c(-2, 2), xlim=range(site_keep$long)+c(-3, 3),main=NA)
#points(pol_cal_count$long, pol_cal_count$lat, pch=19, cex=1,col="gray")
#points(site.data$long, site.data$lat, pch=19, cex=1,col="blue")
points(site_keep_bacon$long,site_keep_bacon$lat,col='blue',lwd=1,pch=19)
points(site_keep$long,site_keep$lat,col='red',lwd=1,pch=19)
legend('topright',c('Baconized','NOT Baconized'),pch=c(19,19),
       col=c('blue','red'))
dev.off()

if(which(table(site_keep$site.name)>1)) {
  print('Site Doubled Up.')
  print(table(site_keep$site.name)[which(table(site_keep$site.name)>1)])
}

n.sites <- nrow(site_keep)
n.betas <- 20

dataID <- data.frame(name = sort(rep(site_keep$site.name,n.betas)), ID = 1:n.sites,
                     sigma = rep(0.12,n.sites*n.betas), beta = rep(1:n.betas,n.sites))




#####
##### Create paleon mip datasets #####
#####

x = new.pollen[new.pollen$age_bacon>=200,]
x = x[x$age_bacon<=2000,]

x.meta = x[,c('.id','lat',"long","dataset","site.name","age_bacon")]
colnames(x.meta)[1] <- c('site.id')

ten.count = matrix(0,nrow(x),length(trees)+3)
ten.count[,1] <- unlist(rowSums(x[,prairie],na.rm = TRUE))
ten.count[,2] <- unlist(rowSums(x[,other.trees],na.rm = TRUE))
ten.count[,3:(length(trees)+2)] <- as.matrix(x[,trees])
ten.count[,(length(trees)+3)] <- rowSums(x[,all.pollen.taxa.names],na.rm = TRUE) - rowSums(ten.count,na.rm = TRUE)
colnames(ten.count)<-c("prairie","other_trees",trees,"other_herbs")

ten.count.save = ten.count
ten.count = round(ten.count.save)
ten.count <- ten.count[,colnames(counts)]


save(u,Z,x.meta,ten.count,file = 'paleon.data.Rdata')


#####
##### Plots #####
#####

if(DRAW==TRUE) pdf(paste0(fig.dir,"all.sites.neotoma",Sys.Date(),".pdf"))
map('state', xlim=range(x.meta$long)+c(-2, 2), ylim=range(x.meta$lat)+c(-1, 1))
points(x.meta$long, x.meta$lat, pch=19, cex=1,col="black")
title("All Pollen Sites")
if(DRAW==TRUE) dev.off()

if(DRAW==TRUE) pdf(paste0(fig.dir,"chrono.10k",Sys.Date(),".pdf"))
par(mfrow=c(1,1))
plot(x.meta$age_bacon,x.meta$lat,pch=19,cex=.5,xlim=c(0,10000),main="All Pollen Records",ylab="Latitude",xlab="Age BP")

plot.seq = seq(0,10000,500)

par(mfrow = c(2,2))
for(i in 2:length(plot.seq)){
  map('state', xlim=range(as.numeric(as.character(x.meta$long)))+c(-2, 2), ylim=range(as.numeric(as.character(x.meta$lat)))+c(-1, 1))
  points(x.meta[x.meta$age>plot.seq[i-1]&x.meta$age<plot.seq[i],]$long, 
         x.meta[x.meta$age>plot.seq[i-1]&x.meta$age<plot.seq[i],]$lat, 
         pch=19, cex=.5)
  title(c(plot.seq[i-1],"-",plot.seq[i]))
}
if(DRAW==TRUE) dev.off()

