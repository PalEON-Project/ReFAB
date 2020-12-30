
#####
##### GAM Smoothing #####
#####

library(ncdf4)
library(lattice)
library(RColorBrewer)
library(mgcv)
library(ggplot2)
library(maptools)
library(reshape2)
library(maps)

#### 
#### Load reconstruction
#### 

nc <- nc_open('ReFAB_site_reconstruction_v1.0.nc')

agwb <- ncvar_get(nc,"AGWB")

nc$dim$lon$vals -> lon
nc$dim$lat$vals -> lat

cutpts <- seq(0,250,length.out = 20)
YBP <- ncvar_get(nc,"time")

#### 
#### Create mean data frame
#### 

agwb_time_slice_mean <- matrix(NA,77,100)

for(tt in 1:100){
  agwb_time_slice_mean[,tt] <- diag(apply(agwb[,,,tt],2,rowMeans))
  print(tt)
}

all.preds1 <- cbind(rep(lat,100),rep(lon,100),sort(rep(YBP,77)),c(agwb_time_slice_mean))
head(all.preds1)
colnames(all.preds1) <- c('lat','lon','age','agwb')

par(mfrow=c(1,2))
plot(all.preds1[,2],all.preds1[,1])
map('state',add=T)

cut_preds <- all.preds1[-which(all.preds1[,1]< 42.6 | all.preds1[,2]> -84.67906),]
cut_preds1 <- cut_preds[-which(cut_preds[,1]< 46 & cut_preds[,2]> -86),]
plot(cut_preds1[,2],cut_preds1[,1])
map('state',add=T)

#### 
#### Fit GAM
#### 

b <- gam(log(agwb) ~ te(lon, lat, age, d = c(2,1),
                        bs = c("tp","cr"), k = 50),
         data = as.data.frame(cut_preds1))

gam.check(b)

summary(b)
vis.gam(b)  

#### 
##### Get prediction coordiantes from PLS map
#### 

nc_pls <- nc_open(file.path('~','Downloads','PLS_biomass_agb_western_point_v1.0rc1.nc'))

x <- nc_pls$dim$x$vals
y <- nc_pls$dim$y$vals
data <- ncvar_get(nc_pls,varid = c('Total'))

rownames(data)  <- x
colnames(data)  <- y

r1 <- raster(list(x=x,y=y,z=data))
plot(r1)

foo <- rasterToPoints(r1)

x <- foo[,1]
y <- foo[,2]

#### 
#### Change from Albers to lat / lon
#### 

library(raster)
albers.df = as.data.frame(cbind(x,y))
colnames(albers.df) = c('x', 'y')

coordinates(albers.df) <- ~ x + y
proj4string(albers.df) <- CRS('+init=epsg:3175')

lat.lon <- spTransform(albers.df, CRS('+proj=longlat +ellps=WGS84'))

r2 <- raster(lat.lon)
e2 <- extent(-86.86787, -82.49863, 41.71928, 45.78518)
r2c <- crop(lat.lon, e2)

coors_dat <- as.matrix(data.frame(lat.lon))
coors_remove <- as.matrix(data.frame(r2c))

row_keep <- list()
for(i in 1:nrow(coors_remove)){
  row_keep[i] <- which(coors_dat[,1]==coors_remove[i,1])
  if(length(row_keep[i])>1){
    row_keep[i] <- which(coors_dat[,1]==coors_remove[i,1]&coors_dat[,2]==coors_remove[i,2])
  }
}

row_remove <- unlist(row_keep)

coors_dat_minus_smi <- coors_dat[-row_remove,]
which(coors_dat_minus_smi[,1]>-83)
coors_dat_minus_smi_1 <- coors_dat_minus_smi[-which(coors_dat_minus_smi[,1]>-83),]


latlong2state <- function(pointsDF) {
  # Prepare SpatialPolygons object with one SpatialPolygon
  # per state (plus DC, minus HI & AK)
  states <- map('state', fill=TRUE, col="transparent", plot=FALSE)
  IDs <- sapply(strsplit(states$names, ":"), function(x) x[1])
  states_sp <- map2SpatialPolygons(states, IDs=IDs,
                                   proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # Convert pointsDF to a SpatialPoints object 
  pointsSP <- SpatialPoints(pointsDF, 
                            proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # Use 'over' to get _indices_ of the Polygons object containing each point 
  indices <- over(pointsSP, states_sp)
  
  # Return the state names of the Polygons object containing each point
  stateNames <- sapply(states_sp@polygons, function(x) x@ID)
  stateNames[indices]
}

states_yay <- latlong2state(data.frame(x = coors_dat_minus_smi_1[,1], y = coors_dat_minus_smi_1[,2]))

coors_dat_keep <- coors_dat_minus_smi_1[states_yay%in%c('minnesota','wisconsin','michigan'),]
plot(coors_dat)
points(coors_dat_keep,col='red')

usShp <- readShapeLines(file.path('Data', '/us_alb.shp'), proj4string=CRS('+init=epsg:3175'))
usShp@data$id <- rownames(usShp@data)
usFortified <- fortify(usShp, region='id')

load('prediction.data_v6.Rdata')
dataID <- read.csv('dataID_v5.csv')

dataID_use <- dataID[-which(dataID$name%in%c('Lily Lake','Mud Lake','Chatsworth Bog')),]

#### 
#### Plotting function for differnt time slices
#### 
plot_count <- 0

pdf('gam_maps_secondpass.pdf',height=12,width = 12,compress = T)
par(mfrow=c(2,2))
for(age_slice in rev(seq(1000,10000,1000))){

  pred_data = cbind(coors_dat_keep,rep(age_slice,nrow(coors_dat_keep)))
  colnames(pred_data)<- c("lon","lat","age")
  pred_biomass_gam = exp(predict(b,newdata = as.data.frame(pred_data)))
  
  full.mat <- cbind(coors_dat_keep,as.vector(pred_biomass_gam))
  colnames(full.mat) <- c("x","y","pred_biomass")

  breaks <-  c(seq(0,200,25), seq(300,600,100))
  colors <- rev(terrain.colors(length(breaks)-1))
  legendName <- c("Biomass (Mg/ha)")#paste0("Biomass at Age = ",age_slice, " BP")

  data_binned <-  cut(pred_biomass_gam, breaks, include.lowest = TRUE, labels = FALSE)
  
  breaklabels <- apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1,  function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })

  plot(full.mat[,1],full.mat[,2],col=colors[data_binned],pch=19)
  maps::map('state',add=T)
  title(age_slice)
  
  points_get <- x.meta[x.meta$age_bacon<(age_slice+100)&x.meta$age_bacon>(age_slice-100)&x.meta$site.name%in%dataID_use$name&x.meta$lat%in%cut_preds1[,1],c('lat','long')]
  points(points_get[,2],points_get[,1],pch=8,col='blue',cex=1.5)
  
  points_cols <- cut_preds1[which(cut_preds1[,'age']==age_slice),c('lon','lat','agwb')]
  
  data_binned_pts <-  cut(points_cols[,'agwb'], breaks, include.lowest = TRUE, labels = FALSE)
  points(points_cols[,1],points_cols[,2],pch=21,bg=colors[data_binned_pts],col='black',cex=1)
  
  legend('topright',c('avg bacon age +- 100 years',
                      'ReFAB avg pt est'),
         pch = c(8,21),
         col = c('blue','black'),
         pt.bg = c(NA,colors[5]))
  
  plot_count <- plot_count + 1
  if(plot_count == 3){
    plot.new()
    legend('center',breaklabels,pch=19,col=colors)
    plot_count <- 0
  }
}
dev.off()
### how much to cut off?
