
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
}

all.preds1 <- cbind(rep(lat,100),rep(lon,100),sort(rep(YBP,77)),c(agwb_time_slice_mean))
head(all.preds1)
colnames(all.preds1) <- c('lat','lon','age','agwb')

#### 
#### Fit GAM
#### 

b <- gam(log(agwb) ~ te(lon, lat, age, d = c(2,1),
                        bs = c("tp","cr"), k=50),
         data = as.data.frame(all.preds1))

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
coors_dat <- as.matrix(data.frame(lat.lon))

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

pdf('gam_maps_firstpass.pdf',compress = T)
par(mfrow=c(2,2))
for(age_slice in rev(seq(1000,10000,1000))){

  pred_data = cbind(coors_dat,rep(age_slice,nrow(coors_dat)))
  colnames(pred_data)<- c("lon","lat","age")
  pred_biomass_gam = exp(predict(b,newdata = as.data.frame(pred_data)))
  
  full.mat <- cbind(coors_dat,as.vector(pred_biomass_gam))
  colnames(full.mat) <- c("x","y","pred_biomass")

  breaks <-  c(seq(0,200,25), seq(300,600,100))
  colors <- rev(terrain.colors(length(breaks)-1))
  legendName <- c("Biomass (Mg/ha)")#paste0("Biomass at Age = ",age_slice, " BP")

  data_binned <-  cut(pred_biomass_gam, breaks, include.lowest = TRUE, labels = FALSE)
  
  breaklabels <- apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1,  function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })

  plot(full.mat[,1],full.mat[,2],col=colors[data_binned],pch=19)
  maps::map('state',add=T)
  title(age_slice)
  
  points_get <- x.meta[x.meta$age_bacon<(age_slice+100)&x.meta$age_bacon>(age_slice-100)&x.meta$site.name%in%dataID_use$name,c('lat','long')]
  points(points_get[,2],points_get[,1],pch=19,col='black',cex=.5)
  
  plot_count <- plot_count + 1
  if(plot_count == 3){
    plot.new()
    legend('center',breaklabels,pch=19,col=colors)
    plot_count <- 0
  }
}
dev.off()
### how much to cut off?
