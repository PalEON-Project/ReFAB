
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
library(raster)

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
#### Load Smaller ReFAB Output
#### 

refab <- read.csv('median_biomass_plus_meta.csv')

library(reshape2)

refab_melt <- melt(refab,id.vars = c('lat','lon','name','site_index','X','cluster'))

head(refab_melt)
levels(refab_melt$variable) <- 1:100
refab_melt$variable<-as.numeric(refab_melt$variable)

#### 
#### Fit GAM - with same data just different input files
#### 

#long data format
if(FALSE)
b <- gam(log(agwb) ~ te(lon, lat, age, d = c(2,1),
                        bs = c("tp","cr"), k=50),
         data = as.data.frame(all.preds1))

#short data format
tictoc::tic()
b <- gam(value ~ te(lon, lat, variable, d = c(2,1),
                        bs = c("tp","cr"), k=70), #can choose k up to 77
         data = as.data.frame(refab_melt),
         control = gam.control(nthreads = 4))
tictoc::toc()

summary(b)
vis.gam(b)  

#### 
##### Get prediction coordiantes from PLS map
#### 

nc_pls <- nc_open(file.path('~','Downloads','PLS_biomass_agb_western_point_v1.0rc2.nc'))

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
#### Crop Illinois and Indiana
#### 
foo <- foo[which(foo[,2]>463528),]
#foo <- foo[-which(foo[,1]<664541.9&foo[,2]<690790.5),]

x <- foo[,1]
y <- foo[,2]

plot(x,y)

#### 
#### Change from Albers to lat / lon
#### 

albers.df = as.data.frame(cbind(x,y))
colnames(albers.df) = c('x', 'y')

coordinates(albers.df) <- ~ x + y
proj4string(albers.df) <- CRS('+init=epsg:3175')

lat.lon <- spTransform(albers.df, CRS('+proj=longlat +ellps=WGS84'))
coors_dat <- as.matrix(data.frame(lat.lon))

#### 
#### Plotting function for differnt time slices
#### 
plot_count <- 0

library(spatstat)

ch <- convexhull.xy(x= refab$lon,y=refab$lat)


pdf('gam_maps_thirdpass.pdf',compress = T,height=15,width = 20)

#layout(matrix(c(1,2,9,10,3,4,11,12,5,6,13,14,7,8,15,16,17,18,19,20),4,5))

layout(matrix(c(1,2,3,10,4,5,6,10,7,8,9,11),3,4,byrow=T))

par(oma=c(2,2,0,0),mar=c(0,0,4,0))

for(age_slice in rev(seq(10,80,10))){
print(age_slice)
  pred_data = cbind(coors_dat,rep(age_slice,nrow(coors_dat)))
  colnames(pred_data)<- c("lon","lat","variable")#c("lon","lat","age")
  pred_biomass_gam = (predict(b,newdata = as.data.frame(pred_data),type='response',se.fit=TRUE))
  
  full.mat <- cbind(coors_dat,as.vector(pred_biomass_gam$fit))
  colnames(full.mat) <- c("x","y","pred_biomass")

  breaks <-  c(seq(0,50,10),seq(75,250,25),435)
  colors <- rev(terrain.colors(length(breaks)-1))
  data_binned <-  cut(pred_biomass_gam$fit, breaks, include.lowest = TRUE, labels = FALSE)
  
  legendName <- c("Biomass (Mg/ha)")#paste0("Biomass at Age = ",age_slice, " BP"
  breaklabels <- apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1,  function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })
  
  
  plot(full.mat[,1],full.mat[,2],
       col=adjustcolor(colors[data_binned],
                       alpha.f = .25),
       pch=19,
       xaxt='n',
       yaxt='n')
  maps::map('state',add=T)
  plot(ch,add=T)
  title(paste('Biomass',age_slice*100,'YBP'),cex=3)
  
  #legend('topright',breaklabels,pch=19,col=colors)
 
#### REFAB POINTS 
  #points_get <- x.meta[x.meta$age_bacon<(age_slice+100)&x.meta$age_bacon>(age_slice-100)&x.meta$site.name%in%dataID_use$name,c('lat','long')]
  pt_data <- refab_melt[refab_melt$variable==age_slice,]
  data_binned_pts <- cut(pt_data[,'value'],breaks=breaks,labels=F)
  points(pt_data[,2],pt_data[,1],pch=21,col='black',bg=colors[data_binned_pts],cex=3)
  
  

  
}


#### VARIANCE PLOTS  
breaks_se <-  c(seq(0,10,1))
white_red <- colorRampPalette(c('white','bisque','red','darkorchid4'))
colors_se <- (white_red(length(breaks_se)))
data_binned_se <-  cut(pred_biomass_gam$se, breaks_se, include.lowest = TRUE, labels = FALSE)
breaklabels_se <- apply(cbind(breaks_se[1:(length(breaks_se)-1)], breaks_se[2:length(breaks_se)]), 1,  function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })

plot(full.mat[, 1], full.mat[, 2],
     col = colors_se[data_binned_se],
     pch = 19,
     xaxt='n',
     yaxt='n')
maps::map('state',add=T)
plot(ch,add=T)
title(paste('GAM Standard Error'))

points(pt_data[,2],pt_data[,1],pch=21,bg='black',cex=1)


plot.new()
legend('center',breaklabels,pch=19,col=colors,cex=2,title = 'Biomass (Mg/ha)')

plot.new()
legend('center',legend=breaklabels_se,col=colors_se,pch=19,cex=2,title = 'SE')

dev.off()
### how much to cut off?
