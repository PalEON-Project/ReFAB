
library(geoR)
library(raster)
library(maptools)
library(ncdf4)
library(maps)

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

#### 
#### Crop Illinois and Indiana
#### 
foo <- foo[which(foo[,2]>561553.9),]
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
#### Load ReFAB Output
#### 

refab <- read.csv('median_biomass_plus_meta.csv')

library(spatstat)

ch <- convexhull.xy(x= refab$lon,y=refab$lat)

#### 
#### Variogram
####  

pdf('refab_krig_exponential.pdf',height=15,width = 15)
par(mfrow=c(3,3),mar=c(0,0,4,0))
for(tt in rev(seq(10,90,10))){

time_do <- colnames(refab)[2:101][tt]
time_slice <- data.frame(refab[, c('lon', 'lat', time_do)],X=rep(1,77))
ut.gd <- as.geodata(time_slice, coords.col =
                      1:2, data.col = 3,covar.col = 4)

ut.v = variog(ut.gd,
              trend =  ~ refab[, c('lon')] + refab[, c('lat')],
              max.dist = max(dist(refab[, c('lon', 'lat')]))/3)

ut.wls = variofit(
  ut.v,
  ini = c(400, 6),
  cov.model = "exponential",
  fix.nug = FALSE,
  nugget = 1500,
  wei = "cressie"
)


if(TRUE){
  #par(mfrow=c(1,1))
  plot(ut.v,pch=19)
  lines(ut.wls) 
  title('Semi-Variogram')
}

#### 
#### Krige
#### 

ut.krig = krige.conv(
  ut.gd,
  locations = coors_dat,
  borders=NULL,
  krige = krige.control(
    type.krige = "OK",
    obj.m = ut.wls
  )
)

grid <- data.frame(coors_dat)
grid$z <- ut.krig$predict

breaks <- breaks <-  c(seq(0,50,10),seq(75,250,25),435)
data_binned <- cut(grid$z,breaks=breaks,labels=F,include.lowest=F)
colors <- rev(terrain.colors(length(breaks)-1))

data_binned_pts <- cut(time_slice[,3],breaks=breaks,labels=F)

plot(grid$x,grid$y,col=adjustcolor(colors[data_binned],alpha.f = .5),pch=19)
plot(ch,add=T)
points(time_slice$lon,time_slice$lat,pch=21,bg = colors[data_binned_pts],cex=2)
map('state',add=T)
title(paste(tt*100,'Year BP Prediction'))
legend('topright',legend=breaks,col=colors,pch=19)

var_breaks <- seq(0,3000,500)
var_data_binned <- cut(ut.krig$krige.var,breaks=var_breaks,labels=F)
white_red <- colorRampPalette(c('white','bisque','red','darkorchid4'))
var_colors <- (white_red(length(var_breaks)))


plot(grid$x,grid$y,col=var_colors[var_data_binned],pch=19)
plot(ch,add=T)
points(time_slice$lon,time_slice$lat,pch=21,bg = 'black')
map('state',add=T)
legend('topright',legend=var_breaks,col=var_colors,pch=19)
title('Variance')


}
dev.off()
