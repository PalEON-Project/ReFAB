library(raster)
library(ncdf4)

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


#### 
#### Change from Albers to lat / lon
#### 

proj4string(r1) <- CRS('+init=epsg:3175')
lat.lon <- projectRaster(r1, crs = ('+proj=longlat +ellps=WGS84'))

coors_dat <- rasterToPoints(lat.lon)

coors_dat <- cbind(coordinates(lat.lon), v=values(lat.lon))#as.matrix(data.frame(lat.lon))

coors_dat[is.na(coors_dat)] <- 0

breaks <- seq(0,1000,100)
data_binned <- as.numeric(cut(coors_dat[,3],breaks))
plot(coors_dat[,1],coors_dat[,2],col=data_binned)

colnames(coors_dat) <- c('lon','lat','TotalBiomass_Mgha')

write.csv(coors_dat,file='paleon_total_biomass_settlement.csv')

#loop over all vars
all_types <- list()
for(i in 1:22){
  q <- ncvar_get(nc_pls,names(nc_pls$var)[i])
  rownames(q)  <- x
  colnames(q)  <- y
  r2 <- raster(list(x=x,y=y,z=q))
  proj4string(r2) <- CRS('+init=epsg:3175')
  r3 <- projectRaster(r2, crs = ('+proj=longlat +ellps=WGS84'))
  all_types[[i]] <- p <- rasterToPoints(r3)
  
  breaks <- seq(0,1000,50)
  data_binned <- as.numeric(cut(p[,3],breaks))
  plot(p[,1],p[,2],col=data_binned)
}

all_mat <- cbind(do.call(cbind,lapply(all_types,FUN=function(x)x[,3])),all_types[[1]][,1:2])
colnames(all_mat) <- c(names(nc_pls$var),'lon','lat')

write.csv(all_mat,file='paleon_byspp_biomass_settlement.csv')



