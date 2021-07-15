
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
library(raster)
library(spatstat)
library(plotrix)

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
if(FALSE){
b <- gam(log(agwb) ~ te(lon, lat, age, d = c(2,1),
                        bs = c("tp","cr"), k = 50),
         data = as.data.frame(cut_preds1))

gam.check(b)

#short data format
b <- gam(log(value) ~ te(lon, lat, variable, d = c(2,1),
                        bs = c("tp","cr"), k=40),
         data = as.data.frame(refab_melt))

#short data format
tictoc::tic()
b <- gam(value ~ te(lon, lat, variable, d = c(2,1),
                        bs = c("tp","cr"), k=40), #can choose k up to 77
         data = as.data.frame(refab_melt))
tictoc::toc()

save(b,file='b.Rdata')
}

summary(b)
vis.gam(b)

#### 
##### Get prediction coordiantes from FIA map
#### 

nc_fia <- nc_open(file.path('~','Downloads','FIA_biomass_agb_point_v1.0rc2.nc'))

x_fia <- nc_fia$dim$x$vals
y_fia <- nc_fia$dim$y$vals
data_fia <- ncvar_get(nc_fia,varid = c('Total'))

rownames(data_fia)  <- x_fia
colnames(data_fia)  <- y_fia

r1_fia <- raster(list(x=x_fia,y=y_fia,z=data_fia))
plot(r1_fia)

foo_fia <- rasterToPoints(r1_fia)

x_fia <- foo_fia[,1]
y_fia <- foo_fia[,2]

albers.df_fia = data.frame(x=x_fia,y=y_fia)

coordinates(albers.df_fia) <- ~ x + y
proj4string(albers.df_fia) <- CRS('+init=epsg:3175')

lat.lon_fia <- spTransform(albers.df_fia, CRS('+proj=longlat +ellps=WGS84'))
coors_dat_fia <- as.matrix(data.frame(lat.lon_fia,fia_biomass=foo_fia[,3]))

###
### GET modern prediction from ESA 
###

esa <- read.csv('~/Downloads/ESAonPEONgrid.csv')

breaks <-  c(0,25, 50,75,100,125,150,175,200,225,250,300,800)#c(0, seq(25, 100, 25), seq(150, 300, 50), max(y$Tot_Biomass))
#breaks <-  c(0,50,100,150,200,250, max(y$Tot_Biomass))#c(0, seq(25, 100, 25), seq(150, 300, 50), max(y$Tot_Biomass))
colors <- rev(terrain.colors(length(breaks))[-4])

esa$ESA_AGWB[esa$PEON_AGWB>999] <- NA
cuts <- cut(esa$ESA_AGWB,breaks=breaks,include.lowest=T,labels=F)

plot(esa$X,esa$Y,pch=19,col=colors[cuts],main = 'ESA 2010 biomass estimate')

albers.df_esa = data.frame(x=esa$X,y=esa$Y)

coordinates(albers.df_esa) <- ~ x + y
proj4string(albers.df_esa) <- CRS('+init=epsg:3175')

lat.lon_esa <- spTransform(albers.df_esa, CRS('+proj=longlat +ellps=WGS84'))
coors_dat_esa <- as.matrix(data.frame(lat.lon_esa,biomass=esa$ESA_AGWB))

cuts <- cut(coors_dat_esa[,3],breaks=breaks,include.lowest=T,labels=F)

plot(coors_dat_esa[,1],coors_dat_esa[,2],pch=19,col=colors[cuts],main = 'ESA 2010 biomass estimate')


# esa$PEON_AGWB[esa$PEON_AGWB>999] <- NA
# cuts <- cut(esa$PEON_AGWB,breaks=breaks,include.lowest=T,labels=F)
# plot(esa$X,esa$Y,pch=19,col=colors[cuts])


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

foo_pls <- rasterToPoints(r1)

x <- foo_pls[,1]
y <- foo_pls[,2]

x_pls <- foo_pls[,1]
y_pls <- foo_pls[,2]

albers.df_pls = data.frame(x=x_pls,y=y_pls)

coordinates(albers.df_pls) <- ~ x + y
proj4string(albers.df_pls) <- CRS('+init=epsg:3175')

lat.lon_pls <- spTransform(albers.df_pls, CRS('+proj=longlat +ellps=WGS84'))
coors_dat_pls <- as.matrix(data.frame(lat.lon_pls,fia_biomass=foo_pls[,3]))

ove <- over(albers.df_esa,albers.df_pls)
ove_pls <- over(albers.df_pls,albers.df_esa)


diff_do <-
  coors_dat_esa[which(!is.na(ove)), 3] - coors_dat_pls[which(!is.na(ove_pls)), 3]
diff_breaks <- seq(-400, 400, 10)
diff_cuts <-
  cut(
    diff_do,
    breaks = diff_breaks,
    include.lowest = T,
    labels = F
  )
diff_colors <-
  colorRampPalette(c('purple', 'blue', 'white', 'red', 'maroon'))(length(diff_breaks) -
                                                                    1)
pdf('ESA_v_PLS.pdf',height = 10,width = 13)
layout(matrix(c(1, 1, 2, 2, 0,
                0, 3, 3, 4, 4), 2, 5, byrow = T))

cuts <-
  cut(
    coors_dat_esa[which(!is.na(ove)), 3],
    breaks = breaks,
    include.lowest = T,
    labels = F
  )
plot(
  coors_dat_esa[which(!is.na(ove)), 1],
  coors_dat_esa[which(!is.na(ove)), 2],
  pch = 15,
  col = colors[cuts],
  main = 'ESA 2010 biomass estimate',
  ylab = NA,
  xlab = NA,
  cex=.75
)
maps::map('state',add=T)

cuts <-
  cut(
    coors_dat_pls[which(!is.na(ove_pls)), 3],
    breaks = breaks,
    include.lowest = T,
    labels = F
  )
plot(
  coors_dat_pls[which(!is.na(ove_pls)), 1],
  coors_dat_pls[which(!is.na(ove_pls)), 2],
  pch = 15,
  col = colors[cuts],
  main = 'PLS biomass estimate',
  ylab = NA,
  xlab = NA,
  cex=.75
)
maps::map('state',add=T)
legend('topright',
       as.character(breaks[seq(1,length(breaks)-1,length.out = 6)]),
       col = colors[seq(1,length(breaks)-1,length.out = 6)],
       pch = 15)


plot(
  coors_dat_esa[which(!is.na(ove)), 1],
  coors_dat_esa[which(!is.na(ove)), 2],
  pch = 15,
  col = diff_colors[diff_cuts],
  main = 'ESA minus PLS',
  ylab = NA,
  xlab = NA,
  cex=.75
)
maps::map('state',add=T)
legend('topright',
       as.character(pretty(diff_breaks, n = 8)),
       col = diff_colors[c(which(diff_breaks %in% pretty(diff_breaks, n = 8))[-9], length(diff_colors))],
       pch = 15)

hist(diff_do, main = 'Histogram of ESA minius PLS', xlab = 'Biomass Difference')


dev.off()

####
#### cut pls and fia to be the same without IN and IL and so they all have 9293 points
####

match_fia <- which(
  outer(coors_dat_fia[, 1], coors_dat[, 1], "==") &
    outer(coors_dat_fia[, 2], coors_dat[, 2], "=="),arr.ind = T
)
coors_dat_fia_sm <- coors_dat_fia[match_fia[,1],]
plot(coors_dat_fia_sm[,1],coors_dat_fia_sm[,2])

match_pls <- which(
  outer(coors_dat_pls[, 1], coors_dat[, 1], "==") &
    outer(coors_dat_pls[, 2], coors_dat[, 2], "=="),arr.ind = T
)
coors_dat_pls_sm <- coors_dat_pls[match_pls[,1],]
plot(coors_dat_pls_sm[,1],coors_dat_pls_sm[,2])

#### Need all iterations to get error bar for gam map figure
nc_pls_iter <- nc_open(file.path('~/Downloads/PLS_biomass_agb_western_v1.0rc2.nc'))

x1 <- nc_pls_iter$dim$x$vals
y1 <- nc_pls_iter$dim$y$vals
data1 <- ncvar_get(nc_pls_iter,varid = c('Total'))

rownames(data1)  <- x1
colnames(data1)  <- y1

ras <- list()
for(ii in 1:250){
  rast1 <- raster(list(x=x1,y=y1,z=data1[,,ii]))
  ras[[ii]] <- rasterToPoints(rast1)
}

albers.df1 = as.data.frame(ras[[1]][,1:2])
colnames(albers.df1) = c('x', 'y')

coordinates(albers.df1) <- ~ x + y
proj4string(albers.df1) <- CRS('+init=epsg:3175')

lat.lon1 <- spTransform(albers.df1, CRS('+proj=longlat +ellps=WGS84'))
coors_dat1 <- as.matrix(data.frame(lat.lon1))

coords <- data.frame(x=refab$lon[chull(x= refab$lon,y=refab$lat)],y=refab$lat[chull(x= refab$lon,y=refab$lat)]) #convex hull coordinates
pts_in1 <- which(point.in.polygon(point.x = coors_dat_pls_sm[,1],point.y = coors_dat_pls_sm[,2],pol.x = coords$x,pol.y = coords$y)==1)

pls_hull <- do.call(cbind,lapply(ras,FUN=function(x) x[,3]))

breaks1 <- c(seq(0,50,10),seq(75,250,25),1000)
cutspls <- cut(coors_dat_pls_sm[pts_in1,3],breaks=breaks1)#cut(rowMeans(pls_hull[pts_in1,]),breaks=breaks1)
colors <- rev(terrain.colors(length(breaks1)-1))

plot(coors_dat_pls_sm[pts_in1,1],coors_dat_pls[pts_in1,2],col=colors[cutspls])

#regional mean PLS
pls_mean <- sum(coors_dat_pls_sm[,3])#mean(colSums(pls_hull[pts_in1,]))
pls_ci05 <- quantile(colSums(pls_hull[pts_in1,]),.05)
pls_ci95 <- quantile(colSums(pls_hull[pts_in1,]),.975)

#regional mean FIA no uncertainty yet
pts_in_fia_regional <- which(point.in.polygon(point.x = coors_dat_fia_sm[,1],point.y = coors_dat_fia_sm[,2],pol.x = coords$x,pol.y = coords$y)==1)

cutsfia <- cut(coors_dat_fia_sm[pts_in_fia_regional,3],breaks=c(seq(0,50,10),seq(75,250,25),435),labels=F)
plot(coors_dat_fia_sm[pts_in_fia_regional,1],coors_dat_fia_sm[pts_in_fia_regional,2],col=colors[cutsfia])


##difference maps

##fia - pls convex hull
breaks <- seq(-300,300,10)
colors <- cm.colors(length(breaks)-1)
diffcuts <- cut(coors_dat_fia_sm[pts_in_fia_regional,3]-coors_dat_pls_sm[pts_in1,3],
    breaks = breaks)
plot(coors_dat_fia_sm[pts_in_fia_regional,1],coors_dat_fia_sm[pts_in_fia_regional,2],col=colors[diffcuts])
plot(ch,add=T)

##fia - pls whole domain
breaks <- seq(-400,400,10)
colors <- colorRampPalette(c('blue','white','red'))(length(breaks)-1)
diffcuts <- cut(coors_dat_pls_sm[,3]-coors_dat_fia_sm[,3],
                breaks = breaks)
plot(coors_dat_fia_sm[,1],coors_dat_fia_sm[,2],col=colors[diffcuts],pch=19)
#plot(ch,add=T)
par(mfrow=c(2,3))
plot(coors_dat_fia_sm[,1],coors_dat_fia_sm[,2],col=colors[diffcuts],pch=19)
plot(ch1,add=T)
plot(coors_dat_fia_sm[,1],coors_dat_fia_sm[,2],col=colors[diffcuts],pch=19)
plot(ch2,add=T)
plot(coors_dat_fia_sm[,1],coors_dat_fia_sm[,2],col=colors[diffcuts],pch=19)
plot(ch3,add=T)

hist(coors_dat_fia_sm[pts_in1_fia,3] - coors_dat_pls_sm[pts_in1_pls,3],xlim=c(-350,350),col='gray')
abline(v=0,lwd=2,col='red')
hist(coors_dat_fia_sm[pts_in2_fia,3] - coors_dat_pls_sm[pts_in2_pls,3],xlim=c(-350,350),col='gray')
abline(v=0,lwd=2,col='red')
hist(coors_dat_fia_sm[pts_in3_fia,3] - coors_dat_pls_sm[pts_in3_pls,3],xlim=c(-350,350),col='gray')
abline(v=0,lwd=2,col='red')

##gam - pls whole domain
pred_biomass_gam <- pred_biomass_gam_list[[1]]
full.mat <-  cbind(coors_dat,as.vector(pred_biomass_gam$fit))
colnames(full.mat) <- c("x","y","pred_biomass")
breaks <- seq(-300,300,1)
colors <- colorRampPalette(c('blue','white','red'))(length(breaks)-1)
diffcuts <- cut(full.mat[,3]-coors_dat_pls_sm[,3],
                breaks = breaks)
plot(coors_dat_fia_sm[,1],coors_dat_fia_sm[,2],col=colors[diffcuts])
plot(ch,add=T)

sum(full.mat[pts_in_fia_regional,3]-coors_dat_pls_sm[pts_in_fia_regional,3])

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
coors_dat <- as.matrix(data.frame(lat.lon))


#### 
#### Plotting function for differnt time slices
#### 

plot_count <- 0

pdf('gam_maps_secondpass.pdf',height=12,width = 12,compress = T)

full.mat <- cbind(coors_dat_keep,as.vector(pred_biomass_gam))
library(spatstat)

ch <- convexhull.xy(x= refab$lon,y=refab$lat)

save_mat <- list()

load('clusters.Rdata')

alpha_all <- .5
title_cex <- 3
pt_cex <- 3
cluster_colors <- colors_tri <- viridis::viridis(3,begin = .75,end = 0,alpha = 1)#c('navy','royalblue','magenta3')


coords <- data.frame(x=refab$lon[chull(x= refab$lon,y=refab$lat)],y=refab$lat[chull(x= refab$lon,y=refab$lat)]) #convex hull coordinates
pts_in <- which(point.in.polygon(point.x = coors_dat[,1],point.y = coors_dat[,2],pol.x = coords$x,pol.y = coords$y)==1)

pred_biomass_gam_list <- list()
for (age_slice in seq(1, 100, 1)) {
  pred_data = cbind(coors_dat, rep(age_slice, nrow(coors_dat)))
  colnames(pred_data) <- c("lon", "lat", "variable")#c("lon","lat","age")
  pred_biomass_gam_list[[age_slice]] = (predict(
    b,
    newdata = as.data.frame(pred_data),
    type = 'response',
    se.fit = TRUE
  ))
  print(age_slice)
}

#save(pred_biomass_gam_list,file='pred_biomass_gam_list.Rdata')

pred_biomass_gam <- pred_biomass_gam_list[[age_slice]]
full_mat <- cbind(coors_dat,
                  do.call(cbind,lapply(pred_biomass_gam_list, FUN=function(x) x$fit)))
#### Save spatial predictions
colnames(full_mat)[-c(1,2)] <- paste0('AgeYBP',seq(100,10000,100))
save(full_mat,file='full_mat.Rdata')

#regional mean
space_mat <- do.call(rbind,lapply(pred_biomass_gam_list, FUN=function(x) quantile(x$fit,c(.25,.5,.75))))
space_mat_2 <- do.call(rbind,lapply(pred_biomass_gam_list, FUN=function(x) quantile(x$fit,c(.375,.625))))

#total mean
coords <- data.frame(x=refab$lon[chull(x= refab$lon,y=refab$lat)],y=refab$lat[chull(x= refab$lon,y=refab$lat)]) #convex hull coordinates
pts_in <- which(point.in.polygon(point.x = coors_dat[,1],point.y = coors_dat[,2],pol.x = coords$x,pol.y = coords$y)==1)
tot_mat <- do.call(rbind,lapply(pred_biomass_gam_list, FUN=function(x) c(sum(x$fit[pts_in]),sum(x$se[pts_in]))))

pdf('convexhull.pdf')
plot(refab$lon,refab$lat,ylab = 'Latitude',xlab='Longitude')
maps::map('state',add=T)
plot(ch,add=T)
dev.off()

source('Workflow_Code/utils/ciEnvelope.R')

pdf('regional_mean_and_total_ts.pdf',width=15,height=10)
layout(matrix(c(1,1,2,
                3,3,4),2,3,byrow=T))
par(cex.axis=2,oma=c(6,7,0,0))

plot(tot_mat[,1],
     xlim = c(100,0),
     ylim=c(3e5,9e5),
     ylab = NA,
     xlab = NA,
     xaxt= 'n',typ='l',lty=2,las=2)
ciEnvelope(x = 1:100,ylo = tot_mat[,1]-tot_mat[,2],yhi = tot_mat[,1]+tot_mat[,2],col=adjustcolor('gray',alpha.f = .5))
points(1:100,tot_mat[,1]-tot_mat[,2],typ='l',lty=1,lwd=2)
points(1:100,tot_mat[,1]+tot_mat[,2],typ='l',lty=1,lwd=2)

points(1,pls_mean,pch=19,col='darkgreen',cex=.5)

mtext('Total Regional Biomass (Mg)', side=2, outer = F,cex=2,line = 7)
axis(side=1,at = seq(0,100,5),labels = seq(0,10000,500),las=2)
axis(side=4,at = seq(3e5,9e5,1e5),labels = seq(3e5,9e5,1e5),las=2)

plot.new()
legend('center',c('Total SE','Total Median'),
       col=c('black','black'),
       lty=c(1,2),
       lwd=3,cex=2)

matplot(t(refab[, 2:101]), typ = 'l',lty=1,lwd=1.5,
        col = cluster_colors[clusters@cluster[refab$name]],
        xlim = c(100,0),
        ylim=c(0,250),
        las = 2,
        ylab = NA,
        xlab = NA,
        xaxt= 'n')
mtext('Years Before Present', side=1, outer = F,line=8,cex=2)
mtext('Biomass (Mg/ha)', side=2, outer = F,cex=2,line = 5)
axis(side=1,at = seq(0,100,5),labels = seq(0,10000,500),las=2)
axis(side=4,at = seq(0,250,50),labels = seq(0,250,50),las=2)

ciEnvelope(x = 1:100,ylo = space_mat[,1],yhi = space_mat[,3],col=adjustcolor('gray',alpha.f = .5))
ciEnvelope(x = 1:100,ylo = space_mat_2[,1],yhi = space_mat_2[,2],col=adjustcolor('gray',alpha.f = .5))
points(1:100,space_mat[,1],typ='l',lty=1,lwd=2)
points(1:100,space_mat[,3],typ='l',lty=1,lwd=2)

points(1:100,space_mat_2[,1],typ='l',lty=2,lwd=2)
points(1:100,space_mat_2[,2],typ='l',lty=2,lwd=2)

points(1:100,space_mat[,2],typ='l',lty=3,lwd=2)

plot.new()
legend('center',c('Regional 50% Quantiles','Regional 25% Quantiles','Regional Median','Cluster 1','Cluster 2','Cluster 3'),
       col=c('black','black','black',(cluster_colors)[c(3,1,2)]),
       lty=c(1,2,3,1,1,1),
       lwd=3,cex=2)
dev.off()


#####
##### this part calculates the regional biomass for each cluster's convex hull
#####

m_1 <- m_2 <- m_3 <- v_1 <- v_2 <- v_3 <- v1_1 <- v1_2 <- v1_3 <- numeric(100)

pt_data <- refab_melt[refab_melt$variable==1,] #1 was t but should be the same no matter what time so taking it out of the loop for now
#Calculating total AGB in the convex hull of each cluster
#Might not be the best because they are all different sizes
#note: chull give you the indices of the corners of the convex hull from whatever dataframe you gave
chull1 <- chull(x=refab[which(clusters@cluster[pt_data$name]==1),c('lon')],y=refab[which(clusters@cluster[pt_data$name]==1),c('lat')])
chull2 <- chull(x=refab[which(clusters@cluster[pt_data$name]==2),c('lon')],y=refab[which(clusters@cluster[pt_data$name]==2),c('lat')])
chull3 <- chull(x=refab[which(clusters@cluster[pt_data$name]==3),c('lon')],y=refab[which(clusters@cluster[pt_data$name]==3),c('lat')])

ch1 <- convexhull.xy(x=refab[which(clusters@cluster[pt_data$name]==1),c('lon')],y=refab[which(clusters@cluster[pt_data$name]==1),c('lat')])
ch2 <- convexhull.xy(x=refab[which(clusters@cluster[pt_data$name]==2),c('lon')],y=refab[which(clusters@cluster[pt_data$name]==2),c('lat')])
ch3 <- convexhull.xy(x=refab[which(clusters@cluster[pt_data$name]==3),c('lon')],y=refab[which(clusters@cluster[pt_data$name]==3),c('lat')])

c1 <- refab[which(clusters@cluster[pt_data$name]==1),c('lon','lat')][chull1,]
c2 <- refab[which(clusters@cluster[pt_data$name]==2),c('lon','lat')][chull2,]
c3 <- refab[which(clusters@cluster[pt_data$name]==3),c('lon','lat')][chull3,]

pts_in1 <- which(point.in.polygon(point.x = full_mat[,1],point.y = full_mat[,2],pol.x = c1[,1],pol.y = c1[,2])==1)
pts_in2 <- which(point.in.polygon(point.x = full_mat[,1],point.y = full_mat[,2],pol.x = c2[,1],pol.y = c2[,2])==1)
pts_in3 <- which(point.in.polygon(point.x = full_mat[,1],point.y = full_mat[,2],pol.x = c3[,1],pol.y = c3[,2])==1)

## adding another part to calculate the difference between pls and fia
## points to take from data products for pls and fia in 3 convex hulls
pts_in1_pls <- which(point.in.polygon(point.x = coors_dat_pls_sm[,1],point.y = coors_dat_pls_sm[,2],pol.x = c1[,1],pol.y = c1[,2])==1)
pts_in2_pls <- which(point.in.polygon(point.x = coors_dat_pls_sm[,1],point.y = coors_dat_pls_sm[,2],pol.x = c2[,1],pol.y = c2[,2])==1)
pts_in3_pls <- which(point.in.polygon(point.x = coors_dat_pls_sm[,1],point.y = coors_dat_pls_sm[,2],pol.x = c3[,1],pol.y = c3[,2])==1)

coors_dat_fia_sm <- coors_dat_esa[which(!is.na(ove)), ]

pts_in1_fia <- which(point.in.polygon(point.x = coors_dat_fia_sm[,1],point.y = coors_dat_fia_sm[,2],pol.x = c1[,1],pol.y = c1[,2])==1)
pts_in2_fia <- which(point.in.polygon(point.x = coors_dat_fia_sm[,1],point.y = coors_dat_fia_sm[,2],pol.x = c2[,1],pol.y = c2[,2])==1)
pts_in3_fia <- which(point.in.polygon(point.x = coors_dat_fia_sm[,1],point.y = coors_dat_fia_sm[,2],pol.x = c3[,1],pol.y = c3[,2])==1)

set_mod_diff1 <- sum(coors_dat_fia_sm[pts_in1_fia,3]) - sum(coors_dat_pls_sm[pts_in1_pls,3])
set_mod_diff2 <- sum(coors_dat_fia_sm[pts_in2_fia,3]) - sum(coors_dat_pls_sm[pts_in2_pls,3])
set_mod_diff3 <- sum(coors_dat_fia_sm[pts_in3_fia,3]) - sum(coors_dat_pls_sm[pts_in3_pls,3])

sum(set_mod_diff1,set_mod_diff2,set_mod_diff3)
fia_mean <- sum(coors_dat_fia_sm[pts_in,3])

for(t in 1:100){
  #Calculates sum and var by refab estimate points. But this doesn't tell the whole landscape story...
  m_1[t] <- sum(full_mat[pts_in1,2+t])#/length(pts_in1)
  m_2[t] <- sum(full_mat[pts_in2,2+t])#/length(pts_in2)
  m_3[t] <- sum(full_mat[pts_in3,2+t])#/length(pts_in3)
  
  v_1[t] <- quantile(full_mat[pts_in1,2+t],.25)
  v_2[t] <- quantile(full_mat[pts_in2,2+t],.25)
  v_3[t] <- quantile(full_mat[pts_in3,2+t],.25)
  
  v1_1[t] <- quantile(full_mat[pts_in1,2+t],.75)
  v1_2[t] <- quantile(full_mat[pts_in2,2+t],.75)
  v1_3[t] <- quantile(full_mat[pts_in3,2+t],.75)
}

mean_clusts <- rbind(m_1,m_2,m_3)
sd_clusts <- rbind(v_1,v_2,v_3)
sd_clusts1 <- rbind(v1_1,v1_2,v1_3)

save(sd_clusts,file='sd_clusts.Rdata')

#####
##### PLOT
#####

pdf(paste0(Sys.Date(),'gam_maps_with_ts.pdf'),compress = T,height=7,width = 12)

#layout(matrix(c(1,2,9,10,3,4,11,12,5,6,13,14,7,8,15,16,17,18,19,20),4,5))
#layout(matrix(c(1,2,3,10,4,5,6,10,7,8,9,11),3,4,byrow=T))

layout(matrix(c(1,1,2,2,3,3,4,4,5,
                1,1,2,2,3,3,4,4,5,
                rep(6,9),
                rep(6,9),
                rep(7,9),
                rep(7,9)),6,9,byrow=T))

par(oma=c(6,10,2,6),mar=rep(0,4))

for(age_slice in c(100,80,50,10)){ #rev(seq(10,80,10))
print(age_slice)
  par(mar=c(0,0,0,0))
  pred_biomass_gam <- pred_biomass_gam_list[[age_slice]]
  
  full.mat <- save_mat[[age_slice]] <- cbind(coors_dat,as.vector(pred_biomass_gam$fit))
  
  colnames(full.mat) <- c("x","y","pred_biomass")

  # breaks <-  c(seq(0,50,10),seq(75,250,25),435)
  # colors <- rev(terrain.colors(length(breaks)-1))
  breaks <-  c(0,25, 50,75,100,125,150,175,200,225,250,300,800)#c(0, seq(25, 100, 25), seq(150, 300, 50), max(y$Tot_Biomass))
  #breaks <-  c(0,50,100,150,200,250, max(y$Tot_Biomass))#c(0, seq(25, 100, 25), seq(150, 300, 50), max(y$Tot_Biomass))
  colors <- rev(terrain.colors(length(breaks))[-4])
  breaklabels <- signif(breaks[2:length(breaks)],digits = 2)
  
  data_binned <-  cut(pred_biomass_gam$fit, breaks, include.lowest = TRUE, labels = FALSE)
  
  legendName <- c("Biomass (Mg/ha)")#paste0("Biomass at Age = ",age_slice, " BP"
  #breaklabels <- apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1,  function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })
  
  coords <- data.frame(x=refab$lon[chull(x= refab$lon,y=refab$lat)],y=refab$lat[chull(x= refab$lon,y=refab$lat)]) #convex hull coordinates
  pts_in <- which(point.in.polygon(point.x = full.mat[,1],point.y = full.mat[,2],pol.x = coords$x,pol.y = coords$y)==1)
  
  biomass_sum <- round(sum(full.mat[pts_in,3]))
  
  plot(full.mat[pts_in,1],full.mat[pts_in,2],
       col=adjustcolor(colors[data_binned[pts_in]],
                       alpha.f = alpha_all),
       pch=19,
       xaxt='n',
       yaxt='n')
  maps::map('state',add=T)
  plot(ch,add=T)
  mtext(paste(age_slice/10,'cal ka BP'),side = 3)

  
  #### REFAB POINTS 
  #points_get <- x.meta[x.meta$age_bacon<(age_slice+100)&x.meta$age_bacon>(age_slice-100)&x.meta$site.name%in%dataID_use$name,c('lat','long')]
  pt_data <- refab_melt[refab_melt$variable==age_slice,]
  data_binned_pts <- cut(pt_data[,'value'],breaks=breaks,labels=F)
  points(pt_data[,2],pt_data[,1],pch=(21:23)[clusters@cluster[pt_data$name]],bg=cluster_colors[clusters@cluster[pt_data$name]],col=NA,cex=1)#,bg=colors[data_binned_pts]
  
  #### legend
  if(age_slice==10) {
    plot.new()
    legend(
      'center',
      c(as.character(breaklabels)[c(1,3,5,7,9)],'>300'),
      pch = 22,
      pt.bg = colors[c(1,3,5,7,9,11)],
      cex = 1.5,
      pt.cex = 3.5,
      col=NA,
      title = 'Biomass \n (Mg/ha)',box.col = 'white'
    )
  }
  
  biomass_sum1 <- round(sum(full.mat[pts_in1,3]))
  biomass_sum2 <- round(sum(full.mat[pts_in2,3]))
  biomass_sum3 <- round(sum(full.mat[pts_in3,3]))
  
  #sum(biomass_sum1,biomass_sum2,biomass_sum3) #doesn't add up because convex hulls are overlapping for clusters
  if(FALSE){
  legend('topright',c(paste('Total AGB =',signif(biomass_sum,digits = 3),'Mg'),
                          paste('Cluster 1 AGB =',signif(biomass_sum1,digits = 3),'Mg'),
                                paste('Cluster 2 AGB =',signif(biomass_sum2,digits = 3),'Mg'),
                                      paste('Cluster 3 AGB =',signif(biomass_sum3,digits = 3),'Mg')))
  }
  
  
  # #### REFAB POINTS 
  # #points_get <- x.meta[x.meta$age_bacon<(age_slice+100)&x.meta$age_bacon>(age_slice-100)&x.meta$site.name%in%dataID_use$name,c('lat','long')]
  # pt_data <- refab_melt[refab_melt$variable==age_slice,]
  # data_binned_pts <- cut(pt_data[,'value'],breaks=breaks,labels=F)
  # points(pt_data[,2],pt_data[,1],pch=21,col=cluster_colors[clusters@cluster[pt_data$name]],bg=colors[data_binned_pts],cex=pt_cex)
  # 
  # points_get <- x.meta[x.meta$age_bacon<(age_slice+100)&x.meta$age_bacon>(age_slice-100)&x.meta$site.name%in%dataID_use$name&x.meta$lat%in%cut_preds1[,1],c('lat','long')]
  # points(points_get[,2],points_get[,1],pch=8,col='blue',cex=1.5)
  # 
  # points_cols <- cut_preds1[which(cut_preds1[,'age']==age_slice),c('lon','lat','agwb')]
  # 
  # data_binned_pts <-  cut(points_cols[,'agwb'], breaks, include.lowest = TRUE, labels = FALSE)
  # points(points_cols[,1],points_cols[,2],pch=21,bg=colors[data_binned_pts],col='black',cex=1)
  # 
  # legend('topright',c('avg bacon age +- 100 years',
  #                     'ReFAB avg pt est'),
  #        pch = c(8,21),
  #        col = c('blue','black'),
  #        pt.bg = c(NA,colors[5]))
  
  #Calculates mean and var by refab estimate points. But this doesn't tell the whole landscape story...
  m1 <- mean(pt_data[which(clusters@cluster[pt_data$name]==1),c('value')])
  m2 <- mean(pt_data[which(clusters@cluster[pt_data$name]==2),c('value')])
  m3 <- mean(pt_data[which(clusters@cluster[pt_data$name]==3),c('value')])
  
  v1 <- sd(pt_data[which(clusters@cluster[pt_data$name]==1),c('value')])
  v2 <- sd(pt_data[which(clusters@cluster[pt_data$name]==2),c('value')])
  v3 <- sd(pt_data[which(clusters@cluster[pt_data$name]==3),c('value')])
  
  # par(mar=c(0,0,0,0)) #removed 12/2
  # plot.new()
  # legend('center',c(paste('Total AGB =',signif(biomass_sum,digits = 3),'Mg'),
  #                     paste('Cluster 1 AGB Mean =',signif(m1,digits = 3),'Mg/ha,','SD =',signif(v1,digits = 3),'Mg/ha,'),
  #                     paste('Cluster 2 AGB Mean =',signif(m2,digits = 3),'Mg/ha,','SD =',signif(v2,digits = 3),'Mg/ha,'),
  #                     paste('Cluster 3 AGB Mean =',signif(m3,digits = 3),'Mg/ha,','SD =',signif(v3,digits = 3),'Mg/ha,')),
  # 
  #        text.col = c('black',cluster_colors[1:3]),
  #        pch = c(NA,21:23),
  #        bg='white',cex=.5)
}

pls_hull <- do.call(cbind,lapply(ras,FUN=function(x) x[,3]))

pls_mean <- mean(colSums(pls_hull[pts_in,]))
pls_ci05 <- quantile(colSums(pls_hull[pts_in,]),.05)
pls_ci95 <- quantile(colSums(pls_hull[pts_in,]),.975)
### regional mean time series
plot(tot_mat[,1],
     xlim = c(100,0),
     ylim=c(2e5,7e5),
     ylab = NA,
     xlab = NA,
     xaxt= 'n',yaxt='n',typ='l',lty=1,lwd=2,las=2)
abline(v=c(100,80,50,10),lty=1,col='lightgray')
abline(h=seq(2e5,7e5,1e5)[-c(1,6)],lty=1,col='lightgray')

axis(side=2,at=seq(200000,700000,100000),labels=prettyNum(seq(200000,700000,100000),big.mark = ','),las=2)

ciEnvelope(x = 1:100,ylo = tot_mat[,1]-tot_mat[,2],yhi = tot_mat[,1]+tot_mat[,2],col=adjustcolor('#3b4c41',alpha.f = .5))
points(1:100,tot_mat[,1]-tot_mat[,2],typ='l',lty=2,lwd=2)
points(1:100,tot_mat[,1]+tot_mat[,2],typ='l',lty=2,lwd=2)
#rect(xleft = .5,ybottom = -18000000,xright = -5,ytop = 18000000,col=adjustcolor('gray',alpha.f = .5))


#segments(x0 = 1,x1=1, y0= pls_ci05,y1=pls_ci95,col='darkred',lwd=4)
points(.5,pls_mean,pch=19,col='darkred',cex=2)
text(-2,pls_mean,'1800s')
#abline(h=c(500000,600000,700000),col='gray')
points(0,fia_mean,pch=19,col='red',cex=2)
text(-2.25,fia_mean,'2000s')

options(scipen=999) #gets rid of scientific notation
mtext('Total Regional\n Biomass (Mg)', side=2, outer = F,line = 5)



par(mar=c(0,0,0,0))
plot(99:1, diff(rev(m_1)),pch=21,bg=cluster_colors[1],ylim=c(-15000,8000),xlim=c(100,0),
     ylab = NA,
     xlab = NA,
     xaxt= 'n',
     yaxt= 'n',
     las=2,cex=1.5,col=NA)
axis(2,at=c(-8000,-4000,0,4000,8000),labels = prettyNum(c(-8000,-4000,0,4000,8000),big.mark = ','),las=2)

points(99:1,diff(rev(m_2)),pch=22,bg=cluster_colors[2],cex=1.5,col=NA)
points(99:1,diff(rev(m_3)),pch=23,bg=cluster_colors[3],cex=1.5,col=NA)
#abline(h=c(-2000,2000),lty=1,col='gray')
abline(h=c(0),lty=1)
abline(v=(c(100,80,50,10)),lty=1,col='lightgray')
mtext('100 Year AGB \n Differences (Mg)', side=2, outer = F,line = 5)

axis.break(2,-9000,style="slash") 


par(new=T)

plot(-1,set_mod_diff1,xlim=c(100,0),ylim=c(-150000,100000),xaxt='n',yaxt='n',ylab=NA,xlab=NA,col='white')
axis.break(4,-75000,style="slash") 

#rect(xleft = .5,ybottom = -18000000,xright = -5,ytop = 18000000,col=adjustcolor('gray',alpha.f = .5))
axis(side=4,at =seq(-140000,-70000,30000),labels = prettyNum(seq(-140000,-70000,30000),big.mark = ','),las=2)
points(-1,set_mod_diff1,pch=21,bg=cluster_colors[1],cex=2.5,col=NA)
points(-1,set_mod_diff2,pch=22,bg=cluster_colors[2],cex=2.5,col=NA)
points(-1,set_mod_diff3,pch=23,bg=cluster_colors[3],cex=2.5,col=NA)
#abline(h=c(set_mod_diff1,set_mod_diff2,set_mod_diff3),col=cluster_colors,lty=1)

legend('bottomleft',c('West','Central','East'),pch=21:23,pt.bg = (cluster_colors),col = NA,ncol = 3,cex=1.5,pt.cex = 1.75,bg='white')

# axis(side=1,at = seq(0,100,5),labels = rev(seq(0,10000,500)),las=2)
# axis(side=4,at = seq(-7000,7000,1000),labels = seq(-7000,7000,1000),las=2)

levs <- c('forest','PINUSX','QUERCUS','TSUGAX','conifer','FAGUS')
levs_cols <- RColorBrewer::brewer.pal(n=length(levs),'Paired')

# for(i in 1:3){
#   linesmat <- cast(cm4[cm4$clusters==i,],variable~bins)
#   
#   matplot(seq(1,9,1),apply(
#     linesmat[, 11:2],
#     1,
#     FUN = function(x)
#       diff(x / max(x))
#   ),
#   lty=1,
#   typ = 'b',
#   lwd=2,
#   col = levs_cols[which(levs%in%linesmat[,1])],pch=19,
#   ylab = NA,xaxt='n',las=2,ylim=c(-.45,.45),xlim=c(0,10))
#   
#   abline(h = 0)
# 
#   box(which = 'plot',col=cluster_colors[i],lwd=4)
#   
#   if(i == 1){
#     legend('topleft',legend=
#              levs,col=levs_cols,pch=19,ncol = 2)
#   }
#   if(i == 2){
#     mtext(text = 'Millenial Normalized Pollen Proportion Change',side = 2,line = 5,outer = F)
#   }
# }



axis(side=1,at = seq(0,100,10),labels = seq(0,10,length.out=11),las=1)
# axis(side=4,at = seq(300000,900000,100000),labels = seq(300000,900000,100000),las=2)

# legend('bottomright',c('Regional AGB SE','Regional AGB Median','PLS AGB 95% CIs'),
#        col=c('black','black','green'),
#        lty=c(2,1,1),
#        lwd=3)

mtext('Age (cal ka BP)', side=1, outer = F,line = 3)

dev.off()



#####
##### END PLOT
#####

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
legend('center',legend=breaklabels_se,col=colors_se,pch=19,cex=2,title = 'SE')

dev.off()
### how much to cut off?




#### Old difference maps from Jason march 2020 visit

layout(matrix(c(1,2,3,0),nrow=2,ncol=2,byrow=T))
par(oma=c(2,2,0,0),mar=c(0,0,4,0))

age_slice = 1

pred_biomass_gam <- pred_biomass_gam_list[[age_slice]]

full.mat <- save_mat[[age_slice]] <- cbind(coors_dat,as.vector(pred_biomass_gam$fit))
colnames(full.mat) <- c("x","y","pred_biomass")

breaks <-  c(seq(0,50,10),seq(75,250,25),435)
colors <- rev(viridis(length(breaks)-1))
data_binned <-  cut(pred_biomass_gam$fit, breaks, include.lowest = TRUE, labels = FALSE)

legendName <- c("Biomass (Mg/ha)")#paste0("Biomass at Age = ",age_slice, " BP"
breaklabels <- apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1,  function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })

coords <- data.frame(x=refab$lon[chull(x= refab$lon,y=refab$lat)],y=refab$lat[chull(x= refab$lon,y=refab$lat)]) #convex hull coordinates
pts_in <- which(point.in.polygon(point.x = full.mat[,1],point.y = full.mat[,2],pol.x = coords$x,pol.y = coords$y)==1)

biomass_sum <- round(sum(full.mat[pts_in,3]))

plot(full.mat[pts_in,1],full.mat[pts_in,2],
     col=adjustcolor(colors[data_binned[pts_in]],
                     alpha.f = alpha_all),
     pch=19,
     xaxt='n',
     yaxt='n')
maps::map('state',add=T)
plot(ch,add=T)
title(paste('Biomass',age_slice*100,'YBP'),cex.main=title_cex)

coords <- data.frame(x=refab$lon[chull(x= refab$lon,y=refab$lat)],y=refab$lat[chull(x= refab$lon,y=refab$lat)]) #convex hull coordinates
pts_in1 <- which(point.in.polygon(point.x = coors_dat1[,1],point.y = coors_dat1[,2],pol.x = coords$x,pol.y = coords$y)==1)

pls_hull <- do.call(cbind,lapply(ras,FUN=function(x) x[pts_in1,3]))

cutspls <- cut(rowMeans(pls_hull),breaks=c(seq(0,50,10),seq(75,250,25),435))
colors <- rev(viridis(length(breaks)-1))

plot(coors_dat1[pts_in1,1],coors_dat1[pts_in1,2],col=colors[cutspls])


pls_hull <- do.call(cbind,lapply(ras,FUN=function(x) x[,3]))

diff_mat1 <- save_mat[[1]][,3] - rowMeans(pls_hull)[1:nrow(save_mat[[1]])]#save_mat[[80]][,3] 

sum(diff_mat1)

diff_breaks <-  c(seq(-200,-50,50),seq(-30,30,20),seq(50,200,50))
diff_colors <- (colorRampPalette(c('darkblue','blue','white','red','darkred'))(length(diff_breaks)-1))
diff_data_binned <-  cut(diff_mat1, diff_breaks, include.lowest = TRUE, labels = FALSE)


plot(full.mat[pts_in,1],full.mat[pts_in,2],
     col=adjustcolor(diff_colors[diff_data_binned[pts_in]],
                     alpha.f = alpha_all),
     pch=19,
     xaxt='n',
     yaxt='n')
maps::map('state',add=T)
plot(ch,add=T)
title(paste('Biomass Difference 100 YBP - PLS'),cex.main=title_cex-.5)

diff_pt_data <- refab_melt[refab_melt$variable==50,'value'] - refab_melt[refab_melt$variable==80,'value']
diff_data_binned_pts <- cut(diff_pt_data,breaks=diff_breaks,labels=F)
points(pt_data[,2],pt_data[,1],
       pch=(21:23)[clusters@cluster[pt_data$name]],
       col=cluster_colors[clusters@cluster[pt_data$name]],
       bg=diff_colors[diff_data_binned_pts],cex=pt_cex)
















diff_mat2 <- save_mat[[10]][,3] - save_mat[[50]][,3]

diff_data_binned2 <-  cut(diff_mat2, diff_breaks, include.lowest = TRUE, labels = FALSE)

plot(full.mat[pts_in,1],full.mat[pts_in,2],
     col=adjustcolor(diff_colors[diff_data_binned2[pts_in]],
                     alpha.f = alpha_all),
     pch=19,
     xaxt='n',
     yaxt='n')
maps::map('state',add=T)
plot(ch,add=T)
title(paste('Biomass Difference 1000 - 5000 YBP'),cex.main=title_cex-.5)

diff_pt_data <- refab_melt[refab_melt$variable==10,'value'] - refab_melt[refab_melt$variable==50,'value']
diff_data_binned_pts <- cut(diff_pt_data,breaks=diff_breaks,labels=F)
points(pt_data[,2],pt_data[,1],
       pch=(21:23)[clusters@cluster[pt_data$name]],
       col=cluster_colors[clusters@cluster[pt_data$name]],
       bg=diff_colors[diff_data_binned_pts],cex=pt_cex)

diff_legendName <- c("Biomass Difference (Mg/ha)")#paste0("Biomass at Age = ",age_slice, " BP"
diff_breaklabels <- apply(cbind(diff_breaks[1:(length(diff_breaks)-1)], diff_breaks[2:length(diff_breaks)]), 1,  function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })

plot.new()
legend('center',diff_breaklabels,pch=19,col=diff_colors,cex=1.5,title = 'Biomass Difference (Mg/ha)')


