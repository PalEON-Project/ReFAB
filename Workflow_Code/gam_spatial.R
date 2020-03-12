
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
<<<<<<< HEAD
library(maps)
=======
library(raster)
>>>>>>> d76c2fe... little update

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

<<<<<<< HEAD
par(mfrow=c(1,2))
plot(all.preds1[,2],all.preds1[,1])
map('state',add=T)

cut_preds <- all.preds1[-which(all.preds1[,1]< 42.6 | all.preds1[,2]> -84.67906),]
cut_preds1 <- cut_preds[-which(cut_preds[,1]< 46 & cut_preds[,2]> -86),]
plot(cut_preds1[,2],cut_preds1[,1])
map('state',add=T)
=======
>>>>>>> d76c2fe... little update

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
<<<<<<< HEAD
b <- gam(log(value) ~ te(lon, lat, variable, d = c(2,1),
                        bs = c("tp","cr"), k=77),
         data = as.data.frame(refab_melt))
=======
tictoc::tic()
b <- gam(value ~ te(lon, lat, variable, d = c(2,1),
                        bs = c("tp","cr"), k=70), #can choose k up to 77
         data = as.data.frame(refab_melt),
         control = gam.control(nthreads = 4))
tictoc::toc()
<<<<<<< HEAD
>>>>>>> daaa8bc... fixing the graphics to make something nice for publication

=======
}
>>>>>>> d6f02b0... making better plots
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

#### 
#### Plotting function for differnt time slices
#### 
plot_count <- 0

<<<<<<< HEAD
pdf('gam_maps_secondpass.pdf',height=12,width = 12,compress = T)
par(mfrow=c(2,2))
for(age_slice in rev(seq(1000,10000,1000))){

<<<<<<< HEAD
  pred_data = cbind(coors_dat_keep,rep(age_slice,nrow(coors_dat_keep)))
  colnames(pred_data)<- c("lon","lat","age")
=======
  pred_data = cbind(coors_dat,rep(1000,nrow(coors_dat)))
  colnames(pred_data)<- c("lon","lat","variable")#c("lon","lat","age")
>>>>>>> d76c2fe... little update
  pred_biomass_gam = exp(predict(b,newdata = as.data.frame(pred_data)))
  
  full.mat <- cbind(coors_dat_keep,as.vector(pred_biomass_gam))
=======
library(spatstat)

ch <- convexhull.xy(x= refab$lon,y=refab$lat)

save_mat <- list()

load('clusters.Rdata')


alpha_all <- .5
title_cex <- 3
pt_cex <- 3
cluster_colors <- wes_palette('Rushmore1',5)[c(5,3,4)]#viridis::plasma(3)

pred_biomass_gam_list <- list()
for(age_slice in c(80,50,10)){
pred_data = cbind(coors_dat,rep(age_slice,nrow(coors_dat)))
colnames(pred_data)<- c("lon","lat","variable")#c("lon","lat","age")
pred_biomass_gam_list[[age_slice]] = (predict(b,newdata = as.data.frame(pred_data),type='response',se.fit=TRUE))
}

pdf('gam_maps_with_diff.pdf',compress = T,height=13,width = 20)

#layout(matrix(c(1,2,9,10,3,4,11,12,5,6,13,14,7,8,15,16,17,18,19,20),4,5))

library(wesanderson)

#layout(matrix(c(1,2,3,10,4,5,6,10,7,8,9,11),3,4,byrow=T))

layout(matrix(c(1,1,2,2,3,3,6,0,4,4,5,5,0,7),2,7,byrow=T))

par(oma=c(2,2,0,0),mar=c(0,0,4,0))

for(age_slice in c(80,50,10)){ #rev(seq(10,80,10))
print(age_slice)
  
<<<<<<< HEAD
  full.mat <- cbind(coors_dat,as.vector(pred_biomass_gam$fit))
>>>>>>> daaa8bc... fixing the graphics to make something nice for publication
=======
  pred_biomass_gam <- pred_biomass_gam_list[[age_slice]]
  
  full.mat <- save_mat[[age_slice]] <- cbind(coors_dat,as.vector(pred_biomass_gam$fit))
>>>>>>> 917342b... messing with colors
  colnames(full.mat) <- c("x","y","pred_biomass")

  breaks <-  c(seq(0,50,10),seq(75,250,25),435)
  colors <- rev(terrain.colors(length(breaks)-1))
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
  
  
  #Calculating total AGB in the convex hull of each cluster
  #Might not be the best because they are all different sizes
  chull1 <- chull(x=refab[which(clusters@cluster[pt_data$name]==1),c('lon')],y=refab[which(clusters@cluster[pt_data$name]==1),c('lat')])
  chull2 <- chull(x=refab[which(clusters@cluster[pt_data$name]==2),c('lon')],y=refab[which(clusters@cluster[pt_data$name]==2),c('lat')])
  chull3 <- chull(x=refab[which(clusters@cluster[pt_data$name]==3),c('lon')],y=refab[which(clusters@cluster[pt_data$name]==3),c('lat')])
  
  c1 <- refab[which(clusters@cluster[pt_data$name]==1),c('lon','lat')][chull1,]
  c2 <- refab[which(clusters@cluster[pt_data$name]==2),c('lon','lat')][chull2,]
  c3 <- refab[which(clusters@cluster[pt_data$name]==3),c('lon','lat')][chull3,]
  
  pts_in1 <- which(point.in.polygon(point.x = full.mat[,1],point.y = full.mat[,2],pol.x = c1[,1],pol.y = c1[,2])==1)
  pts_in2 <- which(point.in.polygon(point.x = full.mat[,1],point.y = full.mat[,2],pol.x = c2[,1],pol.y = c2[,2])==1)
  pts_in3 <- which(point.in.polygon(point.x = full.mat[,1],point.y = full.mat[,2],pol.x = c3[,1],pol.y = c3[,2])==1)
  
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
  
  
  #### REFAB POINTS 
  #points_get <- x.meta[x.meta$age_bacon<(age_slice+100)&x.meta$age_bacon>(age_slice-100)&x.meta$site.name%in%dataID_use$name,c('lat','long')]
  pt_data <- refab_melt[refab_melt$variable==age_slice,]
  data_binned_pts <- cut(pt_data[,'value'],breaks=breaks,labels=F)
  points(pt_data[,2],pt_data[,1],pch=21,col=cluster_colors[clusters@cluster[pt_data$name]],bg=colors[data_binned_pts],cex=pt_cex)
  
<<<<<<< HEAD
<<<<<<< HEAD
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
=======
>>>>>>> daaa8bc... fixing the graphics to make something nice for publication
=======
  #Calculates mean and var by refab estimate points. But this doesn't tell the whole landscape story...
  m1 <- mean(pt_data[which(clusters@cluster[pt_data$name]==1),c('value')])
  m2 <- mean(pt_data[which(clusters@cluster[pt_data$name]==2),c('value')])
  m3 <- mean(pt_data[which(clusters@cluster[pt_data$name]==3),c('value')])
>>>>>>> d6f02b0... making better plots
  
  v1 <- sd(pt_data[which(clusters@cluster[pt_data$name]==1),c('value')])
  v2 <- sd(pt_data[which(clusters@cluster[pt_data$name]==2),c('value')])
  v3 <- sd(pt_data[which(clusters@cluster[pt_data$name]==3),c('value')])
  
  legend('topright',c(paste('Total AGB =',signif(biomass_sum,digits = 3),'Mg'),
                      paste('Cluster 1 AGB Mean =',signif(m1,digits = 3),'Mg/ha,','SD =',signif(v1,digits = 3),'Mg/ha,'),
                      paste('Cluster 2 AGB Mean =',signif(m2,digits = 3),'Mg/ha,','SD =',signif(v2,digits = 3),'Mg/ha,'),
                      paste('Cluster 3 AGB Mean =',signif(m3,digits = 3),'Mg/ha,','SD =',signif(v3,digits = 3),'Mg/ha,')),
         
         text.col = c('black',cluster_colors[1:3]),
         cex = 1,
         bg='white')
  
 
}

diff_mat1 <- save_mat[[50]][,3] - save_mat[[80]][,3] 

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
title(paste('Biomass Difference 5000 - 8000 YBP'),cex.main=title_cex-.5)

diff_pt_data <- refab_melt[refab_melt$variable==50,'value'] - refab_melt[refab_melt$variable==80,'value']
diff_data_binned_pts <- cut(diff_pt_data,breaks=diff_breaks,labels=F)
points(pt_data[,2],pt_data[,1],pch=21,col=wes_palette('Rushmore1',n=5)[c(3,5,4)][clusters@cluster[pt_data$name]],
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
points(pt_data[,2],pt_data[,1],pch=21,col=wes_palette('Rushmore1',n=5)[c(3,5,4)][clusters@cluster[pt_data$name]],
       bg=diff_colors[diff_data_binned_pts],cex=pt_cex)


plot.new()
legend('center',breaklabels,pch=19,col=colors,cex=2,title = 'Biomass (Mg/ha)')

diff_legendName <- c("Biomass Difference (Mg/ha)")#paste0("Biomass at Age = ",age_slice, " BP"
diff_breaklabels <- apply(cbind(diff_breaks[1:(length(diff_breaks)-1)], diff_breaks[2:length(diff_breaks)]), 1,  function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })

plot.new()
legend('center',diff_breaklabels,pch=19,col=diff_colors,cex=1.5,title = 'Biomass Difference (Mg/ha)')

dev.off()

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
