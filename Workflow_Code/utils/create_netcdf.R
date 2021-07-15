
## This script is used to create more widely useable data products from the ReFAB project

####
#### Calibration Dataset ####
####

load('final_datasets/cast.x.Rdata')
load('biomass_draws_v3.0.Rdata')

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
                    all.pollen.taxa.names = colnames(cast.x)[-c(1:4)],
                    prairie.include = T,bigwoods.include=F, other.herbs.include = T,
                    other.trees.include = T, drop.taxa = drop.taxa,
                    PFT.do = F)

all_biom_draws <- do.call(rbind,biomass_draws)
colnames(all_biom_draws) <- paste0('biomassdraw',1:250)

calib_full <- cbind(cast.x[,1:4],Y,all_biom_draws)

calib_meta <- read.csv('calib.meta.csv') 

calib_tog <- dplyr::full_join(calib_meta,calib_full)

#check biomass distribution
breaks <-  c(seq(0,40,10),seq(50,300,50))
colors <- rev(terrain.colors(length(breaks)-1))
data_binned <-  cut(calib_tog$biomassdraw100, c(breaks), include.lowest = FALSE, labels = FALSE)
plot(calib_tog$long,calib_tog$lat,col=colors[data_binned],pch=19)
map('state',add=T)

#check hemlock distribution
breaks <-  seq(0,100,1)
colors <- rev(terrain.colors(length(breaks)-1))
data_binned <-  cut(calib_tog$TSUGAX, c(breaks), include.lowest = FALSE, labels = FALSE)
plot(calib_tog$long,calib_tog$lat,col=colors[data_binned],pch=19)
map('state',add=T)


write.csv(calib_tog,file='ReFAB_calibration_data_v1.0.csv')


#####
##### BETAS #####
#####

path_to_betas <- '~/Downloads/betas_FULL/'

samples.mixed.keep <- list()

for(i in 1:50){
  load(paste0(path_to_betas,list.files(path_to_betas)[i]))
  samples.mixed.keep[[i]] <- samples.mixed[sample(x = 10000:50000,size = 5),]
}

samps.all <- do.call(rbind,samples.mixed.keep)

write.csv(samps.all,'ReFAB_betas_v1.0.csv')


####
#### Reconstruction Dataset ####
####

load('prediction.data_v6.Rdata')

pred_all <- cbind(x.meta,ten.count,x.bacon)

write.csv(pred_all,'ReFAB_prediction_data_v1.0.csv')

#####
##### Biomass site level reconstructions - Follows from [4.5]average_biomass_figure.R #####
#####

# create and write the netCDF file -- ncdf4 version
library(ncdf4)

# define dimensions
londim <- ncdim_def("lon","degrees_east",as.double(unlist(long)))
latdim <- ncdim_def("lat","degrees_north",as.double(unlist(lat)))
timedim <- ncdim_def("time",units = 'years before present (present = 1950)',
                     as.double((seq(100,10000,100))))
iterdim <- ncdim_def("MCMC iterations",units = '',
                     as.double((seq(1,250,1))))

# define variables
fillvalue <- 1e32
dlname <- "aboveground woody biomass"
tmp_def <- ncvar_def("AGWB","Mg/ha",list(londim,latdim,iterdim,timedim),fillvalue,dlname,prec="double")
tmp_def_name <- ncvar_def("sitename","",list(londim),fillvalue,'site name',prec="double")

# create netCDF file and put arrays
ncfname <- "ReFAB_site_reconstruction_v1.0.nc"
ncout <- nc_create(ncfname,list(tmp_def,tmp_def_name))

all_data <- array(NA,dim=c(80,80,250,100))
for(ii in 1:80){
  if(is.null(all.samps.list[[ii]])) {print(ii);next()}
  all_data[ii,ii,,] <- all.samps.list[[ii]]
}

all_data_1 <- all_data[-c(22,27,73),-c(22,27,73),,]

# put variables
ncvar_put(ncout,tmp_def,all_data_1)
ncvar_put(ncout,tmp_def_name,name.keep[-c(22,27,73)])

# put additional attributes into dimension and data variables
ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout,"lat","axis","Y")
ncatt_put(ncout,"time","axis","T")
ncatt_put(ncout,"MCMC iterations","axis","I")

# add global attributes
ncatt_put(ncout,0,"title",'ReFAB Site Reconstruction Data Product')
ncatt_put(ncout,0,"institution",'University of Notre Dame')
history <- paste(c("A.M. Raiho",'C.J. Paciorek','J.S. McLachlan'), date(), sep=", ")
ncatt_put(ncout,0,"history",history)

# close the file, writing data to disk
nc_close(ncout)

#####
##### Checking site level output #####
#####

library(lattice)
library(RColorBrewer)

nc <- nc_open('ReFAB_site_reconstruction_v1.0.nc')

agwb <- ncvar_get(nc,"AGWB");

nc$dim$lon$vals -> lon
nc$dim$lat$vals -> lat

cutpts <- seq(0,250,length.out = 20)
YBP <- ncvar_get(nc,"time")

####
#### Mean Biomass Map ####
####

par(mfrow=c(1,1))
for(tt in rev(seq(1,100,2))){
  agwb_time_slice_mean <- diag(apply(agwb[,,,tt],2,rowMeans))
  
  data_binned <- cut(agwb_time_slice_mean,breaks = cutpts)
  
  plot(lon, lat, bg = rev(terrain.colors(length(cutpts)))[data_binned], pch = 21,cex = 2)
  title(YBP[tt])
  maps::map('state',add=T)
}

####
#### Variance Biomass Map ####
####
cutpts_var <- seq(0,10000,length.out = 20)

agwb_var_keep <- matrix(NA,77,100)

for(ii in 1:77){
  agwb_var_keep[ii,] <- apply(agwb[ii,ii,,],2,var)
}

par(mfrow=c(1,1))
for(tt in rev(seq(1,100,2))){
  
  data_binned_var <- cut(agwb_var_keep[,tt],breaks = cutpts_var)
  
  #could make circles in previous map colored by variance
  plot(lon, lat, bg = rev(gray.colors(length(cutpts)))[data_binned_var], pch = 21,cex = 2)
  title(YBP[tt])
  maps::map('state',add=T)
}

####
#### Time Series Check ####
####

quants <- array(NA,dim = c(3,100,77))

par(mfrow = c(1, 2))
for (ii in 1:77) {
  quants[, , ii] <- apply(agwb[ii, ii, , ], 2, quantile, c(0.025, .5, .975))
  
  matplot(
    t(quants[, , ii]),
    typ = 'l',
    ylim = c(0, 350),
    xlim = c(100, 0),
    xaxt = 'n',
    col = c(1, 2, 1)
  )
  axis(1,
       at = seq(0, 100, 10),
       labels = seq(0, 10000, 1000))
  title(ii)
  
  plot(
    lon,
    lat,
    bg = rev(terrain.colors(length(cutpts)))[data_binned],
    pch = 21,
    cex = 2
  )
  points(lon[ii],
         lat[ii],
         pch = 6,
         col = 'blue',
         cex = 3)
  title(YBP[tt])
  maps::map('state', add = T)
  
}

####
#### Spatial Reconstruction ####
####

load('full_mat.Rdata')
load('pred_biomass_gam_list.Rdata')
full_mat_se <- cbind(full_mat[,1:2],
                  do.call(cbind,lapply(pred_biomass_gam_list, FUN=function(x) x$se)))

colnames(full_mat_se)[-c(1,2)] <- paste0('AgeYBP',seq(100,10000,100),'_SE')
write.csv(full_mat,file = 'ReFAB_spatial_reconstruction_means_v1.0.csv')
write.csv(full_mat_se,file = 'ReFAB_spatial_reconstruction_SEs_v1.0.csv')



      