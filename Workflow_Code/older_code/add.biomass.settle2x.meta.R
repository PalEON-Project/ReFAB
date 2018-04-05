
biomass_dat_est <- read.csv(file.path('Data',"biomass_prediction_v0.9-7_bam.csv"))
biomass_dat_est_uncert <- read.csv(file.path('Data',"biomass_uncertainty_v0.9-7_bam.csv"))

##### Changing pollen coordinates so that we can find the right rows when we find the biomass for each pond.
lat.long.reg <- cbind(as.numeric(as.character(x.meta$long)),as.numeric(as.character(x.meta$lat)))
lat.long.reg.df = data.frame(lat.long.reg)
colnames(lat.long.reg.df) = c('x', 'y')

coordinates(lat.long.reg.df) <- ~ x + y
proj4string(lat.long.reg.df) <- CRS('+proj=longlat +ellps=WGS84')

albers <- spTransform(lat.long.reg.df, CRS('+init=epsg:3175'))
albers <- as.matrix(data.frame(albers))

centers_biomass = cbind(biomass_dat_est$x,biomass_dat_est$y)
idx_cores = vector(length=nrow(x.meta))

for(i in 1:nrow(x.meta)){   
  core_site = albers[i,]
  d = rdist(matrix(core_site, ncol=2), as.matrix(centers_biomass))
  if(min(d)<8000){
    idx_cores[i] = which.min(d) 
  }else{
    idx_cores[i] = NA 
  }
  
}

biomass <- biomass.sd <- list()
for(i in 1:nrow(x.meta)){ 
  biomass[[i]] = c(biomass_dat_est[idx_cores[i],'Total'])
  biomass.sd[[i]] =c(biomass_dat_est_uncert[idx_cores[i],'Total'])
}

x.meta <- cbind(x.meta,unlist(biomass,use.names=FALSE),unlist(biomass.sd,use.names=FALSE))
colnames(x.meta)[(ncol(x.meta)-1):ncol(x.meta)]<-c('SettleBiomassMean','SettleBiomassSD')
