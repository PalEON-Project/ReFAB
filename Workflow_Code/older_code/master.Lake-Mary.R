
library(nimble)
library(RCurl)
library(maps)

ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}
control.pts<-read.csv('~/ReFAB/Data/control.pts.csv')

# nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)

# load in data for all sites, per Ann's original code
if(!file.exists('allPredData.Rda'))
  source('prep_data.R') 
load('allPredData.Rda')

source('~/ReFAB/genPareto/model_dgp_auxil.R')  # BUGS code for model
source('~/ReFAB/genPareto/fit_model.R')        # contains fit() function

#adding settlement biomass to the metadata matrix in the most annoying way possible
##TO DO: put somewhere else
library(fields)
biomass_dat_est <- read.csv(paste0("~/babySTEPPS/Data/biomass_prediction_v0.9-7_bam.csv"))

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

biomass <- list()
for(i in 1:nrow(x.meta)){ 
  biomass[[i]] = biomass_dat_est[idx_cores[i],'Total']
}

x.meta <- cbind(x.meta,unlist(biomass))
colnames(x.meta)[ncol(x.meta)]<-c('SettleBiomass')

#stewerts dark\
#lake O' Pines

for(i in 57:length(how.many)){
  Lake Mary <- 'Martha Lake'#names(how.many)[i]

  locn <- Lake Mary
  site_number = unique(x.meta[x.meta$site.name == locn,1])
  ten_count_use = ten.count[which(x.meta$site.id == site_number), ]
  
  Y = as.matrix(ten_count_use)
 
  smp <- fit(minAge = 0, maxAge = 2000, locn = locn,
             pred_code = pred_code, order = 3, Z = Z,
               u = u, x.meta = x.meta,
               ten.count = ten.count, beta1 =  beta1.est.real,
               beta2 = beta2.est.real,
               nIts = 50000, nItsSave = 10000, seed = 1,
               control.pts = control.pts)
}
