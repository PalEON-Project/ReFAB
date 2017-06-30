library(nimble)
library(RCurl)
library(maps)

ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}
control.pts<-read.csv(text=getURL('https://raw.githubusercontent.com/PalEON-Project/stepps-baconizing/master/data/bacon_age_markers_v1.csv'))


# nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)

# load in data for all sites, per Ann's original code
if(!file.exists('allPredData.Rda'))
  source('prep_data.R') 
load('allPredData.Rda')

site.names <- unique(x.meta$site.name)

source('~/ReFAB/genPareto/model_dgp_auxil.R')  # BUGS code for model
source('~/ReFAB/genPareto/fit_model.R')        # contains fit() function
for(i in 75:length(site.names)){ #44:length(site.names)
  locn <- site.names[i]
  site_number = unique(x.meta[x.meta$site.name == locn,1])
  ten_count_use = ten.count[which(x.meta$site.id == site_number), ]
  
  Y = as.matrix(ten_count_use)
  if(length(Y)>21 & nrow(Y) > 15 &
     max(x.meta[x.meta$site.name == locn,'age_bacon'])>8000 & 
     min(x.meta[x.meta$site.name == locn,'age_bacon'])<2000){
    
    print(paste('Getting there!',signif(i/length(site.names)*100,digits=2),'% complete'))
    
    smp <- fit(locn = locn, pred_code = pred_code, order = 3, Z = Z,
               u = u, x.meta = x.meta,
               ten.count = ten.count, beta1 =  beta1.est.real,
               beta2 = beta2.est.real,
               nIts = 50000, nItsSave = 10000, seed = 1,
               control.pts = control.pts)
  }
}