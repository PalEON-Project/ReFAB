library(nimble)

# nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)

# load in data for all sites, per Ann's original code
if(!file.exists('allPredData.Rda'))
  source('prep_data.R') 
load('allPredData.Rda')

locn <- 'Cub Lake'

source('model_dgp_auxil.R')  # BUGS code for model
source('fit_model.R')        # contains fit() function
smp <- fit(locn, pred_code, order = 3, Z, u, x.meta, ten.count, beta1.est.real, beta2.est.real, nIts = 5000, nItsSave = 1000, seed = 1)

