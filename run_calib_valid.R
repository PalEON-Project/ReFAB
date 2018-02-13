
arg <- commandArgs(trailingOnly = TRUE)
if (is.na(arg[1])) {
  runnum <- NA
} else {
  group_rm <- as.numeric(arg[1])
}

### after get.data

library(nimble)
library(splines)
library(maps)
library(methods)

ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}

load("2018-01-08calibration.data.Rdata") #this needs to have cast.x in it. and the entire biomass and sites_rm
load("cast.x.Rdata")
load("sites_rm.Rdata")

all.pollen.taxa.names <- colnames(cast.x)[5:84]
trees <- c("FAGUS","TSUGAX","QUERCUS","BETULA",
           'PINUSX',"JUGLANSX","ACERX","FRAXINUX",
           "OSTRYCAR","ULMUS","TILIA","ALNUSX",
           "CYPERACE","PICEAX",
           "ABIES","POPULUS","CARYA",
           "LARIXPSEU","TAXUS","NYSSA","CASTANEA","PLATANUS","SALIX",
           "LIQUIDAM","CUPRESSA")
other.trees <- c()#NULL#c()
drop.taxa <- NA#c('other_herbs')

source('taxa_selection.R')
Y <- taxa_selection(trees = trees, other.trees = other.trees,
                    cast.x = cast.x, sites_rm = sites_rm,
                    all.pollen.taxa.names = all.pollen.taxa.names,
                    prairie.include = T, other.herbs.include = T,
                    other.trees.include = T, drop.taxa = drop.taxa,
                    PFT.do = F)

Niters <- 5000
bMax <- 143

#### Setting up 10 fold cross validation
set.seed(5)
sets10 <- replicate(10, sample(size = 10,x = nrow(Y), replace = F))
Y.keep <- Y
biomass.keep <- biomass
Y.calib <- Y[-sets10[,group_rm],]; Y.pred <- Y[sets10[,group_rm],]
biomass.calib <- biomass[-sets10[,group_rm]]; biomass.pred <- biomass[sets10[,group_rm]]

#### Making sure Z.knots and u are the same between calibration and validation
Z.knots = bs(biomass.calib,intercept=TRUE,df=5)
u <- c(rep(attr(Z.knots,"Boundary.knots")[1],1),attr(Z.knots,"knots"),rep(attr(Z.knots,"Boundary.knots")[2],1))

source(file.path('Workflow_Code','calibration.model.R'))
samples.mixed <- calibration_model(Y = Y.calib, biomass = biomass.calib,
                                   Z.knots = Z.knots, u = u, Niters = Niters,
                                   group_rm = group_rm)

source('validation.R')
samples.pred <- validation_model(Y = Y.pred, Z.knots = Z.knots, 
                 samples.mixed = samples.mixed, u = u,
                 Niters = Niters, bMax = bMax)

load(file=paste0('outLik.group.',group_rm,'.Rdata'))

source('calibration.figs.R')
calibration.figs(bMax = bMax, Z.knots = Z.knots, Y = Y.keep,
                 samples.mixed = samples.mixed, outLik = outLik,
                 biomass = biomass.keep, samples.pred = samples.pred,
                 group_rm = group_rm,Y.pred = Y.pred,
                 biomass.pred = biomass.pred)



