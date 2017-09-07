
arg <- commandArgs(trailingOnly = TRUE)
if (is.na(arg[1])) {
  stop()
} else {
  runnum <- as.numeric(arg[1])
}

dataID <- read.csv('dataID.csv')

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
source('~/ReFAB/fit_fix_sigma.R')        # contains fit() function

#adding settlement biomass to the metadata matrix in the most annoying way possible
##TO DO: put somewhere else
library(fields)

locn <- as.character(dataID[dataID$ID==runnum,'name'])
site_number = unique(x.meta[x.meta$site.name == locn,1])
ten_count_use = ten.count[which(x.meta$site.id == site_number), ]

Y = as.matrix(ten_count_use)
minAge = 0
maxAge = 10000
ageInterval = 100

sample_ages <- x.meta[x.meta[,1] == site_number, ]$age_bacon
age_bins <- seq(minAge, maxAge, ageInterval)
age_index <- as.matrix(as.numeric(
  cut(sample_ages, breaks = age_bins, labels=seq(1:(length(age_bins)-1)))
))

tmp <- data.frame(cbind(age_index, Y))
names(tmp)[1] <- 'age_index'

Y2 <- aggregate(tmp, by = list(tmp$age_index), FUN = sum)

set.seed(0)
group.sample <- sample(x = 1:nrow(Y2), size = nrow(Y2), replace = FALSE)
group.mat <- matrix(group.sample[1:(round((nrow(Y2) / 10))*10)],
                    ncol = round((nrow(Y2) / 10)))
group.mat[is.na(group.mat)] <- sample(x = 1:nrow(Y2), size = length(which(is.na(group.mat))))

sigma <- as.numeric(dataID[dataID$ID==runnum,'sigma'])
group <- as.numeric(dataID[dataID$ID==runnum,'group'])
if(is.na(group)){
  group <- NULL
} 
  
smp <- fit_fix_sigma(locn = locn, pred_code_fix_sigma = pred_code_fix_sigma,
                     pred_code_fix_b = pred_code_fix_b, order = 3, Z = Z,
                     u = u, x.meta = x.meta,
                     ten.count = ten.count, beta1 =  beta1.est.real,
                     beta2 = beta2.est.real,
                     nIts = 50000, nItsSave = 10000, seed = 1,
                     control.pts = control.pts, sigma = sigma,
                     group = group, group.mat = group.mat)
