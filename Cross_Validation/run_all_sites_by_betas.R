
dataDir <- c(getwd()) #or wherever allPredData.Rda is located

arg <- commandArgs(trailingOnly = TRUE)
if (is.na(arg[1])) {
  runnum <- NA
} else {
  runnum <- as.numeric(arg[1])
}

## each runnum will be a different beta estimate with a different site so 
## we'll have length(runnum) == nrow(beta.est)*length(sites)

#dataID <- read.csv(file.path('Cross_Validation','dataID.csv')) #for 10K
dataID <- read.csv(file.path('Cross_Validation','beta_run_dataID.csv')) #for paleon mip

load('prediction.data.Rdata')
####
#### Master Setup
####

library(nimble)
library(RCurl)
library(maps)
library(fields)
library(sp)

ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}

control.pts<-read.csv(file.path('Data','control.pts.csv'))

# load in data for all sites, per Ann's original code
#if(!file.exists(file.path(dataDir,'allPredData.Rda')))
#  source(file.path('Workflow_Code','prep_data.R'))
# load(file.path(dataDir,'allPredData.Rda')) # load for 10K run
#load('paleon.data.Rdata') # adding via ssh rather than github

load(file = paste0("nimble.betas_1_2_horiz_plus2017-12-21.Rdata")) #load via ssh

i.beta <- grep("beta",colnames(samples.mixed))
i.beta.pine <- grep("beta.pine",colnames(samples.mixed))
i.beta1 <- i.beta[-i.beta.pine]

beta <- dataID[dataID$ID==runnum,'beta']

Nbeta <- round(seq(3000,nrow(samples.mixed),length.out = 20))[beta]

beta1.est.real = matrix(samples.mixed[Nbeta,i.beta1],5,ncol(ten.count))
beta2.est.real = matrix(samples.mixed[Nbeta,i.beta.pine],5,ncol(ten.count))

source(file.path('genPareto','model_dgp_auxil.R')) # BUGS code for model
source(file.path('Cross_Validation','fit_fix_sigma.R')) # contains fit_fix_sigma() function

####
#### Start One Site Model Run ####
####

if(!is.na(runnum)){
  locn <- as.character(dataID[dataID$ID==runnum,'name'])
}else{
  locn <- readline(prompt = 'Location Name:')
}
#load(file.path(dataDir,'x.meta.w.settle.Rdata')) #load for 10K

site_number = unique(x.meta[x.meta$site.name == locn,1])
ten_count_use = ten.count[which(x.meta$site.id == site_number), ]
ten_count_use[which(is.na(ten_count_use))] <- 0

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

if(!is.na(runnum)){
  sigma <- as.numeric(dataID[dataID$ID==runnum,'sigma'])
  group <- as.numeric(dataID[dataID$ID==runnum,'group'])
}else{
  sigma <- as.numeric(readline(prompt = 'sigma = '))
  group <- NULL
}

smp <- fit_fix_sigma(locn = locn, pred_code_fix_sigma = pred_code_fix_sigma,
                     pred_code_fix_b = pred_code_fix_b, order = 3, Z = Z,
                     u = u, x.meta = x.meta,
                     ten_count_use = ten_count_use,
                     beta1 =  beta1.est.real,
                     beta2 = beta2.est.real,
                     nIts = 20000, nItsSave = 200, seed = 1,
                     control.pts = control.pts, sigma = sigma,
                     group = group, group.mat = group.mat, lik.only = FALSE,
                     maxAge = 10000, Nbeta = beta, ID = runnum)


