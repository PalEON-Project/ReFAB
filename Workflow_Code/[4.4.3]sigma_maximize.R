
dataDir <- c(getwd()) #or wherever allPredData.Rda is located

####
#### This script runs all sites all groups all sigmas to maximize over the log likelihood of the sigma values
#### Only for running job.array.sigma.sh ####
####

arg <- commandArgs(trailingOnly = TRUE)
if (is.na(arg[1])) {
  runnum <- NA
} else {
  runnum <- as.numeric(arg[1])
}

dataID <- read.csv(file.path('Cross_Validation','dataID.csv')) #for 10K #dataID_fine

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

load('prediction.data_v5.Rdata')
group_rm <- dataID[dataID$ID==runnum, 'group']

source(file.path('Workflow_Code','utils','validation_args.R')) #file with constants that should be constant between validation exercises
load(file = paste0(length(u),"beta.est.group.in", group_rm, ".Rdata"))

i.beta1 <- grep("beta1",colnames(samples.mixed))
i.beta2 <- grep("beta2",colnames(samples.mixed))
burnin <- .2*nrow(samples.mixed)

source(file.path('Workflow_Code','models','model_dgp_auxil.R')) # BUGS code for model
source(file.path('Workflow_Code','models','run_prediction.R')) # contains run_prediction() function

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

# for cross validation #removes samples
  set.seed(0)
  group.sample <- sample(x = 1:nrow(Y2), size = nrow(Y2), replace = FALSE)
  group.mat <- matrix(group.sample[1:(round((nrow(Y2) / 10))*10)],
                      ncol = round((nrow(Y2) / 10)))
  group.mat[is.na(group.mat)] <- sample(x = 1:nrow(Y2), size = length(which(is.na(group.mat))))


if(!is.na(runnum)){
  sigma <- as.numeric(dataID[dataID$ID==runnum,'sigma'])
  group <- as.numeric(dataID[dataID$ID==runnum,'group'])
}else{
  sigma <- as.numeric(readline(prompt = 'sigma = '))
  group <- NULL
}
  
  new.biomass <- 1:bMax
  Z = matrix(0,nrow=length(new.biomass),ncol=5)
  source(file.path('Workflow_Code','utils','bs_nimble.R'))
  for(i in 1:length(new.biomass)){
    u_given <- new.biomass[i]
    Z[i,] = bs_nimble(u_given, u=u, N0 = rep(0, (length(u)-1)),
                      N1 = rep(0, (length(u))), 
                      N2 = rep(0, (length(u)+1)), 
                      N3 = rep(0, (length(u)+2)))
  }  
  
  beta <- NA
  
  if(!is.na(beta)){
    Nbeta <- round(seq(8000,nrow(samples.mixed),length.out = 20))[beta]
    
    beta1.est.real = matrix(samples.mixed[Nbeta,i.beta1],ncol(Z),ncol(ten.count))
    beta2.est.real = matrix(samples.mixed[Nbeta,i.beta2],ncol(Z),ncol(ten.count))
  }else{
    beta1.est.real = matrix(colMeans(samples.mixed[burnin:nrow(samples.mixed),i.beta1]),ncol(Z),ncol(Y))
    beta2.est.real = matrix(colMeans(samples.mixed[burnin:nrow(samples.mixed),i.beta2]),ncol(Z),ncol(Y))
  }
  
  
smp <- run_prediction(locn = locn, pred_code_fix_sigma = pred_code_fix_sigma,
                     pred_code_fix_b = pred_code_fix_b, order = 3, Z = Z,
                     u = u, x.meta = x.meta,
                     ten_count_use = ten_count_use, beta1 =  beta1.est.real,
                     beta2 = beta2.est.real,
                     nIts = 5000, nItsSave = 2500, seed = 1,
                     control.pts = control.pts, sigma = sigma,
                     group = group, group.mat = group.mat, lik.only = FALSE,
                     bMax = bMax,
                     maxAge = 10000, ID = runnum, number.save = 250, get.log.prob = TRUE)

