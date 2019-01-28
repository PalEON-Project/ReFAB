
#####
##### This runs all prediction sites with a job array on the cluster
#####

arg <- commandArgs(trailingOnly = TRUE)
if (is.na(arg[1])) {
  runnum <- NA
} else {
  runnum <- as.numeric(arg[1])
}
dataDir <- c(getwd()) #or wherever allPredData.Rda is located

## each runnum will be a different beta estimate with a different site so 
## we'll have length(runnum) == nrow(beta.est)*length(sites)

#dataID <- read.csv(file.path('Cross_Validation','dataID.csv')) #for 10K
#dataID <- read.csv(file.path('Cross_Validation','beta_run_dataID.csv')) #for paleon mip
#dataID <- read.csv('dataID_bacon_v4_means.csv') #for paleon mip

#### These files come from [2]create_prediction_datasets.R

dataID <- read.csv('dataID_v5.csv') #dataID <- read.csv('dataID_bacon_new_v1.csv') #for original preds dataID <- 
load('prediction.data_v5.Rdata') #load('prediction.data_new_v1.Rdata') #

source(file.path('Workflow_Code','utils','validation_args.R')) #file with constants that should be constant between validation exercises

####
#### Master Setup
####

library(nimble)
library(RCurl)
library(maps)
library(fields)
library(sp)

control.pts<-read.csv(file.path('Data','control.pts.csv'))

# load in data for all sites, per Ann's original code
# if(!file.exists(file.path(dataDir,'allPredData.Rda')))
# source(file.path('Workflow_Code','prep_data.R'))
# load(file.path(dataDir,'allPredData.Rda')) # load for 10K run
# load('paleon.data.Rdata') # adding via ssh rather than github
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


beta <- dataID[dataID$ID==runnum,'beta']
load(file = paste0("beta.est.group.in",seq(100,150,1)[beta],"FULL.Rdata")) #load via ssh

i.beta1 <- grep("beta1",colnames(samples.mixed))
i.beta2 <- grep("beta2",colnames(samples.mixed))
last_get <- nrow(samples.mixed)

beta1.est.real = matrix(samples.mixed[last_get,i.beta1],ncol(Z),ncol(ten.count))
beta2.est.real = matrix(samples.mixed[last_get,i.beta2],ncol(Z),ncol(ten.count))


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

minAge = 0
maxAge = 10000
ageInterval = 100

site_number = unique(x.meta[x.meta$site.name == locn,1])
x.meta.use <- x.meta[x.meta$site.name == locn,]

x.bacon.use <- x.bacon[x.meta$site.name == locn,]
if(!is.na(beta)){
  which_bacon <- dataID[dataID$ID==runnum,'beta']
  sample_ages <- x.bacon.use[,which_bacon]#x.meta.use$age_bacon
}else{
  sample_ages <- rowMeans(x.bacon.use)
}

source(file.path('Workflow_Code','utils','test_site.R'))
test_site(x.meta.use)

ten_count_use = ten.count[which(x.meta$site.name == locn), ]
ten_count_use[which(is.na(ten_count_use))] <- 0
Y = as.matrix(ten_count_use)

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
  group.mat <- NULL
}

smp <- fit_fix_sigma(locn = locn, pred_code_fix_sigma = pred_code,
                     pred_code_fix_b = pred_code_fix_b, order = 3, Z = Z,
                     u = u, x.meta = x.meta,
                     ten_count_use = ten_count_use,
                     beta1 =  beta1.est.real,
                     beta2 = beta2.est.real,
                     nIts = Niters, nItsSave = .2*Niters, seed = 1,
		                 control.pts = control.pts, sigma = .03, lik.only = FALSE,
                     maxAge = 10000, Nbeta = beta, ID = runnum,
		                 liks.by.taxa = TRUE, bMax = bMax, number.save = 250,
		                 get.log.prob = FALSE, group = NULL)

if(FALSE){
  ## Simple Site Diagnositics
  plot(colMeans(smp[,1:100]),col='red',pch=19,ylim=c(0,150))
  load(workInfo) #to get out load the workInfo file make in fit_fix_sigma
  points(age_index,seq(5, bMax-5, by = 2)[apply(out,2,which.max)])
}


