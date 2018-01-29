
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

#beta1.est.real = matrix(colMeans(samples.mixed[,i.beta1]),5,ncol(ten.count))
#beta2.est.real = matrix(colMeans(samples.mixed[,i.beta.pine]),5,ncol(ten.count))

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

source('test_site.R')
test_site(x.meta.use)

ten_count_use = ten.count[which(x.meta$site.name == locn), ]
ten_count_use[which(is.na(ten_count_use))] <- 0
Y = as.matrix(ten_count_use)

sample_ages <- x.meta.use$age_bacon
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
                     nIts = 10000, nItsSave = 2000, seed = 1,
		     control.pts = control.pts, sigma = sigma,
                     group = group, group.mat = group.mat, lik.only = FALSE,
                     maxAge = 10000, Nbeta = beta, ID = runnum)

stop()

source('site_diag.R')

## Full Site Diagnostics
for(l in unique(dataID$name)){
  site_diag(locn = l, x.meta = x.meta, ten.count = ten.count,
            control.pts = control.pts, bMax = 150)
}

## Individual Site Diagnostics
site_diag(locn = 'Cub Lake', x.meta = x.meta, ten.count = ten.count,
          control.pts = control.pts, bMax = 150,
          path_to_samps = "~/ReFAB/samplesList_workInfo_1_Cub-Lake_Beta_20.Rda",
          path_to_Info = "~/ReFAB/workInfo_1_Cub-Lake_Beta_20.Rda")

## Simple Site Diagnositics
plot(colMeans(samplesList[,1:100]),col='red',pch=19,ylim=c(0,150))
points(age_index,seq(5, 150-5, by = 2)[apply(out,2,which.max)])

## Drawing Multiple MCMCs
load("~/Downloads/samps/samplesList_workInfo_181_Cub-Lake_Beta_1.Rda")
load("~/Downloads/work/workInfo_181_Cub-Lake_Beta_1.Rda")
out.save <- out
samplesList.save <- samplesList
load("~/Downloads/samps/samplesList_workInfo_184_Cub-Lake_Beta_4.Rda")
load("~/Downloads/work/workInfo_184_Cub-Lake_Beta_4.Rda")
samplesList.save.1 <- samplesList
out.save.1 <- out
load("~/Downloads/samps/samplesList_workInfo_194_Cub-Lake_Beta_14.Rda")
load("~/Downloads/work/workInfo_194_Cub-Lake_Beta_14.Rda")

pdf('Cub-Lake-MCMC.pdf')
par(mfrow=c(2,2))
for(i in 1:100){
  plot(samplesList[,i],typ='l',lwd=1.5,ylim=c(0,bMax))
  points(samplesList.save[,i],typ='l',lwd=1.5,col='red') 
  points(samplesList.save.1[,i],typ='l',lwd=1.5,col='blue') 
  title(i)
  if(any(i==age_index)){
    abline(h = seq(5, 150-5, by = 2)[apply(out,2,which.max)][which(i==age_index)],lwd=2)
    abline(h = seq(5, 150-5, by = 2)[apply(out.save,2,which.max)][which(i==age_index)],col='red',lwd=2)
    abline(h = seq(5, 150-5, by = 2)[apply(out.save.1,2,which.max)][which(i==age_index)],col='blue',lwd=2)
  }
}
dev.off()

