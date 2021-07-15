


dataID <- read.csv('dataID_v5.csv') #dataID <- read.csv('dataID_bacon_new_v1.csv') #for original preds dataID <- 
load('prediction.data_v6.Rdata') #load('prediction.data_new_v1.Rdata') #

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

b = 2

load(file = paste0('~/Dropbox/',length(u),"beta.est.group.in",seq(100,150,1)[b],"FULL.Rdata")) #load via ssh

i.beta1 <- grep("beta1",colnames(samples.mixed))
i.beta2 <- grep("beta2",colnames(samples.mixed))
last_get <- nrow(samples.mixed)

beta1.est.real = matrix(samples.mixed[last_get,i.beta1],ncol(Z),ncol(ten.count))
beta2.est.real = matrix(samples.mixed[last_get,i.beta2],ncol(Z),ncol(ten.count))

####
#### Start One Site Model Run ####
####


#load(file.path(dataDir,'x.meta.w.settle.Rdata')) #load for 10K

beta_mat_keep <- list()

minAge = 0
maxAge = 10000
ageInterval = 100

source(file.path('Workflow_Code','models','model_dgp_auxil.R')) # BUGS code for model
source(file.path('Workflow_Code','models','run_prediction.R')) # contains fit_fix_sigma() function

for(ss in 1:80){
  
  locn <- as.character(unique(dataID$name)[ss])
  if(locn == 'Lily Lake' | locn == 'Mud Lake' | locn == 'Seidel') next()
 
  site_number = unique(x.meta[x.meta$site.name == locn, 1])
  x.meta.use <- x.meta[x.meta$site.name == locn, ]
  
  x.bacon.use <- x.bacon[x.meta$site.name == locn, ]
  if (!is.na(beta)) {
    which_bacon <- dataID[dataID$ID == b, 'beta']
    sample_ages <- x.bacon.use[, which_bacon]#x.meta.use$age_bacon
  } else{
    sample_ages <- rowMeans(x.bacon.use)
  }
  
  source(file.path('Workflow_Code', 'utils', 'test_site.R'))
  test_site(x.meta.use)
  
  ten_count_use = ten.count[which(x.meta$site.name == locn),]
  ten_count_use[which(is.na(ten_count_use))] <- 0
  Y = as.matrix(ten_count_use)
  
  age_bins <- seq(minAge, maxAge, ageInterval)
  age_index <- as.matrix(as.numeric(cut(
    sample_ages, breaks = age_bins, labels = seq(1:(length(age_bins) - 1))
  )))
  
  
  ID <- dataID[dataID$name == as.character(locn), 'ID'][b]
  path_to_samps <-
    c('~/Downloads/samps_again_FULL/')#c('~/ReFAB/samps_final/')
  locnClean <- gsub(' ', '-', locn)
  load(paste0(path_to_samps,
         'samplesList_workInfo_',
         ID,
         '_',
         locnClean,
         '_Beta_',
         b,
         '.Rdata'))
  
  
  tmp_y <- data.frame(cbind(age_index, Y))
  names(tmp_y)[1] <- 'age_index'
  
  Y2 <- aggregate(tmp_y, by = list(tmp_y$age_index), FUN = sum)
  
  Y <- as.matrix(Y2[ , -c(1,2)])
  age_index <- Y2[,1]
  
  biomass <- colMeans(samplesList[, age_index])
  
  betabin_mat <- matrix(NA, nrow(Y), ncol(Y))
  for (rr in 1:nrow(Y)) {
    betabin_mat[rr, ] <-
      get_dbetabin(
        biomass[rr],
        u = u,
        N0 = rep(0, (length(u) - 1)),
        N1 = rep(0, (length(u))),
        N2 = rep(0, (length(u) + 1)),
        N3 = rep(0, (length(u) + 2)),
        beta1 = beta1.est.real,
        beta2 = beta2.est.real,
        n = rowSums(Y)[rr],
        Y = as.matrix(Y[rr, ])
      )
  }
  colnames(betabin_mat) <- paste0('dens_', colnames(Y))
  
  beta_mat_keep[[ss]] <-
    data.frame(betabin_mat, age = age_index, site = rep(locn, nrow(Y)), Y,lat = rep(x.meta.use$lat[1],nrow(Y)),lon=rep(x.meta.use$long[1],nrow(Y)))
}


all_betabin <- do.call(rbind,beta_mat_keep)

par(mfrow=c(1,3))
for(cc in 1:3){
  clust_get <- names(which(clusters@cluster==cc))
  
  bar_mat <- rbind(colMeans(all_betabin[which(all_betabin$age>80 & all_betabin$site %in% clust_get),1:22]),
        colMeans(all_betabin[which(all_betabin$age<=80 & all_betabin$age>=50& all_betabin$site %in% clust_get),1:22]),
        colMeans(all_betabin[which(all_betabin$age<50 & all_betabin$site %in% clust_get),1:22]))
  
  barplot(abs(bar_mat),beside=T,las=2)
  title(paste('Cluster',cc))
}

col_get <- numeric(nrow(all_betabin))
for(ii in 1:nrow(all_betabin)) col_get[ii] <- adjustcolor('black',abs(all_betabin[,"dens_TSUGAX"])[ii])

keep_min <- numeric(78)
for(ss in 1:78){
  keep_min[ss] <- order(colSums(all_betabin[all_betabin$site==unique(all_betabin$site)[ss] & all_betabin$age<=50,1:22]))[3]
}
table(colnames(Y)[keep_min])

plot(all_betabin[,'lon'],
     all_betabin[,'lat'],
     col = keep_min,pch=19)
maps::map('state',add=T)
