
library(dplyr)

### loading nimble pull
load("~/ReFAB/nimble_pull2018-10-31.Rdata")
pol_cal_count <- pol_cal_count %>% mutate(.id = as.numeric(.id))
pol_cal_count_save <- pol_cal_count

### loading bacon data created by Andria
bacon <- read.csv(file.path('Data','sediment_ages_v1.0_varves.csv')) #andria bacon

### merging bacon dating with neotoma raw data
new.pollen <- bacon %>% 
  #dplyr::select(id, age_default, depth, starts_with('bacon')) %>% 
  left_join(pol_cal_count,
            by = c('id' = '.id', 'PINE' = 'PINUSX', 'depth'))

### selecting only the data we want to use for our 10000 year prediction
age_bacon <- apply(new.pollen[,grep('bacon',colnames(new.pollen))],1,range)

dataset_ids <- unique(new.pollen$id)
pol_keep <- pol_not_used <- list()
for(i in 1:length(dataset_ids)){
  if(any(which(new.pollen$id==dataset_ids[i]))){
    mat_use <- new.pollen[which(new.pollen$id==dataset_ids[i]),]
    #if(any(which(mat_use$lat.x!=mat_use$lat.y))) stop()
    age_bacon <- apply(mat_use[,grep('bacon',colnames(new.pollen))],1,quantile,c(.25,.75))
    if(min(age_bacon[1,])<=2000&max(age_bacon[2,])>=9000&nrow(mat_use)>10){
      pol_keep[[i]] <- mat_use
    }else{
      pol_not_used[[i]] <- mat_use
    }
  }
}

pol_pred <- do.call(rbind,pol_keep)
pol_no <- do.call(rbind,pol_not_used)

x = pol_pred

quantile(table(x$dataset))
hist(table(x$dataset))
length(which(table(x$dataset) > 10))/length(unique(x$dataset))*100
x <- x[which(x$dataset%in%names(which(table(x$dataset)>10))),]

#####
##### Find which sites in MN, WI, and MI need to be baconized
#####

no_bacon <- pol_cal_count[-which(pol_cal_count$.id%in%x$id),]

ids_no_bacon <- unique(no_bacon$.id)

no_bacon_keep <- list()

for(i in 1:length(ids_no_bacon)){
  
  no_bacon_look <- no_bacon[which(no_bacon$.id==ids_no_bacon[i]),]

  if(min(no_bacon_look$age) <=2000 &
     max(no_bacon_look$age) >= 9000 &
     nrow(no_bacon_look) >= 10 & !is.na(no_bacon_look$age)){
    no_bacon_keep[[i]] <- no_bacon_look
  }
}

no_bacon_make <- do.call(rbind,no_bacon_keep)

###Getting rid of IN and IL sites
pollen_v1 <- readRDS("~/bulk-baconizing/data/pollen_v1.rds")
ids <- names(pollen_v1)

no_bacon_make <- no_bacon_make[-which(no_bacon_make$.id%in%ids),]

maps::map('state', ylim=range(new.pollen$lat.x)+c(-2, 2), xlim=range(new.pollen$long.x)+c(-2, 2),main=NA)
points(no_bacon_make$long, no_bacon_make$lat, pch=19, cex=1,col="gray")


new.ids <- unique(no_bacon_make$.id)

save(new.ids,file='~/bulk-baconizing/new.ids.Rdata')

##loading new MN, MI, and WI sites with bacon posts

pollen_v2 <- readRDS("~/bulk-baconizing/data/pollen_v2.rds")
source('~/bulk-bacon-my-fork/R/agg_posts_counts.R')

new_mnwimi <- make_posts_counts(pollen_v1 = pollen_v2,n.samps = 50)

new_mnwimi1 <- new_mnwimi[-which(is.na(new_mnwimi[,ncol(new_mnwimi)])),]

new_mnwimi2 <- new_mnwimi1[-which(new_mnwimi1$site.name%in%names(which(table(new_mnwimi1$site.name)<10))),]

colnames(new_mnwimi2)[ncol(new_mnwimi2)] <- 'age_bacon'

x.meta = new_mnwimi2[,c('.id','lat',"long","dataset","site.name","age_bacon")]
colnames(x.meta) <- c('site.id','lat','long','dataset','site.name','age_bacon')

x.bacon <- new_mnwimi2[,grep(pattern = 'bacon_draw',colnames(new_mnwimi2))]

### need just the pollen data to give to organizational function
all.pollen.taxa.names <- colnames(new_mnwimi2)[11:78]
pred.x <- new_mnwimi2[,which(colnames(new_mnwimi2)%in%all.pollen.taxa.names)]

### organizing and aggregating pollen data
trees <- c("JUGLANSX","FRAXINUX","OSTRYCAR","ULMUS","TILIA","CARYA",
           "FAGUS","TSUGAX","QUERCUS","BETULA",
           'PINUSX',"ACERX","ALNUSX",
           "CYPERACE","PICEAX","ABIES","POPULUS",
           "LARIXPSEU","CUPRESSA") #
other.trees <- c("CASTANEA","PLATANUS","SALIX","LIQUIDAM","TAXUS","NYSSA")#NULL#c()
drop.taxa <- NA#c('other_herbs')

source(file.path('Workflow_Code','utils','taxa_selection.R'))
ten.count <- taxa_selection(trees = trees, other.trees = other.trees,
                            cast.x = pred.x, sites_rm = 0, bigwoods.include = F,
                            all.pollen.taxa.names = all.pollen.taxa.names,
                            prairie.include = T, other.herbs.include = T,
                            other.trees.include = T, drop.taxa = drop.taxa,
                            PFT.do = F)
### Needs to be in same order as calibration so betas match with the right columns
load("threethirds_v2.0.Rdata")
ten.count <- ten.count[,colnames(Y)] #would it be better to sort originally based off of the prediction datasets?

save(x.meta,x.bacon,ten.count,file='prediction.data_MN_WI_MI_v1.Rdata')

#####
##### Add new baconized sites from IN and IL
#####

pollen_v1 <- readRDS("~/bulk-baconizing/data/pollen_v1.rds")
comp.tax <- neotoma::compile_taxa(pollen_v1, 'WhitmoreSmall')
pol_cal_count2 <- neotoma::compile_downloads(comp.tax)
length(unique(pol_cal_count2$dataset))

bacon_df <- matrix(NA,nrow = nrow(pol_cal_count2),ncol=50)

handles <- sapply(pollen_v1, function(x) { x$dataset$dataset.meta$collection.handle })
ids <- names(pollen_v1)

for(i in 1:length(handles)){
  dir_path <- file.path('~/bulk-baconizing', 'Cores', handles[i])
  files_get <- list.files(dir_path)
  pick_file <- grep(files_get, pattern = 'posteriorout.csv')
  
  if(any(pick_file)){
    posts <- read.csv(file.path(dir_path, files_get[pick_file]))
    bacon_df[which(pol_cal_count2$dataset==ids[i]),] <- as.matrix(posts[,sample(x = 1:1000,size = 50)])
  }else{
    print(paste(handles[i],'not run in bacon'))
  }
}

pol_hold <- cbind(pol_cal_count2,bacon_df,rowMeans(bacon_df))
new.pol <- list()

for(i in 1:length(ids)){
  if(min(bacon_df[which(pol_hold$dataset==ids[i]),],na.rm = T) <=2000 &
  max(bacon_df[which(pol_hold$dataset==ids[i]),],na.rm = T) >= 9000 &
  length(which(!is.na(bacon_df[which(pol_hold$dataset==ids[i]),1]))) >= 10){
    new.pol[[i]] <- pol_hold[which(pol_hold$dataset==ids[i]),]
  }
}

np <- do.call(rbind,new.pol)

colnames(np)[87:137] <- c(paste0('bacon_draw',1:50),'age_bacon')

x.meta = np[,c('.id','lat',"long","dataset","site.name","age_bacon")]
colnames(x.meta) <- c('site.id','lat','long','dataset','site.name','age_bacon')

x.bacon <- np[,grep(pattern = 'bacon_draw',colnames(np))]

## should probably do something different than this
for(i in 1:nrow(x.meta)){
  if(is.na(x.meta$age_bacon[i])){
    x.meta$age_bacon[i] <- np$age[i]
    x.bacon[i,] <-  np$age[i]
  } 
}

### need just the pollen data to give to organizational function
all.pollen.taxa.names <- colnames(pol_cal_count)[11:length(colnames(pol_cal_count))]
pred.x <- np[,which(colnames(np)%in%all.pollen.taxa.names)]

### organizing and aggregating pollen data
trees <- c("JUGLANSX","FRAXINUX","OSTRYCAR","ULMUS","TILIA","CARYA",
           "FAGUS","TSUGAX","QUERCUS","BETULA",
           'PINUSX',"ACERX","ALNUSX",
           "CYPERACE","PICEAX","ABIES","POPULUS",
           "LARIXPSEU","CUPRESSA") #
other.trees <- c("CASTANEA","PLATANUS","SALIX","LIQUIDAM","TAXUS","NYSSA")#NULL#c()
drop.taxa <- NA#c('other_herbs')

source(file.path('Workflow_Code','utils','taxa_selection.R'))
ten.count <- taxa_selection(trees = trees, other.trees = other.trees,
                            cast.x = pred.x, sites_rm = 0, bigwoods.include = F,
                            all.pollen.taxa.names = all.pollen.taxa.names,
                            prairie.include = T, other.herbs.include = T,
                            other.trees.include = T, drop.taxa = drop.taxa,
                            PFT.do = F)
### Needs to be in same order as calibration so betas match with the right columns
load("threethirds_v2.0.Rdata")
ten.count <- ten.count[,colnames(Y)] #would it be better to sort originally based off of the prediction datasets?

save(x.meta,x.bacon,ten.count,file='prediction.data_IL_IN_v1.Rdata')

rm_nobacon <- which(x.bacon[,49]==x.bacon[,50])

x.meta <- x.meta[-rm_nobacon,]
x.bacon <- x.bacon[-rm_nobacon,]
ten.count <- ten.count[-rm_nobacon,]

save(x.meta,x.bacon,ten.count,file='prediction.data_IL_IN_v2.Rdata')

load('prediction.data_IL_IN_v2.Rdata')

x.meta.s <- x.meta
x.bacon.s <- x.bacon
ten.count.s <- ten.count

load('prediction.data_MN_WI_MI_v1.Rdata')

colnames(x.bacon) <- colnames(x.bacon.s)

x.meta <- rbind(x.meta,x.meta.s)
x.bacon <- rbind(x.bacon,x.bacon.s)
ten.count <- rbind(ten.count,ten.count.s)

plot(x.meta$age_bacon,x.meta$lat)

save(x.meta,x.bacon,ten.count,file='prediction.data_new_v1.Rdata')

n.betas <- 20
n.sites <- length(unique(x.meta$site.name))

dataID <- data.frame(name = sort(rep(unique(x.meta$site.name),n.betas)), ID = 1:(n.sites*n.betas),
                     sigma = rep(0.12,n.sites*n.betas), beta = rep(1:n.betas,n.sites))
write.csv(dataID, file='dataID_bacon_new_v1.csv')

#####
##### Look at plots for full prediction dataset
#####

pdf(file.path(fig.dir,paste('chrono_',Sys.Date(),'.pdf')))
par(mfrow=c(1,1),mar=c(4,4,4,10),xpd=FALSE,pty="s")
plot(new.pollen$age_bacon,new.pollen$lat.x,xlim=c(0,10000),
     col='black',ylab='Latitude',xlab='Median Bacon Age')
points(x$age_bacon,x$lat.x,col='black',cex=1,pch=21,bg='red')
#points(pol_no$age_bacon,pol_no$lat.x,col='blue',cex=.5,pch=19)
legend(11000,46,xpd=TRUE,c('All Pollen \n Samples','Prediction\n Dataset'),
       pch=c(1,19),col=c('black','red'),y.intersp=2)
dev.off()

### making the meta data table
x.meta = x[,c('id','lat.x',"long.x","dataset","site.name","age_bacon")]
colnames(x.meta) <- c('site.id','lat','long','dataset','site.name','age_bacon')

### mapping to look at how many points we have (blue) compared to everything in neotoma (gray)
map('state', ylim=range(new.pollen$lat.x)+c(-2, 2), xlim=range(new.pollen$long.x)+c(-2, 2),main=NA)
points(new.pollen$long.x, new.pollen$lat.x, pch=19, cex=1,col="gray")
points(x.meta$long, x.meta$lat, pch=19, cex=1,col="blue")
legend('topright',c('All Sites','Prediction Sites'),pch=19,col=c('gray','blue'))

### need just bacon draws for parts of analysis in same order as x.meta
x.bacon <- x[,grep(pattern = 'bacon_draw',colnames(x))]

### need just the pollen data to give to organizational function
all.pollen.taxa.names <- colnames(pol_cal_count)[11:length(colnames(pol_cal_count))]
pred.x <- x[,which(colnames(x)%in%all.pollen.taxa.names)]

### organizing and aggregating pollen data
trees <- c("JUGLANSX","FRAXINUX","OSTRYCAR","ULMUS","TILIA","CARYA",
           "FAGUS","TSUGAX","QUERCUS","BETULA",
           'PINUSX',"ACERX","ALNUSX",
           "CYPERACE","PICEAX","ABIES","POPULUS",
           "LARIXPSEU","CUPRESSA") #
other.trees <- c("CASTANEA","PLATANUS","SALIX","LIQUIDAM","TAXUS","NYSSA")#NULL#c()
drop.taxa <- NA#c('other_herbs')

source(file.path('Workflow_Code','utils','taxa_selection.R'))
ten.count <- taxa_selection(trees = trees, other.trees = other.trees,
                            cast.x = pred.x, sites_rm = 0, bigwoods.include = F,
                            all.pollen.taxa.names = all.pollen.taxa.names,
                            prairie.include = T, other.herbs.include = T,
                            other.trees.include = T, drop.taxa = drop.taxa,
                            PFT.do = F)
### IMPORTANT! Needs to be in same order as calibration so betas match with the right columns
load("threethirds_v2.0.Rdata")
ten.count <- ten.count[,colnames(Y)] #would it be better to sort originally based off of the prediction datasets?

### Save MAKE SURE TO CHANGE VERSION NUMBER IF YOU CHANGE ANYTHING
save(x.meta,x.bacon,ten.count,file='prediction.data_v4.Rdata')

######
###### Writing Job Array Script
######
IDs <- unique(x.meta$site.id)
name.list <- list()
site.data <- as.data.frame(matrix(NA,length(IDs),5))

source(file.path('Workflow_Code','utils','test_site.R'))
for(i in 1:length(IDs)){ 
  which_rows <- which(x.meta$site.id == IDs[i])
  x.meta.use <- x.meta[which_rows,]
  
  n.samps <- length(which(x.meta.use$age_bacon<10000&x.meta.use$age_bacon>2000))
  
  test_site(x.meta.use)
  
  ten_count_use = ten.count[which_rows, ]
  
  site.data[i,] <- data.frame(max.age = max(x.meta[which_rows,'age_bacon']),
                              min.age = min(x.meta[which_rows,'age_bacon']),
                              lat = x.meta[which_rows,'lat.x'][1],
                              long = x.meta[which_rows,'long.x'][1],
                              n.samps = n.samps)
  name.list[[i]] <-  x.meta[which_rows,'site.name'][1]
  
}

site.data <- cbind(site.data,unlist(name.list),unlist(IDs))
colnames(site.data) <- c('max.age','min.age','lat','long','n.samps','site.name','site.id')

pp <- site.data[site.data$max.age>8000,]
pp <- pp[pp$min.age<=2000,]
pp <- pp[pp$n.samps>=10,]

site.data[-which(site.data$site.name%in%pp$site.name),]

site.data <- pp

site_keep <- site.data

if(which(table(site.data$site.name)>1)) {
  print('Site Doubled Up.')
  print(table(site_keep$site.name)[which(table(site_keep$site.name)>1)])
}

site.data <- site.data[-which(site.data$site.name=='Kellys Hollow')[1],]
site.data <- site.data[-which(site.data$site.name=='Lily Lake')[2],]

n.sites <- nrow(site.data)
n.betas <- 1

# Means
dataID <- data.frame(name = sort(rep(site.data$site.name,n.betas)),
                     ID = 1:n.sites,
                     sigma = rep(0.03,n.sites*n.betas),
                     beta = rep(NA,n.sites))
write.csv(dataID, file='dataID_bacon_v4_means.csv')

n.betas <- 50
# Beta Draws
dataID <- data.frame(name = sort(rep(site.data$site.name,n.betas)), ID = 1:n.sites,
                     sigma = rep(0.03,n.sites*n.betas), beta = rep(1:n.betas,n.sites))
write.csv(dataID, file='dataID_v5.csv')
