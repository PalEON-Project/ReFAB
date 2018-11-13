
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
##### Add new baconized sites
#####

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
### Needs to be in same order as calibration so betas match with the right columns
load("~/ReFAB/2018-11-07twothirds.calibration.data.Rdata")
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
  test_site(x.meta.use)
  
  ten_count_use = ten.count[which_rows, ]
  
  site.data[i,] <- data.frame(max.age = max(x.meta[which_rows,'age_bacon']),
                              min.age = min(x.meta[which_rows,'age_bacon']),
                              lat = x.meta[which_rows,'lat'][1],
                              long = x.meta[which_rows,'long'][1],
                              n.samps = length(which_rows))
  name.list[[i]] <-  x.meta[which_rows,'site.name'][1]
  
}

site.data <- cbind(site.data,unlist(name.list),unlist(IDs))
colnames(site.data) <- c('max.age','min.age','lat','long','n.samps','site.name','site.id')

site_keep_bacon <- site.data

if(which(table(site.data$site.name)>1)) {
  print('Site Doubled Up.')
  print(table(site_keep$site.name)[which(table(site_keep$site.name)>1)])
}

n.sites <- nrow(site.data)
n.betas <- 1

# Means
dataID <- data.frame(name = sort(rep(site.data$site.name,n.betas)), ID = 1:n.sites,
                     sigma = rep(0.12,n.sites*n.betas), beta = rep(NA,n.sites))
write.csv(dataID, file='dataID_bacon_v4_means.csv')

n.betas <- 20
# Beta Draws
dataID <- data.frame(name = sort(rep(site.data$site.name,n.betas)), ID = 1:n.sites,
                     sigma = rep(0.12,n.sites*n.betas), beta = rep(1:n.betas,n.sites))
write.csv(dataID, file='dataID_bacon_v4.csv')



