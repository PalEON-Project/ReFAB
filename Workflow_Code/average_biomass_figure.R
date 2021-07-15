

###### This script creates average biomass figure from the main text.
###### Also, makes all site diagnostic pdfs


dataID <- read.csv('dataID_v5.csv') #dataID <- read.csv('dataID_bacon_new_v1.csv') #for original preds dataID <- 
load('prediction.data_v6.Rdata') #load('prediction.data_new_v1.Rdata') #
source(file.path('Workflow_Code', 'utils', 'validation_args.R'))
control.pts <- read.csv(file.path('Data', 'control.pts.csv'))

n.sites <- length(unique(dataID$name))

path_to_samps <- c('~/Downloads/samps_again_FULL/')#c('~/ReFAB/samps_final/')
path_to_Info <- c('~/Downloads/outs_again_FULL/')#c('~/ReFAB/workInfo_final/')

blue       <- col2rgb("blue")
alphablue  <- rgb(blue[1], blue[2], blue[3], 75, max = 255)

biomassCI <- diff.median <-  list()
prob.of.inc <- matrix(NA, n.sites, 99)
all.samps <- numeric(100)
nItsSave <- 250
minAge <- 0
maxAge <- 10000
ageInterval <- 100
big_keep <- hem.is <- out.list <- LS.is <- Y.keep <- list()
lat <- long <- all.samps.list <- info <- name.keep <- age.keep <- site_num_keep <- list()
do.samps = TRUE
#####
##### Load Data #####
#####

for(i in 1:length(unique(dataID$name))){
  locn <- as.character(unique(dataID$name)[i])
  locnClean <- gsub(' ', '-', unique(dataID$name)[i])
  site_number <- site_num_keep[[i]]<- unique(x.meta[x.meta$site.name == locn,1])
  if(locn == 'Lily Lake' | locn == 'Mud Lake' | locn == 'Seidel') next()
  
  x.meta.use <- x.meta[x.meta$site.name == locn,]
  
  source(file.path('Workflow_Code','utils','test_site.R'))
  test_site(x.meta.use)
  
  ten_count_use = ten.count[which(x.meta$site.name == locn), ]
  ten_count_use[which(is.na(ten_count_use))] <- 0
  
  big_keep[[i]] <- cbind(x.meta.use,ten_count_use)
  
  Y = as.matrix(ten_count_use)
  hem.is[[i]]<-max(prop.table(Y,1)[,'TSUGAX']) + max(prop.table(Y,1)[,'FAGUS'])#sum(ten_count_use[,11])
  
  sample_ages <- x.meta.use$age_bacon
  age_bins <- seq(minAge, maxAge, ageInterval)
  age_index <- as.matrix(as.numeric(
    cut(sample_ages, breaks = age_bins, labels=seq(1:(length(age_bins)-1)))
  ))
  
  tmp <- data.frame(cbind(age_index, Y))
  names(tmp)[1] <- 'age_index'
  
  Y2 <- aggregate(tmp, by = list(tmp$age_index), FUN = sum)
  
  Y <- as.matrix(Y2[ , -c(1,2)])
  age_index <- Y2[,1] 
  
  Y.keep[[i]] <- Y2
  
  LS.is[[i]] <- cbind(age_index,prop.table(Y,1)[,'TSUGAX'] + prop.table(Y,1)[,'FAGUS'])
  
  if(length(Y)>21 & nrow(Y) > 10 &
     max(x.meta[x.meta$site.name == locn,'age_bacon'])>8000 & 
     min(x.meta[x.meta$site.name == locn,'age_bacon'])<2000 &do.samps == TRUE){
    
    samples.keep <- numeric(300)
    
    for(b in 1:50){
      ID <- dataID[dataID$name==as.character(locn),'ID'][b]
      file_name <- paste0(path_to_samps,'samplesList_workInfo_',ID,'_',locnClean,'_Beta_',b,'.Rdata')
      if(!file.exists(file_name)) next()
      load(file_name)
      samples.keep <- rbind(samples.keep, samplesList)
    }
    
    samplesList <- samples.keep
    
    if(file.exists(file_name)){
      out.list[[i]] <- list()
      
      for(b in 1:20){
        ID <- dataID[dataID$name==as.character(locn),'ID'][b]
        file_name1 <- paste0(path_to_Info,'workInfo_',ID,'_',locnClean,'_Beta_',b,'.Rdata')
        
        if(file.exists(file_name1)){
          load(file = file_name1)
          
          out.list[[i]][[b]] <- out
          
          
        } else{
          out.list[[i]] <- NA
          print(paste('missing',locn))
        }
      }

      
      #Takes out modern data estimates where data stop# i.e. 'cut' approach
      #samplesList[,1:min(age_index)] <- NA
      
      not_burn <- seq(1,nrow(samplesList),length.out = 250)
      biomassCI[[i]] <- apply(samplesList[not_burn,1:100],2,quantile,c(0.025,0.5,0.975),na.rm=TRUE)
      names(biomassCI)[i] <- x.meta[x.meta$site.name == locn,'site.id'][1]
      
      all.samps.list[[i]] <- samplesList[not_burn,1:100]
      all.samps <- rbind(all.samps,samplesList[not_burn,1:100])

      diff.median[[i]] <- apply(diff(t(samplesList[not_burn,1:100])),1,quantile,c(0.5),na.rm=TRUE)
      test <- diff(t(samplesList[not_burn,1:100]))
      lat[[i]] <- x.meta[x.meta$site.name == locn,'lat'][1]
      long[[i]] <- x.meta[x.meta$site.name == locn,'long'][1]
      name.keep[[i]] <- locn
      
      #consider taking more max_ests from other beta draws
      info[[i]] <-
        data.frame(
          prop.table(Y, 1),
          max_est = colMeans(do.call(rbind,lapply(out.list[[i]],FUN=function(x)seq(5, bMax - 5, by = 2)[apply(x,2,which.max)]))),
          age = age_index,
          site = rep(locn[1], nrow(Y)),
          lat = rep(lat[[i]][1],nrow(Y)),
          lon = rep(long[[i]][1],nrow(Y))
        )
      
      sample_ages <- x.meta[x.meta[,1] == site_number, ]$age_bacon
      age_bins <- seq(minAge, maxAge, ageInterval)
      age_index <- as.matrix(as.numeric(
        cut(sample_ages, breaks = age_bins, labels=seq(1:(length(age_bins)-1)))
      ))
      
      
      
      age.keep[[i]] <- age_index
      for(t in 1:99){
        prob.of.inc[i,t] <- length(which(test[t,]<0))/length(not_burn)
      }
      
      }
  }
  
  
}

save(info,file='info.Rdata')
#save(all.samps.list,lat,long,file = 'refab_for_stability_v3.Rdata')

#####
##### Write out netcdf for biomass reconstructions #####
#####

big_mat <- do.call(rbind,big_keep)

write.csv(big_mat,file = paste0(Sys.Date(),'prediction_samples_ReFAB.csv'))

#####
##### Calculate derived quantities for plotting #####
#####

bio.quant <- apply(all.samps,2,quantile,c(0.025,0.5,0.975),na.rm=TRUE)
diff.mat.all<- t(diff(t(all.samps[,100:1])))

## Find cutoff for time series by last data point
stop.spot <- start.spot <- list()
for(i in 1:n.sites){
  stop.spot[[i]] <- min(age.keep[[i]],na.rm = TRUE)
  start.spot[[i]] <- max(age.keep[[i]],na.rm = TRUE)
}

## Looking for where to cut all.samps matrix -- changes with number of sites
site_count <- 78-1 #number of sites minus 1

## calculation for proportion sites increasing -- barplot
how.much <- matrix(NA,length(not_burn),99)
for(i in 1:99){
  for(r in 1:length(not_burn)){
    div.num <- site_count - length(which(is.na(diff.mat.all[seq(1+r,site_count*length(not_burn),length(not_burn)),i])))
    if(div.num==0) div.num <- site_count
    how.much[r,i] <- length(which(diff.mat.all[seq(1+r,site_count*length(not_burn),length(not_burn)),i]>0))/div.num
  }
}
#calculate the probability of increase for each site then make a boxplot of all of the sites probability of increase at each time.
#proportion of medians that are greater than zero
save.how.much <- apply(how.much,2,quantile,c(0.025,0.5,0.975),na.rm=TRUE)


## calculation for median line over time series -- purple line
median.bucket <- matrix(NA,length(not_burn),100)

for(r in 1:length(not_burn)){
    samps.vec.list <- lapply(all.samps.list,function(x){return(x[r,])})
    means.iter <- do.call(rbind,samps.vec.list)
    median.bucket[r,] <- apply(means.iter, 2, median)#colMeans(means.iter,na.rm = TRUE)#apply(means.iter,2,quantile,probs=c(.5),na.rm =TRUE)
}

save.median.bucket <- apply(median.bucket,2,quantile,c(0.025,0.5,0.975),na.rm =TRUE)
save.median.bucket[2,] <- colMeans(median.bucket)

## biomass breaks
breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)))

## late successional breaks
breaks.LS <- seq(0,.55,.05)
colors.LS.palette <- colorRampPalette(c('#99d8c9','#66c2a4','#41ae76','#238b45','#005824')) #colors found using Rcolorbrewer
colors.LS <- c('lightgrey',colors.LS.palette(length(breaks.LS)))

## binning for colors
data_binned <-  cut(save.median.bucket[2,], c(breaks), include.lowest = FALSE, labels = FALSE)
data_binned_LS <-  cut(unlist(hem.is), c(breaks.LS), include.lowest = TRUE, labels = FALSE)
  
##### 
##### Creating Figure ##### 
##### 

bMax.use <- max( unlist(lapply(biomassCI,function(x) x[2,]))) +10

pdf(paste0('average.biomass_',Sys.Date(),'.pdf'))
#quartz()
zones=matrix(c(2,1), ncol=1, byrow=TRUE)
layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
par(mar=c(5,4,1,1))
plot(seq(100,10000,100),save.median.bucket[2,],
     col=colors[data_binned],pch=19,xlim=c(10000,-10),
     ylim=c(0,bMax.use),xlab = 'Year Before Present',ylab='Biomass (Mg/ha)')
for(i in 1:length(biomassCI)){  
  if(length(biomassCI[[i]])>1){
    #if(hem.is[[i]]>10){
      #stop.spot <- min(age.keep[[i]],na.rm = TRUE)
    #stop.spot <- 1
     # if(stop.spot=='Inf') stop.spot = 1
      interp.LS <- approx(LS.is[[i]][,1],LS.is[[i]][,2],xout=1:100,rule=2)
      data_binned_LS <-  cut(interp.LS$y, c(breaks.LS), include.lowest = TRUE, labels = FALSE)
      
     # points(seq(stop.spot*100,10000,100)+rnorm(1,10,10),biomassCI[[i]][2,stop.spot:100],
      #       col=colors.LS[data_binned_LS],type='b',pch=19,cex=data_binned/10,lwd=data_binned/10)
    time.bin <- seq(100,10000,100)
    for(t in 2:length(time.bin)){
      segments(x0 = time.bin[t-1],y0 = biomassCI[[i]][2,t-1],
               x1 = time.bin[t],y1 = biomassCI[[i]][2,t],
               #col='lightgray',lwd = 1)
               col=colors.LS[data_binned_LS[t]],lwd = data_binned_LS[t])
    }
      
      #}else{
     # lines(seq(100,10000,100),biomassCI[[i]][2,],col='red')
    #green       <- col2rgb("green")
    #alphagreen  <- rgb(green[1], green[2], green[3], alpha = 30, max = 255)
    
    #ciEnvelope(x = seq(100,10000,100),ylo = biomassCI[[i]][1,],yhi =biomassCI[[i]][3,],col=alphagreen)
  
    #}
    
  }
}
abline(h=bMax)
ciEnvelope(x = seq(100,10000,100),ylo = save.median.bucket[1,],yhi =save.median.bucket[3,],col=alphablue)
points(seq(100,10000,100),save.median.bucket[2,],
       col=colors[data_binned],pch=19)
text(x = 9900,y = 175,'B',cex=2)

#points(seq(100,10000,100),bio.quant[1,], pch = 19)
#points(seq(100,10000,100),bio.quant[3,], pch=19)
par(mar=c(0,4,1,1))
barplot(height=save.how.much[2,], ylim=c(.2, .8), 
        space=0,xlim=c(0,100),ylab = 'Proportion Increasing',
        cex.lab=.75, cex.axis = .6, col='white', border = 'white')
text(x = 1,y = .7,'A',cex=2)
segments(seq(.5,98.5,1),save.how.much[1,],seq(.5,98.5,1),save.how.much[3,],col='darkgrey',lwd=1)
points(seq(.5,98.5,1),save.how.much[2,],pch=19,cex=.5)

abline(h=.5)
dev.off()


#################################### END FIGURE 
source(file.path('Workflow_Code','utils','site_diag.R'))

#### Site Diagnostic for original 62 prediction sites
for(i in 1:length(unique(dataID$name))){
  locn <- as.character(unique(dataID$name)[i])
  if(locn == 'Lily Lake' | locn == 'Mud Lake' | locn == 'Seidel') next()
  site_diag(bMax = bMax, locn = locn, x.meta = x.meta,
            minAge = 0, maxAge = 10000, 
            ageInterval = 100, path_to_samps = path_to_samps,
            path_to_Info = path_to_Info, control.pts = control.pts,
            ten.count = ten.count, dataID = dataID, out.dir = '~/ReFAB/sites_diag_again_final/')
}
