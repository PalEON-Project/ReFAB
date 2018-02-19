blue       <- col2rgb("blue")
alphablue  <- rgb(blue[1], blue[2], blue[3], 75, max = 255)

biomassCI <- diff.median <-  list()
prob.of.inc <- matrix(NA,62,99)
all.samps<- numeric(100)
nItsSave = 10000
hem.is <- out.list <- LS.is <- list()
lat <- long <- all.samps.list <- name.keep <- age.keep <- list()
for(i in 39:length(unique(dataID$name))){
  locn <- as.character(unique(dataID$name)[i])
  locnClean <- gsub(' ', '-', unique(dataID$name)[i])
  site_number = unique(x.meta[x.meta$site.name == locn,1])
  
  x.meta.use <- x.meta[x.meta$site.name == locn,]
  
  source('test_site.R')
  test_site(x.meta.use)
  
  ten_count_use = ten.count[which(x.meta$site.name == locn), ]
  ten_count_use[which(is.na(ten_count_use))] <- 0
  
  Y = as.matrix(ten_count_use)
  hem.is[[i]]<-max(prop.table(ten_count_use,1)[,11]) + max(prop.table(ten_count_use,1)[,19])#sum(ten_count_use[,11])
  
  sample_ages <- x.meta.use$age_bacon
  age_bins <- seq(minAge, maxAge, ageInterval)
  age_index <- as.matrix(as.numeric(
    cut(sample_ages, breaks = age_bins, labels=seq(1:(length(age_bins)-1)))
  ))
  
  tmp <- data.frame(cbind(age_index, Y))
  names(tmp)[1] <- 'age_index'
  
  Y2 <- aggregate(tmp, by = list(tmp$age_index), FUN = sum)
  
  if(!is.null(group) | FALSE){
    Y2 <- Y2[-group.mat[group,],]
  }
  
  Y <- as.matrix(Y2[ , -c(1,2)])
  age_index <- Y2[,1]
  
  LS.is[[i]] <- cbind(age_index,prop.table(Y,1)[,'TSUGAX'] + prop.table(Y,1)[,'FAGUS'])
  
  if(length(Y)>21 & nrow(Y) > 10 &
     max(x.meta[x.meta$site.name == locn,'age_bacon'])>8000 & 
     min(x.meta[x.meta$site.name == locn,'age_bacon'])<2000){
    
    samples.keep <- numeric(300)
    
    for(b in 1:20){
      ID <- dataID[dataID$name==as.character(locn),'ID'][b]
      file_name <- paste0(path_to_samps,'samplesList_workInfo_',ID,'_',locnClean,'_Beta_',b,'.Rda')
      if(!file.exists(file_name)) next()
      load(file_name)
      samples.keep <- rbind(samples.keep, samplesList[1000:2000,])
    }
    
    samplesList <- samples.keep
    
    if(file.exists(file_name)){
      b <- 1
      ID <- dataID[dataID$name==as.character(locn),'ID'][b]
      file_name1 <- paste0(path_to_Info,'workInfo_',ID,'_',locnClean,'_Beta_',b,'.Rda')

      if(file.exists(file_name1)){
        load(file = file_name1)
        out.list[[i]] <- out
      } else{
        out.list[[i]] <- NA
        print(paste('missing',locn))
      }
      
      
      #Takes out modern data estimates where data stop# i.e. 'cut' approach
      #samplesList[,1:min(age_index)] <- NA
      
      not_burn <- 1:nrow(samplesList)
      biomassCI[[i]] <- apply(samplesList[not_burn,1:100],2,quantile,c(0.025,0.5,0.975),na.rm=TRUE)
      names(biomassCI)[i] <- x.meta[x.meta$site.name == locn,'site.id'][1]
      
      all.samps.list[[i]] <- samplesList[not_burn,1:100]
      all.samps <- rbind(all.samps,samplesList[not_burn,1:100])

      diff.median[[i]] <- apply(diff(t(samplesList[not_burn,1:100])),1,quantile,c(0.5),na.rm=TRUE)
      test <- diff(t(samplesList[not_burn,1:100]))
      lat[[i]] <- x.meta[x.meta$site.name == locn,'lat'][1]
      long[[i]] <- x.meta[x.meta$site.name == locn,'long'][1]
      name.keep[[i]] <- locn
      
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

### For Paleon MIP
biomass.mean.df <- data.frame(lat=unlist(lat)[-35],lon=unlist(long)[-35],biomassMean = unlist(lapply(biomassCI,function(x){mean(x[2,1:11])}))[-c(2,35)])
write.csv(biomass.mean.df,file='biomass.means.csv')


biomass.settle <- last.index <-numeric(62)
for(i in 1:62){
  if(!is.null(name.keep[[i]])&i!=35){
    biomass.settle[i] <- x.meta[x.meta$site.name==name.keep[[i]],'SettleBiomassMean'][1]
    if(!is.na(biomass.settle[i])){
      last.index[i] <- biomassCI[[i]][2,1]
    }else{
      last.index[i] <- NA
    }
  }
  }


### Creating correlation plot between settlement PLS estimates and settlement ReFAB estiamtes
last.index <- last.index[-c(which(is.na(last.index)),which(last.index==0))]
biomass.settle <- biomass.settle[-c(which(is.na(biomass.settle)),which(biomass.settle==0))]

pdf('settle.biomass.corr.pdf')
par(mfrow=c(1,1))
plot(biomass.settle,last.index,pch=19,xlab='Settlement Biomass Mean',ylab='Last Prediction Biomass Mean',xlim=c(0,155),ylim=c(0,155))
abline(a=0,b=1)
legend('topleft',paste('cor =', cor(x=biomass.settle,y=last.index)),pch=NA)
dev.off()

### Mapping sites by the last data point to make sure there is no spatial correlation in end of data
map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(unlist(long),unlist(lat), pch=19,
       cex=1.1,lwd=.2,col=rainbow(20,start = 0,end = .8)[unlist(stop.spot)])
legend('topright',as.character(1:20),col=rainbow(20,start = 0,end = .8),pch=rep(19,20))
title('Map Colored by Last Sample Age Index')

#####
##### Calculate derived quantities for plotting #####
#####

bio.quant <- apply(all.samps,2,quantile,c(0.025,0.5,0.975),na.rm=TRUE)
diff.mat.all<- t(diff(t(all.samps[,100:1])))

## Find cutoff for time series by last data point
stop.spot <- list()
for(i in 1:62){
  stop.spot[[i]] <- min(age.keep[[i]],na.rm = TRUE)
}

## Looking for where to cut all.samps matrix -- changes with number of sites
site_count <- 60-1 #number of sites minus 1

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
##### Creating Figure 1 ##### 
##### 

pdf(paste0('average.biomass_',Sys.Date(),'.pdf'))
#quartz()
zones=matrix(c(2,1), ncol=1, byrow=TRUE)
layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
par(mar=c(5,4,1,1))
plot(seq(100,10000,100),save.median.bucket[2,],
     col=colors[data_binned],pch=19,xlim=c(10000,-10),
     ylim=c(0,150),xlab = 'Year Before Present',ylab='Biomass (Mg/ha)')
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
abline(h=150)
ciEnvelope(x = seq(100,10000,100),ylo = save.median.bucket[1,],yhi =save.median.bucket[3,],col=alphablue)
points(seq(100,10000,100),save.median.bucket[2,],
       col=colors[data_binned],pch=19)
#points(seq(100,10000,100),bio.quant[1,], pch = 19)
#points(seq(100,10000,100),bio.quant[3,], pch=19)
par(mar=c(0,4,1,1))
barplot(height=save.how.much[2,], ylim=c(.2, .8), 
        space=0,xlim=c(0,100),ylab = 'Proportion Increasing',
        cex.lab=.6, cex.axis = .5, col='white', border = 'white')

segments(seq(.5,98.5,1),save.how.much[1,],seq(.5,98.5,1),save.how.much[3,],col='darkgrey',lwd=1)
points(seq(.5,98.5,1),save.how.much[2,],pch=19,cex=.5)

abline(h=.5)
dev.off()










### Old code for showing RW and genpareto outputs together
blue       <- col2rgb("blue")
alphablue  <- rgb(blue[1], blue[2], blue[3], 75, max = 255)
orange       <- col2rgb("orange")
alphaorange  <- rgb(orange[1], orange[2], orange[3], 85, max = 255)

pdf('compare.RW.v.genPareto.Sigma0.12.pdf')
par(mfrow=c(2,2))
plot.new()
legend('center',c('GenPareto Sigma = .12','RW'),pch=c(19,19),col=c('darkorange','blue'))
for(i in 1:62){
  locn <- names(how.many)[i]
  site_number = unique(x.meta[x.meta$site.name == locn,1])
  keep.dataset.id <- unique(x.meta[x.meta$site.id==site_number,4])
  
  ten_count_use = ten.count[which(x.meta$site.id == site_number), ]
  
  Y = as.matrix(ten_count_use)
  
  sample_ages <- x.meta[x.meta[,1] == site_number, ]$age_bacon
  age_bins <- seq(minAge, maxAge, ageInterval)
  age_index <- as.matrix(as.numeric(
    cut(sample_ages, breaks = age_bins, labels=seq(1:(length(age_bins)-1)))
  ))
  
  tmp <- data.frame(cbind(age_index, Y))
  names(tmp)[1] <- 'age_index'
  
  Y2 <- aggregate(tmp, by = list(tmp$age_index), FUN = sum)
  
  if(!is.null(group) | FALSE){
    Y2 <- Y2[-group.mat[group,],]
  }
  
  Y <- as.matrix(Y2[ , -c(1,2)])
  age_index <- Y2[,1]
  
  plot(biomassCI[[i]][2,],col='blue',typ='l',xlim=c(100,0),
       ylim=c(0,150),lwd=2,main=names(how.many)[i])
  ciEnvelope(x = 1:100, yhi = biomassCI[[i]][3,],ylo=biomassCI[[i]][1,],col=alphablue)
  lines(biomassCI.save[[i]][2,],col='darkorange',lwd=2)
  ciEnvelope(x = 1:100, yhi = biomassCI.save[[i]][3,],ylo=biomassCI.save[[i]][1,],col=alphaorange)
  rug(x.meta[x.meta[,1]==site_number,]$age_bacon/100,lwd=2)
  
  if(length(out.list[[i]])!=0){
    points(age_index,seq(5, 145, by = 2)[apply(out.list[[i]],2,which.max)])
  }
  
  mean.keep <- x.meta[x.meta$site.name == locn,'SettleBiomassMean'][1]
  sd.keep <- x.meta[x.meta$site.name == locn,'SettleBiomassSD'][1]
  points(0,mean.keep,col='blue',cex=1.5,pch=19)
  segments(x0=0,y0=mean.keep-sd.keep,x1 = 0,y1=mean.keep+sd.keep)
  #rug(control.pts[which(control.pts[,2]%in%keep.dataset.id),]$geo_age/100,lwd=3,col="red")
}
dev.off()



