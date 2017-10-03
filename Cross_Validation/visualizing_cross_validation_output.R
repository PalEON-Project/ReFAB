sigma.vals <- c(.01,.03,.09,.27,.81)

# load in data for all sites, per Ann's original code
if(!file.exists('allPredData.Rda'))
  source('prep_data.R') 
load('allPredData.Rda')

site.names <- unique(x.meta$site.name)
how.many <- list()

for(i in 1:length(site.names)){ #44:length(site.names)
  locn <- site.names[i]
  site_number = unique(x.meta[x.meta$site.name == locn,1])
  ten_count_use = ten.count[which(x.meta$site.id == site_number), ]
  
  Y = as.matrix(ten_count_use)
  if(length(Y)>21 & nrow(Y) > 15 &
     max(x.meta[x.meta$site.name == locn,'age_bacon'])>8000 & 
     min(x.meta[x.meta$site.name == locn,'age_bacon'])<2000){
    
    how.many[[i]]<- locn
  }
}

names(how.many) <- site.names[1:177]
how.many <- unlist(how.many)

sigma.vals <- c(seq(.03,.27,.03),.81)

log.mat <- site.rm.mat <- array(NA,dim=c(10,length(sigma.vals),length(how.many)))

max.list <- min.list <- list()

pdf('Group_sum_by_site.pdf')
par(mfrow=c(3,3))
for(i in 1:length(how.many)){

locnClean <- gsub(' ', '-', names(how.many)[i])

locn <- names(how.many)[i]
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


max.val <- max(group.mat)
max.list[[i]] <- which(group.mat==max.val,arr.ind = TRUE)[1]

min.list[[i]] <- which(group.mat==1,arr.ind = TRUE)[1]

for(g in 1:10){
for(s in 1:length(sigma.vals)){

    if(file.exists(file.path('~/Downloads/log.prob',paste0('log.prob.',locnClean, 'Sigma',
                            sigma.vals[s], 'Group', g,'.Rdata')))){
      load(file.path('~/Downloads/log.prob',paste0('log.prob.',locnClean, 'Sigma',
                                                   sigma.vals[s], 'Group', g,'.Rdata')))
      print('loading')
      #site.rm.mat[g,] <- log.prob.list$samples.rm
      log.mat[g,s,i] <- sum(log.prob.list$log.prob.data)

    }

  }
}
if(file.exists(file.path('~/Downloads/log.prob',paste0('log.prob.',locnClean, 'Sigma',
                                                       sigma.vals[s], 'Group', g,'.Rdata')))){
#plot(sigma.vals,colSums(log.mat[,,i],na.rm = TRUE), pch=19,
 #    cex=2, main =locnClean)
}

#log.mat[c(max.list[[i]],min.list[[i]]),,i] <- NA


}
dev.off()

length(log.mat[which(is.na(log.mat))])



log.list <- numeric(length(sigma.vals))
for(i in 1:length(sigma.vals)){
  log.list[i]<- sum(log.mat[,i,],na.rm = TRUE) 
}

matrix(log.mat)
boxplot(log.mat[,1,])


which(is.na(log.mat))
###all sites plot that doesn't include first and last data
#new.sig.vals <- seq(.01,.09,by = .01) #take out first and last




pdf('Sum_all_sites_all_data.pdf')
par(mfrow=c(1,1))
plot(sigma.vals, log.list, pch=19, cex=2, main = 'Sum All Sites')
dev.off()

locn <- 'Wolsfeld Lake'
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

pdf(paste0(locnClean,'.fixed.sigma.log.prob.pdf'))

plot(sigma.vals,colSums(log.mat[,,2],na.rm = TRUE),pch=19,cex=2)
dev.off()


pdf('Wolsfeld_cross_validation_diag.pdf')
par(mfrow=c(1,1))
boxplot(log.mat[,,39],xlab='sigma',ylab='log likelihood')
boxplot(t(log.mat[,,39]),xlab='group',ylab='log likelihood')
load("~/Downloads/wols/samplesList_workInfo_Wolsfeld-LakeSigma0.81Group7.Rda.Rda")
samplesList.7.5<-samplesList
load("~/Downloads/wols/samplesList_workInfo_Wolsfeld-LakeSigma0.27Group7.Rda.Rda")
samplesList.7.4<-samplesList
load("~/Downloads/wols/samplesList_workInfo_Wolsfeld-LakeSigma0.81Group10.Rda.Rda")

plot(colMeans(samplesList.7.4[,1:100]),pch=19)
points(colMeans(samplesList.7.5[,1:100]),col='red',pch=19)
points(colMeans(samplesList[,1:100]),col='blue',pch=19)
legend('topleft',c('Group 7 Sigma 0.27', 'Group 7 Sigma 0.81', 'Group 10 Sigma 0.81'),pch=c(19,19,19),col=c('black','red','blue'))
#rug(age_index)
rug(Y2$age_index[group.mat[7,]],col='red',lwd=2)
rug(Y2$age_index[group.mat[10,]],col='blue',lwd=2)

par(mfrow=c(3,3))
for(i in 1:100){
  plot(samplesList.7.4[,i],typ='l',main=i)
  lines(samplesList[,i],col='blue')
  lines(samplesList.7.5[,i],col='red')
}
dev.off()

out.keep<-array(NA,dim=c(71,4,10,5))
for(g in 1:10){
  for(s in 1:5){
load(paste0('~/Downloads/left.out.wols/left.out_workInfo.Wolsfeld-LakeSigma',sigma.vals[s],'Group',g,'.Rda'))
print(paste0('~/Downloads/left.out.wols/left.out_workInfo.Wolsfeld-LakeSigma',sigma.vals[s],'Group',g,'.Rda'))
        out.keep[,,g,s] <- out
  }
}

plot(colSums(out.keep[,1,1,]))
for(i in 1:4){
  for(g in 1:10){
    lines(colSums(out.keep[,i,g,]))
  }
}

log.keep<-array(NA,dim=c(4,10,5))
for(g in 1:10){
  for(s in 1:5){
    load(paste0('~/Downloads/log.probs.wols/log.prob.Wolsfeld-LakeSigma',sigma.vals[s],'Group',g,'.Rdata'))
    print(paste0('~/Downloads/log.probs.wols/log.prob.Wolsfeld-LakeSigma',sigma.vals[s],'Group',g,'.Rdata'))
    log.keep[,g,s] <- log.prob.list$log.prob.data
  }
}

diff.keep<-array(NA,dim=c(4,10,4))
plot(sigma.vals,log.keep[1,1,],typ='l',ylim=range(log.keep))
for(g in 1:10){
  for(i in 1:4){
    lines(sigma.vals,log.keep[i,g,])
    if(g==8 & i==3){
      lines(sigma.vals,log.keep[i,g,],col='blue',lwd=3)
    }
    if(g==7 & i==4){
      lines(sigma.vals,log.keep[i,g,],col='red',lwd=3)
    }
    diff.keep[i,g,]<-diff(log.keep[i,g,])
  }
}

diff.keep[which.max(diff.keep)]

#data 3,group 7 #age 60
#data 2,group 3 #age 11

samps.keep <- array(NA,dim=c(200,300,7,6))
out.keep <- age_index_save <- list()
name.keep <- c('Cub-Lake','Gass-Lake','Kirchner-Marsh','Penegor-Lake',
               'Tower-Lake','Wintergreen-Lake','Lake-Mendota') 
sigma.vals <- c(.01,.03,.09,.27,.81,2.43)
for(g in 1:7){
  locn <- gsub('-',' ',name.keep[g])
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
  load(paste0('~/Downloads/workinfo.all/workInfo_',name.keep[g],'Sigma0.81Group','.Rda'))
  out.keep[[g]] <- out
  age_index_save[[g]] <- Y2[,1]
  
  for(s in 1:6){
    load(paste0('~/Downloads/all.samps/samplesList_workInfo_',name.keep[g],'Sigma',sigma.vals[s],'Group','.Rda.Rda'))
    samps.keep[,,g,s] <- samplesList
  }
}

pdf('seven_lakes_six_sigmas.pdf')
par(mfrow=c(6,1),mar=rep(0,4))
for(g in 1:7){
  for(s in 1:6){
    breaks <-  c(seq(0,50,10),seq(75,200,25))
    colors <- rev(terrain.colors(length(breaks)))
    
    #browser()
    
    bio.quants <- apply(samps.keep[,1:100,g,s],2,quantile,c(0.025,0.5,0.975))
    
    data_binned <-  cut(rev(bio.quants[2,]), c(breaks), include.lowest = FALSE, labels = FALSE)
    
    plot(bio.quants[2,], pch=19, main = NA,ylim=c(0,155),col='white',xlim=c(100,0))
    ciEnvelope(x=1:100,ylo = bio.quants[1,],yhi = bio.quants[3,],col = 'grey')
    points(bio.quants[2,],cex=1.1,pch=16,col = rev(colors[data_binned]))
    
    legend('top',paste(name.keep[g],sigma.vals[s]),cex = .8)
    
    points(age_index_save[[g]],seq(5, bMax-5, by = 2)[apply(out.keep[[g]],2,which.max)])
  }
}
dev.off()




