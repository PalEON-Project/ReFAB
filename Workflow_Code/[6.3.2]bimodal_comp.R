load("twothirds_v1.0.Rdata")
load('10foldLiks.Rdata')

outLik.mat <- do.call(rbind,outLik)
rownames(outLik.mat) <- as.vector(sets10)

bimodal_sites <- c(3,98,62,100,73,92,63,55, 76, 4, 87, 79, 51)

locs <- cbind(cast.x$lat,cast.x$long)[-sites_rm,]
near_site <- numeric(length(bimodal_sites))

for(i in 1:length(bimodal_sites)){
  core_site <- locs[bimodal_sites[i],]
  d <- rdist(matrix(core_site, ncol=2), as.matrix(locs))
  near_site[i] <- order(d)[2]
}

prop_all <- Y/rowSums(Y)

pdf('Near_v_bimodal.pdf')
par(mfrow=c(2,2))
for(i in 1:length(bimodal_sites)){
  prop_range_calc <- c(as.numeric(prop_all[near_site[i],]),as.numeric(prop_all[bimodal_sites[i],]))
  plot(as.numeric(prop_all[near_site[i],]),typ='o',ylim = range(prop_range_calc),
       main = paste('Bimod Site',bimodal_sites[i],"v Near Site",near_site[i]),
       ylab = 'Pollen Prop', xlab = NA,xaxt = 'n')
  axis(side = 1,at = 1:22,labels = colnames(Y),las=2,cex=1)
  points(as.numeric(prop_all[bimodal_sites[i],]),type = 'o',col='red')
  legend('topright',
         legend = c(paste('sum Near site =', sum(Y[near_site[i],])),
         paste('sum Bimod site =', sum(Y[bimodal_sites[i],]))),
         pch = c(19,19),col=c('black','red'))
  range_calc <- c(outLik.mat[which(rownames(outLik.mat)==near_site[i]),],outLik.mat[which(rownames(outLik.mat)==bimodal_sites[i]),])
  plot(outLik.mat[which(rownames(outLik.mat)==near_site[i]),],typ='l',
       ylim=range(range_calc), ylab = 'Likelihood',xlab = 'Biomass',
       main = "Likelihoods Compared")
  points(outLik.mat[which(rownames(outLik.mat)==bimodal_sites[i]),],typ='l',col='red')
  }
dev.off()

out.keep <- list()
for(i in 1:20){
  load(paste0('~/Downloads/workInfo-2/workInfo_',1160+i,
       '_Wintergreen-Lake_Beta_',i,'.Rdata'))
  out.keep[[i]] <- out
}
load('prediction.data.Rdata')

locn <- 'Wintergreen Lake'

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
dim(Y2)

pdf('wintergreen_prop_liks.pdf')
bMax <- 150
par(mfrow=c(2,3))
for(i in 1:(ncol(out)-1)){
  plot(as.numeric(Y2[i,3:24]/sum(Y2[i,3:24])),typ='o', ylab = 'Pollen Prop', xlab = NA, xaxt = 'n')
  points(as.numeric(Y2[i+1,3:24]/sum(Y2[i+1,3:24])),typ='o', ylab = 'Pollen Prop',
         xlab = NA,xaxt = 'n',col='red')
  axis(side = 1,at = 1:22,labels = colnames(Y),las=2,cex=1)
  
  plot(seq(5, bMax-5, by = 2),
       exp(out.keep[[b]][,i]-max(out.keep[[b]][,i]))/-sum(out.keep[[b]][,i]),
       typ='l',main=paste('age',Y2$age_index[i]),ylab='Likelihood')
  for(b in 1:20){
    points(seq(5, bMax-5, by = 2),
           exp(out.keep[[b]][,i]-max(out.keep[[b]][,i]))/-sum(out.keep[[b]][,i]),
           typ='l')
  }
  
  plot(seq(5, bMax-5, by = 2),
       exp(out.keep[[b]][,i+1]-max(out.keep[[b]][,i+1]))/-sum(out.keep[[b]][,i+1]),
       typ='l',main=paste('age',Y2$age_index[i+1]),col='red',ylab='Likelihood')
  for(b in 1:20){
   points(seq(5, bMax-5, by = 2),
          exp(out.keep[[b]][,i+1]-max(out.keep[[b]][,i+1]))/-sum(out.keep[[b]][,i+1]),
          typ='l',col='red')
  }
}
dev.off()
