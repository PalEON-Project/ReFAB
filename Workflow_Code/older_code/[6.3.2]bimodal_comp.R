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


source(file.path('Workflow_Code','utils','compare_workInfo.R'))
compare_workInfo(path_to_workInfo = '~/Downloads/work_age/',
                 locn = 'Wintergreen Lake')

