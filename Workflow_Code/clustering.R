library(dtwclust)
library(tidyverse)
library(maps)
library(dendextend)
library(vioplot)
library(ggfortify)
library(gridExtra)
library(plotrix)

#####
##### load biomass and metadata
#####

agb.mat <- read.csv('median_biomass_plus_meta.csv')
agb.list <- split(x = agb.mat[,1:100],f = agb.mat$name) %>%
            lapply(.,as.numeric)

lon <- agb.mat$lon[order(agb.mat$name)]
lat <- agb.mat$lat[order(agb.mat$name)]


#####
##### Make clusters
#####
clusters <- tsclust(agb.list,type = 'hierarchical',k=3)

clusters@cluster

dend <- as.dendrogram(clusters)

dend <- color_labels(dend,k=3, col = c(3,1,2))

#####
##### Make cluster dendrogram and map
#####

pdf('cluster_3.pdf',height=20,width = 18)
layout(matrix(c(1,2,3,3,4,4,5,5),4,2,byrow=T))

plot(dend) # selective coloring of branches AND labels :)

plot(lon,lat,col=clusters@cluster,pch=19,cex=2)
map('state',add=T)
text(lon,lat+.2, labels = names(clusters@cluster),cex=.5)

#par(mfrow=c(3,1))
for(i in 1:3){
  plot_me <- do.call(cbind,agb.list[which(clusters@cluster==i)])
  matplot(plot_me,xlim=c(100,0),typ='l',col=i,ylab = 'Biomass (Mg/ha)',xlab='Time',lwd=2)
  points(rowMeans(plot_me),type = 'l',lwd=4,col='blue')
}
dev.off()

#####
##### divide pollen into cluster groups
#####

pollen_all <- read.csv('refab_pollen_counts.csv')

cluster_all <- numeric(nrow(pollen_all))
for(ii in 1:nrow(pollen_all)){
  if(any(which(pollen_all$site.name[ii]==names(agb.list)))) cluster_all[ii] <- clusters@cluster[which(pollen_all$site.name[ii]==names(agb.list))]  
}

pollen_clusters <- list()
for(i in 1:3){
  pollen_clusters[[i]] <- as.data.frame(prop.table(as.matrix(pollen_all[which(cluster_all==i),2:23]),margin = 1))
}

#####
##### violin plots
#####

pdf('pollen_clusters_3.pdf')
par(mfrow=c(1,1))
vioplot(pollen_clusters[[1]],at = 1:22,las=3,col=adjustcolor(1,alpha.f = .5),ylim=c(0,1))
vioplot(pollen_clusters[[2]],at = 1:22 +.1,las=3,col=adjustcolor(2,alpha.f = .5),add=T)
vioplot(pollen_clusters[[3]],at = 1:22 +.2,las=3,col=adjustcolor(3,alpha.f = .5),add=T)

#####
##### PCAs
#####

pr1 <- prcomp(pollen_clusters[[1]])
p1 <-
  autoplot(
    pr1,
    loadings = T,
    colour = 'black',
    loadings.colour = 'blue',
    loadings.label.colour = 'blue',
    loadings.label = TRUE
  ) + ggtitle('Cluster 1')

pr2 <- prcomp(pollen_clusters[[2]])
p2 <-
  autoplot(
    pr2,
    loadings = T,
    colour = 'red',
    loadings.colour = 'blue',
    loadings.label.colour = 'blue',
    loadings.label = TRUE
  ) + ggtitle('Cluster 2')

pr3 <- prcomp(pollen_clusters[[3]])
p3 <-
  autoplot(
    pr3,
    loadings = T,
    colour = 'green',
    loadings.colour = 'blue',
    loadings.label.colour = 'blue',
    loadings.label = TRUE
  ) + ggtitle('Cluster 3')

grid.arrange(p1,p2,p3,nrow=1)
dev.off()

#####
##### average time series by cluster
#####

pdf('average_pollen_ts.pdf',height=20,width = 18)
pol_ages <- list()
layout(matrix(c(1,1,1,4,2,2,2,4,3,3,3,0),3,4,byrow=T))
for(i in 1:3){
  pol_ages[[i]] <- pollen_all$age_bacon[which(cluster_all==i)]


breaks <- c(seq(0,10000,100),15000)
cuts <- cut(pol_ages[[i]],breaks,labels = 1:(length(breaks)-1),include.lowest = T)

time_bin <- matrix(NA,100,7)
for(tt in 1:100){
 time_bin[tt,] <-  colMeans(pollen_clusters[[i]][which(cuts == tt),c('PINUSX','prairie','QUERCUS','TSUGAX','PICEAX','ACERX','FAGUS')])
}

barplot(t(time_bin),space = 0,col = rainbow(7),main = paste('Cluster',i),xlim=c(100,0))
axis(side = 1, at = 0:100, labels=seq(0,10000,100))
}
plot.new()
legend('center',c('PINUSX','prairie','QUERCUS','TSUGAX','PICEAX','ACERX','FAGUS'),col=rainbow(7),pch=19,cex=2)

for(i in 1:3){
  
  pol_ages[[i]] <- pollen_all$age_bacon[which(cluster_all==i)]
  
  breaks <- c(seq(0,10000,100),15000)
  cuts <- cut(pol_ages[[i]],breaks,labels = 1:(length(breaks)-1),include.lowest = T)
  
  
  time_bin <- array(NA,dim=c(100,7,3))
  for(tt in 1:100){
    time_bin[tt,,] <-  t(apply(pollen_clusters[[i]][which(cuts == tt),c('PINUSX','prairie','QUERCUS','TSUGAX','PICEAX','ACERX','FAGUS')],2,quantile,c(.025,.5,.975),na.rm=T))
  }
  par(mfrow=c(7,1),mar = rep(3,4))
  
  for(s in 1:7){
    plot((time_bin[,s,2]),col = rainbow(7)[s],main = c('PINUSX','prairie','QUERCUS','TSUGAX','PICEAX','ACERX','FAGUS')[s],xlim=c(100,0),typ='l',lwd=3,ylim=range(time_bin[,s,]))
    if(s ==1)  {mtext(paste('Cluster',i),side = 3,at = c(80,1),cex=3)}
    lines((time_bin[,s,1]))
    lines((time_bin[,s,3]))
  }
  
}

dev.off()

