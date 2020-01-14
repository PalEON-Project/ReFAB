library(dtwclust)

agb.list <- lapply(biomassCI,FUN = function(x) x[2,])

names(agb.list) <- name.keep

agb.list <- agb.list[-c(22,27,73)]


length(agb.list)


clusters <- tsclust(agb.list,type = 'hierarchical',k=3)

clusters@cluster

dend <- as.dendrogram(clusters)

library(dendextend)
dend <- color_labels(dend,k=3, col = c(1,2,3))

pdf('cluster_3.pdf',height=20,width = 18)
layout(matrix(c(1,2,3,3,4,4,5,5),4,2,byrow=T))

plot(dend) # selective coloring of branches AND labels :)

library(maps)

plot(unlist(long),unlist(lat),col=clusters@cluster,pch=19,cex=2)
map('state',add=T)

#par(mfrow=c(3,1))
for(i in 1:3){
  plot_me <- do.call(cbind,agb.list[which(clusters@cluster==i)])
  matplot(plot_me,xlim=c(100,0),typ='l',col=i,ylab = 'Biomass (Mg/ha)',xlab='Time',lwd=2)
  points(rowMeans(plot_me),type = 'l',lwd=4,col='blue')
}
dev.off()

cluster_all <- numeric(nrow(x.meta))
for(ii in 1:nrow(x.meta)){
  if(any(which(x.meta$site.name[ii]==names(agb.list)))) cluster_all[ii] <- clusters@cluster[which(x.meta$site.name[ii]==names(agb.list))]  
}

pollen_clusters <- list()
for(i in 1:3){
  pollen_clusters[[i]] <- ten.count[which(cluster_all==i),]
}

pdf('pollen_clusters_3.pdf')
library(vioplot)
par(mfrow=c(1,1))
vioplot(pollen_clusters[[1]],at = 1:22,las=3,col=adjustcolor(1,alpha.f = .5),ylim=c(0,1000))
vioplot(pollen_clusters[[2]],at = 1:22 +.1,las=3,col=adjustcolor(2,alpha.f = .5),add=T)
vioplot(pollen_clusters[[3]],at = 1:22 +.2,las=3,col=adjustcolor(3,alpha.f = .5),add=T)

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

library(gridExtra)
grid.arrange(p1,p2,p3,nrow=1)
dev.off()


pdf('average_pollen_ts.pdf',height=20,width = 18)
library(plotrix)
pol_ages <- list()
layout(matrix(c(1,1,1,4,2,2,2,4,3,3,3,0),3,4,byrow=T))
for(i in 1:3){
  pol_ages[[i]] <- x.meta$age_bacon[which(cluster_all==i)]


breaks <- c(seq(0,10000,100),15000)
cuts <- cut(pol_ages[[i]],breaks,labels = 1:(length(breaks)-1),include.lowest = T)

time_bin <- matrix(NA,100,22)
for(tt in 1:100){
 time_bin[tt,] <-  colMeans(prop.table(as.matrix(pollen_clusters[[i]][which(cuts == tt),]),margin = 1))
}

barplot(t(time_bin),col = rainbow(22),main = paste('Cluster',i))
}
plot.new()
legend('center',colnames(pollen_clusters[[1]]),col=rainbow(22),pch=19,cex=2)
dev.off()
