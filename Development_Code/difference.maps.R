


biomass.mat <- matrix(NA, 80, 102)
nItsSave = 10000
for(i in 1:62){
  locn <- names(how.many)[i]
  site_number = unique(x.meta[x.meta$site.name == locn,1])
  ten_count_use = ten.count[which(x.meta$site.id == site_number), ]
  Y = as.matrix(ten_count_use)
  if(length(Y)>21 & nrow(Y) > 15 &
     max(x.meta[x.meta$site.name == locn,'age_bacon'])>8000 & 
     min(x.meta[x.meta$site.name == locn,'age_bacon'])<2000){
    if(file.exists(paste0('~/ReFAB/MCMC Samples/samplesList_',locn,'.Rda'))){
      load(file = paste0('~/ReFAB/MCMC Samples/samplesList_',locn,'.Rda'))
      biomass.mat[i,1:100] <- as.numeric(apply(samplesList[round(nItsSave/5):nItsSave,1:100],2,quantile,c(0.5)))
      biomass.mat[i,101:102] <- as.numeric(x.meta[x.meta$site.name == locn,c('lat',"long")][1,])
    }
  }
}

biomass_meds <- lapply(biomassCI,FUN=function(x) c(x[2,],lat))

for(i in 1:80){
  biomass.mat[i,] <- c(biomass_meds[[i]],lat[[i]],long[[i]])
}

do.call(rbind,biomass_meds)

breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)))

x <- biomass.mat[,1]
data_binned <-  cut(x, c(breaks), include.lowest = FALSE, labels = FALSE)

map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(biomass.mat[,102],biomass.mat[,101],pch=19,col=colors[data_binned])

breaks <- c(-100,-30,-10,-2.5,2.5,10,30,100)
colors <- colorRampPalette(c('red',"white",'blue'))(length(breaks)-1)

#quartz()
pdf(paste0('difference.maps.',Sys.Date(),'.pdf'))
par(mfrow=c(3,4), oma = rep(.1,4))
b <- 10
for(r in rev(seq(10,100,b))){ #just thousand year time bins
    if(r-b==0) b = b-5
    x <- biomass.mat[,r] - biomass.mat[,(r-b)]
    data_binned <-  cut(x, c(breaks), include.lowest = FALSE, labels = FALSE)
      map('state', xlim=c(-97,-83), ylim=c(41.75,49),mar=rep(0,4))
      points(biomass.mat[,102],biomass.mat[,101],pch=21,bg=colors[data_binned],col='lightgray',cex=1.25)
      title(paste(r*100, '-', (r-b)*100))
}
plot.new()
legend('center',legend=as.character(breaks[1:length(breaks)-1]),pch=rep(19,length(breaks)-1),col=colors)
dev.off()

library(corrplot)

#quartz()
#pdf(paste0('correlation.maps1',Sys.Date(),'.pdf'))
par(mfrow=c(3,4), oma = rep(.1,4))
for(i in rev(seq(10,100,10))){
  c <- cor(t(biomass.mat[,1:i]))
  #corrplot(c,order='AOE')
  
  map('state', xlim=c(-97,-83), ylim=c(41.75,49),mar=rep(0,4))
  points(biomass.mat[,102],biomass.mat[,101],pch=21,bg='gray',col='gray')
  points(biomass.mat[which(rowSums(c)<0),102],
         biomass.mat[which(rowSums(c)<0),101],
         pch=21,bg='black',col='black')
  title(paste(i*100,':', 100))
  
  #matplot(t(biomass.mat[which(rowSums(c)<0),(i-10):i]),typ='l',col='black',lwd=.5,ylim=c(0,150))
  #matplot(t(biomass.mat[-which(rowSums(c)<0),(i-10):i]),typ='l',add = TRUE,col='red',lwd=.5)
  plot(x = 1:i,y=colMeans(biomass.mat[-which(rowSums(c)<0),1:i]),ylim=c(0,150),xlim=rev(range(1:i)),col='white')
  lines(x = 1:i,y=colMeans(biomass.mat[-which(rowSums(c)<0),1:i]),lwd=5,col='gray')
  lines(x = 1:i,y=colMeans(biomass.mat[which(rowSums(c)<0),1:i]),lwd=5)
  
}
dev.off()


P <- rda(biomass.mat[,1:100])
plot(P)
pdf(paste0('pca.scores',Sys.Date(),'.pdf'))
x <- scores(P)$sites[,1]
breaks <- quantile(x,c(.975,.75,.625,.375,.25,.025))
colors <- colorRampPalette(c('red',"white",'blue'))(length(breaks)-1)
data_binned <-  cut(x, c(breaks), include.lowest = FALSE, labels = FALSE)
map('state', xlim=c(-97,-83), ylim=c(41.75,49),mar=rep(0,4))
points(biomass.mat[,102],biomass.mat[,101],pch=21,bg=colors[data_binned],col='lightgray')
dev.off()
  




