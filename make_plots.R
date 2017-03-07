

plot_biomass_ts <- function(site_number, biomassCI){
  fig.mat <- matrix(1,27,1)
  fig.mat[1:6,]<-1
  fig.mat[7:27,]<-seq(2,22,1)
  #control.pts<-read.csv(text=getURL('https://raw.githubusercontent.com/PalEON-Project/stepps-baconizing/master/data/bacon_age_markers_v1.csv'))
  
  layout(fig.mat)
  par(mar = c(0, 4, 0, 3), oma = c(5, 0, 4, 0),tcl=.5)
  
  breaks <-  c(seq(0,50,10),seq(75,200,25))
  colors <- rev(terrain.colors(length(breaks)))
  
  
  data_binned <-  cut(biomassCI[2,], c(breaks), include.lowest = FALSE, labels = FALSE)
  
  plot(seq(100,10000,100),biomassCI[2,],
       cex=.1,ylim=c(0,150),xlim=c(10000,-10),ylab="Biomass (Mg / Ha)",
       xlab="Years Before Present",main=NA,, xaxt='n')
  
  title(x.meta[x.meta$site.id==site_number,'site.name'][1],outer=TRUE)
  axis(3,at=seq(0,10000,1000),labels=seq(0,10000,1000),padj=1)
  
  ciEnvelope(seq(100,10000,100),biomassCI[1,],biomassCI[3,],col="gray")
  
  points(seq(100,10000,100),biomassCI[2,],cex=.8,pch=16,col = colors[data_binned])
  
  keep.dataset.id <- unique(x.meta[x.meta$site.id==site_number,4])
  
  rug(x.meta[x.meta[,1]==site_number,]$age_bacon,lwd=2)
  #rug(control.pts[which(control.pts[,1]%in%keep.dataset.id),]$geo_age,lwd=3,col="red")
  
  ten.count.use = ten.count[which(x.meta[,1]==site_number),]
  prop.use <- prop.table(as.matrix(ten.count.use),margin=1)    
  
  for(p in rev(match(names(sort(colMeans(prop.use))),colnames(prop.use)))){
    prop.plot<- cbind(as.vector(x.meta[which(x.meta$site.id==site_number),]$age_bacon),as.matrix(prop.use[,p]))      	
    prop.plot<-prop.plot[order(prop.plot[,1]),]
    plot(x=prop.plot[,1],y=prop.plot[,2],type="l",xlim=c(10000,-10),
         ylim=c(0,max(prop.use[,p])),ylab=NA,yaxt='n', xaxt='n')
    #axis(2,at=signif(max(prop.use)/2),labels=signif(max(prop.use[,p])/2,digits=1))
    ciEnvelope(prop.plot[,1],rep(0,nrow(prop.use)),prop.plot[,2],col="darkblue")
    legend('topleft',colnames(prop.use)[p])
    #legend('topleft',paste(signif(mean(prop.use[,p]),digits=2)),bty="n")
  } 
  
}

biomassCI <- apply(as.matrix(samplesList[[1]][,1:100],samplesList[[2]][,1:100],samplesList[[3]][,1:100]),2,quantile,c(0.025,0.5,0.975))
biomassCI1 <- apply(as.matrix(samples.pred.save[[31]][[1]][,1:100],),2,quantile,c(0.025,0.5,0.975))

default.par<- par()

pdf('norm.50K.iters.pdf')
for(s in 1:31){ 
  print(paste('working on number',s,'of 183',(s/183)*100,'% complete'))
  site_number = unique(x.meta[,1])[s]
  
  ten.count.use = ten.count[which(x.meta$site.id==site_number),]
if(length(ten.count.use)>25*20&	min(x.meta[x.meta[,1]== site_number,]$age_bacon)<1000 & 
   max(x.meta[x.meta[,1]== site_number,]$age_bacon)>9000){
site_number = unique(x.meta[,1])[s]
biomassCI1 <-  apply(as.data.frame(rbind(samples.pred.save[[s]][[1]][,1:100],samples.pred.save[[s]][[2]][,1:100],samples.pred.save[[s]][[3]][,1:100])),2,quantile,c(0.025,0.5,0.975))
plot_biomass_ts(site_number = site_number, biomassCI = biomassCI1)
par(default.par)
par(mfrow=c(3,3))
for(i in sample.int(size = 9,n = 100)){
  plot(samples.pred.save[[s]][[1]][,i],typ='l',main=i,ylim=c(0,150),
       col = rainbow(3)[1])
  points(samples.pred.save[[s]][[2]][,i],typ='l',col= rainbow(3,alpha = 0.8)[2])
  points(samples.pred.save[[s]][[3]][,i],typ='l',col= rainbow(3,alpha = 0.4)[3])
}

}
}

dev.off()