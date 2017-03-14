#####
##### Trace Plots #####
#####

for(t in 1:100){
plot(samplesList[[1]][,t],typ='l',main=t,ylim=c(0,150),col=rainbow(3,alpha = 1)[1])
points(samplesList[[2]][,t],typ='l',col=rainbow(3,alpha = 0.6)[2])
points(samplesList[[3]][,t],typ='l',col=rainbow(3,alpha = 0.6)[3])

}

for(i in 101){
  plot(samplesList[[1]][,i],typ='l',main='SIGMA',
       ylim=c(0,max(c(samplesList[[1]][,i],samplesList[[2]][,i],samplesList[[3]][,i]))),col=rainbow(3,alpha = 1)[1])
  points(samplesList[[2]][,i],typ='l',col=rainbow(3,alpha = 0.6)[2])
  points(samplesList[[3]][,i],typ='l',col=rainbow(3,alpha = 0.6)[3])
}

#####
##### Time Series #####
#####

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
       xlab="Years Before Present",main=NA, xaxt='n')
  
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

plot_biomass_ts(site_number=unique(x.meta$site.id)[182], biomassCI=biomassCI[[182]])
#####
##### biomassCI exploration #####
#####

##### biomass CI is ordered by unique(x.meta$site.id)

pca=princomp(biomass_dat_est[,c(5:19,21:24)]) 
summary(pca) 
biplot(pca)

#####
##### Second Derivative #####
#####

calc.second.deriv <- function(biomassCI,h,second.deriv){
  for(i in 1:length(biomassCI)){
    if(length(biomassCI[[i]])>10){
      T <- dim(biomassCI[[i]])[2] - h
      t <- h + 1 
      biomassCI[[i]]<-(biomassCI[[i]]) #log or no log here?
      second.deriv[[i]]<-sum(((biomassCI[[i]][2,(t:T)+h]-2*biomassCI[[i]][2,(t:T)]+biomassCI[[i]][2,(t:T)-h])/((h*100)^2))^2)
      
    }else{
      second.deriv[[i]]<-NA
    }
  }
  return(second.deriv)
}

pdf(paste('second.deriv.map',Sys.Date(),'.pdf'))
for(i in c(1,2,5,10,20)){
  h=i
  second.deriv<-calc.second.deriv(biomassCI=biomassCI,h=h,second.deriv=list())
  
  #hist(unlist(second.deriv))
  names(second.deriv) <- unique(x.meta[,1])[1:182]
  second.deriv.unlist <- unlist(second.deriv)
  
  second.deriv.unlist <- na.omit(second.deriv.unlist)
  
  breaks <-  c(0,quantile(probs=c(.05,.25,.5,.75,.975),second.deriv.unlist,na.rm=TRUE),max(second.deriv.unlist)+1)
  
  long.keep <- list()
  lat.keep <- list()
  for(i in 1:length(second.deriv.unlist)){
    long.keep[[i]] <- x.meta[x.meta[,1]==names(second.deriv.unlist)[i],'long'][1]
    lat.keep[[i]] <- x.meta[x.meta[,1]==names(second.deriv.unlist)[i],'lat'][1]
  }
  
  colors <- colorRampPalette(c("white","yellow",'green','blue'))(length(breaks)-1)
  
  data_binned <-  cut(second.deriv.unlist, c(breaks), include.lowest = FALSE, labels = FALSE)
  
  map('state', xlim=c(-98,-81), ylim=c(41.5,50))
  points(long.keep, lat.keep, pch=19,
         cex=1.3, col=colors[data_binned],lwd=.2)
  title(paste("Second Derivative Sum Point Estimates h = ",h*100))
  
  plotInset(-90,47,-82.5,50,
            expr={
              hist(data_binned,col=colors,xaxt="n",xlab=NA,
                   ylab=NA,main=NA,cex.lab=.5,cex.axis=.5,breaks=0:length(unique(data_binned)))
              axis(side=1,signif(breaks,digits=2),at = 0:length(unique(data_binned)),cex.axis = .5,las=2,line=0)
              mtext(side = 1, "Squared Sum", line = 1.5,cex=.5)
              mtext(side = 2, "Frequency", line = 1.7,cex=.5)
            })
  
}
dev.off()



