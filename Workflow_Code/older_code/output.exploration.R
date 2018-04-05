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

name.vec <- seq(100,10000,100)
only.means.all<-lapply(biomassCI,function(x){return(x[seq(2,299,3)])})
get.lat.long<-matrix(0,183,2)
for(i in 2:183){
  get.lat.long[i,1]<- unique(x.meta[x.meta[,1]==unique(x.meta[,1])[i],c(3)])
  get.lat.long[i,2]<- unique(x.meta[x.meta[,1]==unique(x.meta[,1])[i],c(2)])
}


names(only.means.all)<-unique(x.meta$site.id)[1:182]

breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)))

pdf(paste0('pred.points.map',Sys.Date(),'.pdf'))
for(r in seq(1,99,1)){
  
  only.means <- unlist(lapply(only.means.all,function(x){return(x[r])}))
  
  data_binned <-  cut(only.means, c(breaks), include.lowest = FALSE, labels = FALSE)
  
  long.keep <- list()
  lat.keep <- list()
  for(i in 1:length(only.means)){
    long.keep[[i]] <- x.meta[x.meta[,1]==names(only.means)[i],'long'][1]
    lat.keep[[i]] <- x.meta[x.meta[,1]==names(only.means)[i],'lat'][1]
  }
  
  
  map('state', xlim=c(-98,-81), ylim=c(41.5,50))
  points(unlist(long.keep),unlist(lat.keep), pch=21,
         cex=1.1, bg=colors[data_binned],lwd=.2)
  plotInset(-90,47,-82.5,50,
            expr={
              keep.col<-unique(data_binned)
              keep.col<-keep.col[!is.na(keep.col)]
              keep.col<-sort(keep.col)
              is.na.vec <- rep(NA,10)
              is.na.vec[keep.col]<-colors[keep.col]
              
              hist(data_binned,col=is.na.vec,xaxt="n",xlab=NA,
                   ylab=NA,main=NA,cex.lab=.5,cex.axis=.5,
                   xlim=c(0,length(breaks)),ylim=c(0,20),breaks=seq(0,12,1))
              
              axis(side=1,breaks,at = seq(0,11,1),cex.axis = .5,las=2,line=0)
              mtext(side = 1, "Biomass (Mg/ha)", line = 1.5,cex=.5)
              mtext(side = 2, "Frequency", line = 1.7,cex=.5)
            })
  title(paste("Biomass @",name.vec[r]))
}
dev.off()

#####
##### Second Derivative #####
#####

calc.second.deriv <- function(biomassCI,h,second.deriv){
  for(i in 1:length(biomassCI)){
    if(length(biomassCI[[i]])>10){
      biomassCI[[i]] <- biomassCI[[i]]
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
for(i in c(1)){
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

pdf('stability.group.ts.100.pdf')
par(mfrow=c(2,3))
for(count in 1:6){
  plot(seq(100,10000,100),seq(100,10000,100),
       cex=.1,ylim=c(0,150),xlim=c(10000,-10),ylab="Biomass (Mg / Ha)",
       xlab="Years Before Present",main=paste('Stability Bin',count))
  for(r in names(second.deriv.unlist[data_binned==count])){
    biomassCI.use <- biomassCI[[which(unique(x.meta$site.id)==as.numeric(r))]]
    print(x.meta[x.meta$site.id==r,'site.name'])[[1]]
    points(seq(100,10000,100),biomassCI.use[2,],
           main=(x.meta[x.meta$site.id==r,'site.name'])[[1]],
           ylim=c(0,150), pch = 21,bg=colors[count])
  }
}
dev.off()

#####
##### Principle Component Analysis #####
#####

prop.ten.count = prop.table(ten.count,margin = 1)

#whole sample #oak,pine,prairie
pca=princomp(prop.ten.count) 
summary(pca) 
biplot(pca)

#pca by stability group #number is time step
pdf('stability.group.pca.2000.pdf')
for(count in 1:6){
  pca.new <- numeric(ncol(prop.ten.count))
  for(b in as.numeric(names(second.deriv.unlist[data_binned==count]))){
    pca.use <- prop.ten.count[x.meta$site.id==b,]
    #print(pca.use)
    pca.new <- rbind(pca.use,pca.new)
  }
pca=princomp(pca.new) 
#summary(pca)
biplot(pca,main=paste('Stability Bin',count))
}
dev.off()

save.cats <- list()
for(i in 1:nrow(prop.ten.count)){
  #save.cats[[i]] <- length(which(prop.ten.count[i,] > .1))
  save.cats[[i]] <- diversity(prop.ten.count[i,])
}

x.meta.new <- cbind(x.meta, unlist(save.cats))
x.meta.new <- cbind(x.meta.new, numeric(nrow(x.meta)))
colnames(x.meta.new) <- c(colnames(x.meta),c('divIndex','estBiomass'))

names(biomassCI) <- unique(x.meta$site.id)[1:182]

pdf('diversity.biomass.pdf')
plot(save.cats[which(x.meta$site.id==site.num)],
     biomassCI[[182]][2,x.meta[x.meta$site.id==site.num,'age_bacon']/100],
     pch = 19,ylab = 'Biomass',xlab = 'Diversity Index',
     ylim = c(0,150),xlim=c(min(unlist(save.cats)),max(unlist(save.cats))))

for(i in 1:182){
  if(length(biomassCI[[i]])>1){
    site.num <- names(biomassCI[i])
    points(save.cats[which(x.meta$site.id==site.num)],
         biomassCI[[i]][2,x.meta[x.meta$site.id==site.num,'age_bacon']/100],
         pch = 19)
    x.meta.new[which(x.meta$site.id==site.num),'estBiomass'] <- biomassCI[[i]][2,x.meta[x.meta$site.id==site.num,'age_bacon']/100]
    
 }
}
dev.off()


map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(x.meta[which(x.meta.new$divIndex>1.75&x.meta.new$estBiomass>100),'long'],
       x.meta[which(x.meta.new$divIndex>1.75&x.meta.new$estBiomass>100),'lat'],
       pch=19,
       cex=1.3,lwd=.2)
hist(x.meta[which(x.meta.new$divIndex>1.75&x.meta.new$estBiomass>100),'age_bacon'])

map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(x.meta[which(x.meta.new$divIndex>1.75&x.meta.new$estBiomass<40),'long'],
       x.meta[which(x.meta.new$divIndex>1.75&x.meta.new$estBiomass<40),'lat'],
       pch=19,
       cex=1.3,lwd=.2)
hist(x.meta[which(x.meta.new$divIndex>1.75&x.meta.new$estBiomass<40),'age_bacon'])

map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(x.meta[which(x.meta.new$divIndex<1.25),'long'],
       x.meta[which(x.meta.new$divIndex<1.25),'lat'],
       pch=19,
       cex=1.3,lwd=.2)


