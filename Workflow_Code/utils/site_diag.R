site_diag <- function(bMax = 150, locn, x.meta, minAge = 0, maxAge = 10000, 
                      ageInterval = 100, path_to_samps = '~/Downloads/samps2zip/',
                      path_to_Info = '~/Downloads/work_3/', control.pts,
                      ten.count){
  
locnClean <- gsub(' ', '-', locn)
site_number = unique(x.meta[x.meta$site.name == locn,1])

x.meta.use <- x.meta[x.meta$site.name == locn,]

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

Y <- as.matrix(Y2[ , -c(1,2)])
age_index <- Y2[,1]
samples.keep <- numeric(300)

out.list <- out.keep <- list()

for(b in 1:20){
  ID <- dataID[dataID$name==as.character(locn),'ID'][b]
  file_name <- paste0(path_to_samps,'samplesList_workInfo_',ID,'_',locnClean,'_Beta_',b,'.Rdata') #Sigma0.12Group
  if(!file.exists(file_name)) next()
  load(file_name)
  samples.keep <- rbind(samples.keep, samplesList)
  file_name1 <- paste0(path_to_Info,'workInfo_',ID,'_',locnClean,'_Beta_',b,'.Rdata')
  if(!file.exists(file_name1)){
    out <- NA
  }else{
    load(file = file_name1)
    out.keep[[b]] <- out
    out.list[[b]] <- seq(5, bMax-5, by = 2)[apply(out,2,which.max)]
  }
  
}

samplesList <- samples.keep

#### Plotting One Site
pdf(paste0('SiteDiagnositcs-Age',locnClean,'.pdf'))
#quartz()
site_number = unique(x.meta[x.meta$site.name == locn,1])
keep.dataset.id <- unique(x.meta[x.meta$site.id==site_number,4])

breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)))

bio.quants <- apply(samplesList[,1:(maxAge/100)],2,quantile,c(0.025,0.5,0.975))

data_binned <-  cut(rev(bio.quants[2,]), c(breaks), include.lowest = FALSE, labels = FALSE)

fig.mat <- matrix(1,28,1)
fig.mat[1:6,]<-1
fig.mat[7:28,]<-seq(2,23,1)
#control.pts<-read.csv(text=getURL('https://raw.githubusercontent.com/PalEON-Project/stepps-baconizing/master/data/bacon_age_markers_v1.csv'))

def.par <- par(no.readonly = TRUE)
layout(fig.mat)
par(mar = c(0, 4, 0, 3), oma = c(5, 0, 4, 0),tcl=.5)

plot(bio.quants[2,], xlim = c(maxAge/100,0), ylim = c(0,bMax), xaxt='n',
     xlab = 'Years BP', ylab = 'Biomass Mg/ha', col = 'white')
axis(side = 3, at = rev(seq(0,maxAge/100, round((maxAge/100)/6))), labels = rev(seq(0,maxAge,round(maxAge/6))))
ciEnvelope(x=1:(maxAge/100), ylo = bio.quants[1,],yhi = bio.quants[3,],col = 'grey')
points(bio.quants[2,],cex=1.1,pch=16,col = rev(colors[data_binned]))
rug(age_index,lwd=2)
rug(control.pts[which(control.pts[,2]%in%keep.dataset.id),]$geo_age/100,lwd=3,col="red")
for(b in 1:20) {
  points(age_index, out.list[[b]], cex = .5)
}
#points(0,unique(x.meta[x.meta$site.name == locn,'SettleBiomass']),pch=19,col='purple',cex=2)
legend('topleft','Mx.Lik.',pch=1)

ten_count_use = ten.count[which(x.meta$site.name == locn), ]

Y = as.matrix(ten_count_use)
prop.use <- prop.table(as.matrix(Y),margin=1)    

for(p in rev(match(names(sort(colMeans(prop.use))),colnames(prop.use)))){
  prop.plot<- cbind(sample_ages,as.matrix(prop.use[,p]))      	
  prop.plot<-prop.plot[order(prop.plot[,1]),]
  plot(x=prop.plot[,1],y=prop.plot[,2],type="l",xlim=c(maxAge,-10),
       ylim=c(0,max(prop.use[,p])),ylab=NA,yaxt='n', xaxt='n')
  #axis(2,at=signif(max(prop.use)/2),labels=signif(max(prop.use[,p])/2,digits=1))
  ciEnvelope(prop.plot[,1],rep(0,nrow(prop.use)),prop.plot[,2],col="darkblue")
  legend('topleft',colnames(prop.use)[p])
  #legend('topleft',paste(signif(mean(prop.use[,p]),digits=2)),bty="n")
} 
par(def.par)

par(mfrow=c(3,3))
for(i in 1:(maxAge/100)){
  plot(samplesList[,i],ylab = 'Biomass Estimate', xlab = 'MCMC iteration', main = i, typ='l',ylim=c(0,150))
  if(any(i==age_index)){
    for(b in 1:20){
      abline(h=out.list[[b]][which(i==age_index)],col='purple',lwd=1)
    }
    #abline(h=seq(5, bMax-5, by = 2)[apply(out,2,which.max)][which(i==age_index)],col='purple',lwd=3)
  }
}
#plot(samplesList[,grep('sigma',colnames(samplesList))],ylab = 'Sigma Estimate', xlab = 'MCMC iteration', main = 'Sigma', typ='l')

map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(unique(x.meta[x.meta$site.name == locn,'long']),
       unique(x.meta[x.meta$site.name == locn,'lat']),
       pch=19,cex=1.5)
title(locn)



for(i in 1:ncol(out)){
 for(b in c(1,5,10,15)){
    plot(seq(5, bMax-5, by = 2),
          exp(out.keep[[b]][,i]-max(out.keep[[b]][,i]))/-sum(out.keep[[b]][,i]),
          typ='l',main=paste('beta=',b,'age',age_index[i]))
  }
}

dev.off()
}
