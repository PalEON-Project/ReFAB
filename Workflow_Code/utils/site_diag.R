site_diag <- function(bMax = 150, locn, x.meta=NULL, minAge = 0, maxAge = 10000, 
                      ageInterval = 100, path_to_samps = '~/Downloads/samps2zip/',
                      path_to_Info = NULL, control.pts = NULL,
                      ten.count = NULL, dataID = NULL,out.dir = getwd()){

names_use <- c('Pinus','Prairie','Quercus','Betula',
               'Other Herb.','Cyperace','Alnus','Ostrycar',
               'Ulmus','Tsugax','Picea','Acer','Fraxinus',
               'Populus','Cupressa','Other Trees','Larix',
               'Fagus','Carya','Abies','Tilia','Juglanus')
  
if(is.null(x.meta)){
  print('Error: Need meta data for prediction dataset. Look in create_prediction_datasets.R for file name.')
  stop()
}
if(is.null(ten.count)){
    print('Error: Need pollen count data for prediction dataset. Look in create_prediction_datasets.R for file name.')
    stop()
}
if(is.null(dataID)){
    print('Error: Need run info for prediction dataset. Look in create_prediction_datasets.R for file name.')
    stop()
}
  
locnClean <- gsub(' ', '-', locn)
site_number = unique(x.meta[x.meta$site.name == locn,1])

x.meta.use <- x.meta[x.meta$site.name == locn,]

source(file.path('Workflow_Code','utils','test_site.R'))
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

for(b in 1:50){
  ID <- dataID[dataID$name==as.character(locn),'ID'][b]
  file_name <- paste0(path_to_samps,'samplesList_workInfo_',ID,'_',locnClean,'_Beta_',b,'.Rdata') #Sigma0.12Group
  if(!file.exists(file_name)) next()
  load(file_name)
  samples.keep <- rbind(samples.keep, samplesList[sample(x = 1:nrow(samplesList),size = 100),])
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
pdf(file.path(out.dir,paste0('SiteDiagnositcs-Age',locnClean,'.pdf')))
#quartz()
site_number = unique(x.meta[x.meta$site.name == locn,1])
keep.dataset.id <- unique(x.meta[x.meta$site.id==site_number,4])

breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)))

bio.quants <- apply(samplesList[,1:(maxAge/100)],2,quantile,c(0.025,0.5,0.975))

data_binned <-  cut(rev(bio.quants[2,]), c(breaks), include.lowest = FALSE, labels = FALSE)

load('dist2bigwoods.Rdata') #from no analog assessment
d2b_site <- dist2bigwoods[which(x.meta$site.name == locn)]
breaks_d2b <- seq(0,1,.05)
red2white <- colorRampPalette(colors=c('black','lightyellow','white'))
colors_d2b <- red2white(length(breaks_d2b))#c('darkred','red','orange','yellow','white')#(heat.colors(length(breaks)))

data_binned_d2b <- cut(d2b_site, c(breaks_d2b), include.lowest = TRUE, labels = FALSE)

prop.use <- prop.table(as.matrix(Y),margin=1)

prop.max <- apply(prop.use[!is.na(age_index),],2,max)+.05

prop.vec <- round(140*prop.max/sum(prop.max))
prop.vec[prop.vec==0] <- 1

prop.list <- list()
for(i in 1:22){
  prop.list[[i]] <- rep(i+1,prop.vec[i])
}

if(length(unlist(prop.list))!=140){
  print(length(unlist(prop.list)))
}

fig.mat <- matrix(1,length(unlist(prop.list))+40,1)
fig.mat[1:40,]<-1
fig.mat[41:(length(unlist(prop.list))+40),]<-unlist(prop.list)
#control.pts<-read.csv(text=getURL('https://raw.githubusercontent.com/PalEON-Project/stepps-baconizing/master/data/bacon_age_markers_v1.csv'))

def.par <- par(no.readonly = TRUE)
layout(fig.mat)
par(mar = c(0, 4, 0, 3), oma = c(5, 0, 4, 0),tcl=.5)

plot(bio.quants[2,], xlim = c(maxAge/100,0), ylim = c(0,bMax), xaxt='n',
     xlab = 'Years BP', ylab = 'Biomass Mg/ha', col = 'white')
axis(side = 3, at = rev(seq(0,maxAge/100, round((maxAge/100)/10))), labels = rev(seq(0,maxAge,round(maxAge/10))))
#abline(v=age_index,col=colors_d2b[data_binned_d2b],lwd=2)
ciEnvelope(x=1:(maxAge/100), ylo = bio.quants[1,],yhi = bio.quants[3,],col = 'grey')
points(bio.quants[2,],cex=1.1,typ='l',col='darkgray')#pch=16,col = rev(colors[data_binned]))
rug(age_index,lwd=2)
rug(control.pts[which(control.pts[,2]%in%keep.dataset.id),]$geo_age/100,lwd=3,col="red")
if(!is.null(path_to_Info)){
  for(b in 1:50) {
    points(age_index, out.list[[b]], cex = .5,col='black')
  }
}
#points(0,unique(x.meta[x.meta$site.name == locn,'SettleBiomass']),pch=19,col='purple',cex=2)
legend('topleft','Mx.Lik.',pch=1)
for(p in 1:22){#rev(match(names(sort(colMeans(prop.use))),colnames(prop.use)))
  prop.plot<- cbind(age_index*100,as.matrix(prop.use[,p]))      	
  prop.plot<-prop.plot[order(prop.plot[,1]),]
  plot(x=prop.plot[,1],y=prop.plot[,2],type="l",xlim=c(maxAge,-10),
       ylim=c(0,max(prop.use[!is.na(age_index),p])),ylab=NA,yaxt='n', xaxt='n')
  max.me <- max(prop.use[!is.na(age_index),p])
  if(max.me>.05){
    axis(2,at = c(max.me,max.me/2),
         labels=c(signif(max.me,digits=1), names_use[p]),
         las=1,cex.axis=.65)
  }else{
    axis(2,at = c(max.me/2),
         labels=c(names_use[p]),
         las=1,cex.axis=.4)
  }
  ciEnvelope(prop.plot[,1],rep(0,nrow(prop.use)),prop.plot[,2],col="darkblue")
  #legend('topright',colnames(prop.use)[p])
  #legend('topleft',paste(signif(mean(prop.use[,p]),digits=2)),bty="n")
} 
par(def.par)

#####
##### Map #####
#####

maps:::map('state', xlim=c(-98,-81), ylim=c(41.5,50))
if(is.null(x.meta[x.meta$site.name == locn,'long.x'])){
  points(unique(x.meta[x.meta$site.name == locn,'long']),
         unique(x.meta[x.meta$site.name == locn,'lat']),
         pch=19,cex=1.5)
  
}
points(unique(x.meta[x.meta$site.name == locn,'long.x']),
       unique(x.meta[x.meta$site.name == locn,'lat.x']),
       pch=19,cex=1.5)
title(locn)


#####
##### MCMC #####
#####

par(def.par)

par(mfrow=c(4,4))
for(i in 1:(maxAge/100)){
  plot(
    samplesList[, i],
    ylab = 'Biomass Estimate',
    xlab = 'MCMC iteration',
    main = paste(i * 100, 'Age BP MCMC'),
    typ = 'l',
    ylim = c(0, bMax)
  )
  if(any(i==age_index)){
    if(!is.null(path_to_Info)){
    for(b in 1:50){
      abline(h=out.list[[b]][which(i==age_index)],
             col=adjustcolor('purple',alpha.f = .5),
             lwd=1)
    }
      get.max <- list()
      for(b in 1:50){
        get.max[[b]] <- exp(out.keep[[b]][, which(i==age_index)] - max(out.keep[[b]][, which(i==age_index)])) / -sum(out.keep[[b]][, which(i==age_index)])
      }
      max.me <- max(unlist(get.max))
      
          for(b in 1) {
            plot(
              seq(5, bMax - 5, by = 2),
              exp(out.keep[[b]][, which(i==age_index)] - max(out.keep[[b]][, which(i==age_index)])) / -sum(out.keep[[b]][, which(i==age_index)]),
              typ = 'l',
              main = paste(i * 100, 'Age BP Log Lik.'),
              xlab = 'Biomass',
              ylab = 'Log Lik.',
              ylim = c(0,max.me)
            )
          }
            for(b in c(10,20,30,40,50)){
              points(seq(5, bMax - 5, by = 2),
                     exp(out.keep[[b]][, which(i==age_index)] - max(out.keep[[b]][, which(i==age_index)])) / -sum(out.keep[[b]][, which(i==age_index)]),
                     typ = 'l',col = fields::tim.colors(50)[b])
            }
      rug(samplesList[,i])
          }
    #abline(h=seq(5, bMax-5, by = 2)[apply(out,2,which.max)][which(i==age_index)],col='purple',lwd=3)
  }else{
    plot.new()
    legend('center',paste('No data during',i*100,'\nyears BP time interval'),cex = .75)
    }
}
#plot(samplesList[,grep('sigma',colnames(samplesList))],ylab = 'Sigma Estimate', xlab = 'MCMC iteration', main = 'Sigma', typ='l')




dev.off()
}
