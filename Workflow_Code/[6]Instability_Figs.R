
#### Do load data portion of average biomass figure


library(analogue)
library(viridis)

dist.keep <- dist.all <- list()
for(i in 1:80){
  
  if(!is.null(Y.keep[[i]])&&!is.null(biomassCI[[i]])){
    
    Y.keep[[i]] <- Y.keep[[i]][order(Y.keep[[i]][,1]),]
    
    dist.keep[[i]] <- numeric(nrow(Y.keep[[i]]))
    for(r in 2:nrow(Y.keep[[i]])){
      dist.keep[[i]][r] <- as.numeric(distance(Y.keep[[i]][r-1,3:24],
                                               Y.keep[[i]][r,3:24],
                                               method = 'SQchord'))
    } 
    
    time.diff <- diff(Y.keep[[i]][,1])
    dist.all[[i]] <- data.frame(#dist.keep[[i]][2:nrow(Y.keep[[i]])]/time.diff,
      which_max = seq(5, bMax-5, by = 2)[apply(out.list[[i]][[1]],2,which.max)],
      max = apply(out.list[[i]][[1]],2,max),
      biomass = (biomassCI[[i]][2,Y.keep[[i]][,1]]),
      #time.diff,rep(i,length(time.diff)),
      age_index = Y.keep[[i]][,1],
      #Y.keep[[i]][,3:24])
      prop.table(as.matrix(Y.keep[[i]][,3:24]),1))
  }else{dist.all[[i]]<-NA}
  
}

distA <- do.call(rbind,dist.all)
distA <- as.data.frame(distA[-which(is.na(distA[,3])),])

pdf('[15]max_lik_by_age.pdf')
par(mfrow=c(1,1))
plot(distA$age_index,distA$max,pch=19,
     cex=.5,ylab='Maximum Likelihood Value',
     xlab = 'Age of Sample (YBP1950)',xaxt='n', xlim = c(100,0))
axis(side = 1, at = seq(0,100,10), labels = seq(0,10000,1000))
dev.off()

comp_mat <- b_mat <- matrix(NA,nrow(distA),length( distance(distA[i,5:ncol(distA)],distA[,5:ncol(distA)],method='SQchord') ))

for(i in 1:nrow(distA)){
  comp_mat[i,] <- distance(distA[i,5:ncol(distA)],distA[,5:ncol(distA)],method='SQchord') #distance between pollen samples
  b_mat[i,] <- distA[i,"which_max"] - distA[,"which_max"] #difference in max lik est of biomass
  print(i)
}

add_legend <- function(){
  xm <- get('xm', envir = parent.frame(1))
  ym <- get('ym', envir = parent.frame(1))
  z  <- get('dens', envir = parent.frame(1))
  colramp <- get('colramp', parent.frame(1))
  fields::image.plot(xm,ym,z, col = colramp(256), legend.only = T, add =F)
}

pdf('[12-B]compvbiom_pred_smooth.pdf')
par(mar = c(5,4,4,5) + .1)
smoothScatter(c(comp_mat), abs(c(b_mat)),
              colramp = plasma, 
              postPlotHook = add_legend,
              xlab='Distance between compositions',
              ylab='Difference between biomass Max. Lik.',
              nrpoints = 0,main=NA,
              ylim = c(0,300),xlim=c(0,2))
dev.off()
##### Calibration dataset

load('threethirds_v3.0.Rdata')
load(file='out_calc.Rdata')

which_max <- apply(out_calc[[20]],1,which.max)
comp_mat_c <- b_mat_c <- matrix(NA,nrow(Y),nrow(Y))

P <- as.data.frame(prop.table(as.matrix(Y),1))

for(i in 1:nrow(P)){
  comp_mat_c[i,] <- analogue::distance(P[i,1:ncol(P)],P[,1:ncol(P)],method='SQchord') #distance between pollen samples
  b_mat_c[i,] <- which_max[i] - which_max #difference in max lik est of biomass
  print(i)
}

pdf('[12-A]compvbiom_calib_smooth.pdf')
par(mar = c(5,4,4,5) + .1)
smoothScatter(c(comp_mat_c), abs(c(b_mat_c)),
              colramp = plasma, 
              postPlotHook = add_legend,
              xlab='Distance between compositions',
              ylab='Difference between biomass max. lik.',
              nrpoints = 0, main = NA,
              ylim = c(0,300),xlim=c(0,2))
dev.off()


################### Map showing sites with different problems

####### Calibration
calibration_underpredicted <- c(3,  10,  11,  13,  18,  29,  32,  54,  58,  59,  71,  73,  82,  83,
                                84, 105, 109)#c(2,10,11,13,14,18,29,31,32,54,59,82,83)
calibration_bimodal <- c(1,9,10,13,37,50,101,105,115,149) #counting likelihoods that are different between cv and 2/3s

breaks <- c(0,10,25,50,100,150,200,250,300)
load('threethirds_v3.0.Rdata')
color_num <- as.numeric(cut(biomass, breaks=breaks,labels=1:(length(breaks)-1)))
colors <- rev(terrain.colors(length(breaks)-1))

load("~/ReFAB/2019-05-22two.thirds.cast.x.Rdata")
load('cast.x.Rdata')

pdf('[13-A]map_unstable_calibration.pdf')
maps::map('state', xlim=c(-98,-81), ylim=c(40,50))
#mtext(side=3,'Calibration Dataset',cex=2)

points(cast.x$long,
       cast.x$lat,
       pch=21,cex=1,bg=colors[color_num])

legend(
  'topright',
  c('0-10', '10-25', '25-50', '50-100', '100-150', '150-200','200-250','250-300',
    'Under Predicted','Bimodal Likelihood'),
  pch = c(rep(21,8),6,0),
  cex = .75,
  pt.bg = c(colors),
  col = c(rep('black',8),'green','blue'),
  title = 'Settlement Biomass'
)

points(ag.two.thirds.cast.x$long[calibration_underpredicted],
       ag.two.thirds.cast.x$lat[calibration_underpredicted],
       pch=6,col='green',cex=3)
points(ag.two.thirds.cast.x$long[calibration_bimodal],
       ag.two.thirds.cast.x$lat[calibration_bimodal],
       pch=0,col='blue',cex=3)

dev.off()


biom_8000 <- lapply(biomassCI,FUN=function(x)x[2,80])

breaks <- c(seq(0,300,25))
colors <- rev(terrain.colors(length(breaks)-1))

pdf('[13-B]map_prediction_unstable.pdf')
maps::map('state', xlim=c(-98,-81), ylim=c(40,50))
#mtext(side=3,'Prediction Dataset 8000 Years BP',cex=2)

for(i in 1:80){
  locn <- as.character(unique(dataID$name)[i])
  if(locn=='Mud Lake' | locn == 'Lily Lake' | locn == 'Chatsworth Bog') next
  print(paste(locn,biom_8000[i]))
  color_num <- as.numeric(cut(biom_8000[[i]], breaks=breaks,labels=1:(length(breaks)-1)))
  points(x.meta[x.meta$site.name == locn,'long'][1],
         x.meta[x.meta$site.name == locn,'lat'][1],
         pch=21,cex=1.5,bg=colors[color_num])
}

legend(
  'topright',
  c('0-10', '10-25', '25-50', '50-100', '100-150', '150-200','200-250','250-300','Bimodal'),
  pch = c(rep(21,8),0),
  cex = .75,
  pt.bg = c(colors),
  col = c(rep('black',8),'blue','green'),
  title = 'Predicted Biomass'
)

bimodal_change <- c('Radtke Lake', 'Devils Lake', 'Lake Mendota',
                    'Wells Mastodon Site', 'Volo Bog', 'Gass Lake',
                    'Fox Lake', 'Reidel Lake', 'Emrick Lake',
                    'Lima Bog', 'Kimble Pond', 'Disterhaft Farm Bog',
                    'Kellners Lake', 'Clear Lake')

for(i in 1:length(bimodal_change)){
  locn <- bimodal_change[i]
  print(locn)
  site_number = unique(x.meta[x.meta$site.name == locn,1])
  points(x.meta[x.meta$site.name == locn,'long'][1],
         x.meta[x.meta$site.name == locn,'lat'][1],
         pch=0,col='blue',
         cex = 4)
}
dev.off()









bin_c <- hexbin(c(comp_mat_c), abs(c(b_mat_c)))
plot(bin_c) 

pdf('compvbiom_calib.pdf')
# plot  - look at str(p)
p_c <- plot(bin_c, xlab='Distance between compositions',
            ylab='Difference between biomass Max. Lik.',
            main = 'Calibration Dataset')

#takes too long to draw
#plot(c(comp_mat), abs(c(b_mat)), col = rgb(0, 0, 0, .1), pch = 16, 
 #    cex = .5, main = 'transparency') 


plot(distA$biomass,distA$which_max)

distB <- distA[which(distA$biomass<50&distA$which_max>50),]

breaks <- c(seq(-126,0,25))
color_num <- as.numeric(cut(distA$max, breaks=breaks,labels=1:(length(breaks)-1)))
colors <- plasma(6)

plot(distB$age_index,distB$max)
pc <- princomp(distA[,5:ncol(distA)])

pdf('biplot_prediction.pdf')
biplot(pc)

raw_pred <- pc$scores[,1:2]
plot(raw_pred[,1],raw_pred[,2],bg=colors[color_num],pch=21,cex=2)
dev.off()

comp_mat <- b_mat <- matrix(NA,nrow(distA),nrow(distA))


lt <- sample(size=10000,x=1:length(c(comp_mat)))

png('compvbiom2.png')
plot(c(comp_mat)[lt],abs(c(b_mat)[lt]))
dev.off()

#plot(c(comp_mat), abs(c(b_mat)), col = rgb(0, 0, 0, .1), pch = 16, 
#    cex = .5, main = 'transparency') 
library(hexbin)
bin <- hexbin(c(comp_mat), abs(c(b_mat)))
plot(bin) 

pdf('compvbiom_pred.pdf')
# plot  - look at str(p)
p <- plot(bin, xlab='Distance between compositions',ylab='Difference between biomass Max. Lik.',
          main='Prediction Dataset')

# push plot viewport
pushHexport(p$plot.vp)

hem <- which(distA[,14]>.1)[1:50]
prar <- which(distA[,6]>.5)[1:50]

grid.points(comp_mat[hem,prar], abs(b_mat[hem,prar]), pch=1,
            gp=gpar(col=adjustcolor("blue",alpha.f = .25),cex=.5))
dev.off()

##### Calibration dataset

load('threethirds_v2.0.Rdata')
load(file='out_calc.Rdata')

which_max <- apply(out_calc[[20]],1,which.max)
comp_mat_c <- b_mat_c <- matrix(NA,nrow(Y),nrow(Y))

P <- as.data.frame(prop.table(as.matrix(Y),1))

for(i in 1:nrow(P)){
  comp_mat_c[i,] <- analogue::distance(P[i,1:ncol(P)],P[,1:ncol(P)],method='SQchord') #distance between pollen samples
  b_mat_c[i,] <- which_max[i] - which_max #difference in max lik est of biomass
  print(i)
}

bin_c <- hexbin(c(comp_mat_c), abs(c(b_mat_c)))
plot(bin_c) 

pdf('compvbiom_calib.pdf')
# plot  - look at str(p)
p_c <- plot(bin_c, xlab='Distance between compositions',
            ylab='Difference between biomass Max. Lik.',
            main = 'Calibration Dataset')

# push plot viewport
pushHexport(p_c$plot.vp)

hem <- which(P[,10]>.1)
prar <- which(P[,2]>.1)

#grid.points(comp_mat_c[hem,prar], abs(b_mat_c[hem,prar]), pch=1,
#            gp=gpar(col=adjustcolor("blue",alpha.f = .25),cex=.5))
dev.off()


### All Sites
library(lattice)
pdf('compvbiom_allsites.pdf')
for(k in 1:length(unique(dataID$name))){
  
  if(k==34) next
  
  if(!is.na(dist.all[[k]])){
    distLM <- dist.all[[k]]
    
    comp_matLM <- b_matLM <- matrix(NA,nrow(distLM),nrow(distLM))
    
    for(i in 1:nrow(distLM)){
      comp_matLM[i,] <- distance(distLM[i,5:ncol(distLM)],distLM[,5:ncol(distLM)],method='SQchord') #distance between pollen samples
      b_matLM[i,] <- distLM[i,"which_max"] - distLM[,"which_max"] #difference in max lik est of biomass
      
    }
    print(k)
    
    p_LM <-hexbinplot(abs(c(b_matLM))~c(comp_matLM), xlab='Distance between compositions',
                      ylab='Difference between biomass Max. Lik.',
                      main = as.character(unique(dataID$name)[k]),
                      ylim = c(-5,228),xlim=c(-.025,1),aspect = 1)
    
    # push plot viewport
    print(p_LM)
  }
}
dev.off()

library(gridExtra)    


lat_lon_mat <- as.matrix(ag.two.thirds.cast.x[,c('lat','long')])
rownames(lat_lon_mat) <- 1:nrow(lat_lon_mat)
bimodal_min <- numeric(length(calibration_bimodal))
for(i in 1:length(calibration_bimodal)){
  bimodal_min[i] <- as.numeric(names(sort(as.matrix(dist(lat_lon_mat))[calibration_bimodal[i],])[2]))
}

bimodal_min <- c(8,36,11,102,2)

box_df <- data.frame(prop.table(as.matrix(rbind(Y[calibration_bimodal,],Y[bimodal_min,])),1),
                     sites =c(rep('bimodal',length(calibration_bimodal)),rep('near',length(bimodal_min))))

box_melt <- melt(data = box_df, id.vars = 'sites')

pdf('bimodal_v_near_boxplot.pdf')
ggplot(data = box_melt, aes(x = variable, y = value)) + geom_boxplot(aes(fill=sites))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y = 'Pollen Proportion',
       x = 'Taxa',
       title = 'Calibration Dataset Bimodal v. Near Bimodal Sites')
dev.off()

pdf('bimodal_v_near_scatter.pdf')
plot(biomass,colMeans(samps2_3))
points(biomass[calibration_bimodal],colMeans(samps2_3)[calibration_bimodal],col='blue',pch=19)
points(biomass[bimodal_min],colMeans(samps2_3)[bimodal_min],col='green',pch=19)
abline(a=0,b=1)
legend('bottomright',c('bimodal','near'),pch=19,col=c('blue','green'))

boxplot(biomass[calibration_bimodal],biomass[bimodal_min],col=c('blue','green'),main='Settlement Biomass')
boxplot(cbind(samps2_3[4000:5000,calibration_bimodal],samps2_3[4000:5000,bimodal_min]),col=c(rep('blue',length(calibration_bimodal)),rep('green',length(bimodal_min))),main='Predicted Biomass')
points(c(biomass[calibration_bimodal],biomass[bimodal_min]),col='red',pch=19)
legend('topright',c('bimodal','near','settlement'),pch=19,col=c('blue','green','red'))


dev.off()

load('twothirds_v2.0.Rdata')

rownames(Y)[calibration_bimodal] <- paste0('AAAA',calibration_bimodal)
pc_calib <- princomp(prop.table(as.matrix(Y[,]),1))

pdf('calibration.biplot.pdf')
biplot(pc_calib)
raw <- pc_calib$scores[,1:2]
#points(raw[calibration_underpredicted,1], raw[calibration_underpredicted,2],
#       col='blue', pch=20,
#     cex=2)
points(raw[calibration_bimodal,1], raw[calibration_bimodal,2],
       col='green', pch=20,
       cex=4)
legend('topright',c('underpredicted','bimodal'),pch=20,col=c('blue','green'))

breaks <- c(seq(0,250,25))
color_num <- as.numeric(cut(biomass, breaks=breaks,labels=1:(length(breaks)-1)))
colors <- terrain.colors(length(breaks)-1)

plot(raw[,1], raw[,2],
     bg=colors[color_num], pch=21,
     cex=2,main='biplot colored by biomass')
points(raw[calibration_bimodal,1], raw[calibration_bimodal,2],
       col='green', pch=20,
       cex=4)

legend('topright',legend = breaks[1:(length(breaks)-1)],
       col=colors,pch=19)

dev.off()



if(FALSE){ #not sure if this is useful
  p_c <- plot(bin_c, xlab='Distance between compositions',
              ylab='Difference between biomass Max. Lik.',
              main = 'Calibration Dataset')
  pushHexport(p_c$plot.vp)
  grid.points(c(comp_mat_c[calibration_underpredicted,calibration_underpredicted]),
              c(abs(b_mat_c[calibration_underpredicted,calibration_underpredicted])),
              pch=1,
              gp = gpar(col=adjustcolor(colors[1],alpha.f = .25),
                        cex=.5))
  grid.points(c(comp_mat_c[calibration_bimodal,calibration_bimodal]),
              c(abs(b_mat_c[calibration_bimodal,calibration_bimodal])),
              pch=8,
              gp = gpar(col=adjustcolor(colors[2],alpha.f = .25),
                        cex=.5))
  grid.points(c(comp_mat_c[which(out_sd>quantile(out_sd,.95)),which(out_sd>quantile(out_sd,.95))]),
              c(abs(b_mat_c[which(out_sd>quantile(out_sd,.95)),which(out_sd>quantile(out_sd,.95))])),
              pch=5,
              gp = gpar(col=adjustcolor(colors[3],alpha.f = .25),
                        cex=.5))
}


##### Prediction
pred_which_max <- pred_which_sd <- look <- list()
for(i in 1:80){
  pred_which_max[[i]] <- list()
  if(!is.null(out.list[[i]])){
    for(b in 1:20){
      pred_which_max[[i]][[b]] <- apply(out.list[[i]][[b]],2,which.max)
    }
    look[[i]] <- apply(do.call(rbind,pred_which_max[[i]]),2,sd)
    pred_which_sd[[i]] <- mean(look[[i]]) 
    look[[i]] <- cbind(look[[i]],Y.keep[[i]])
  }else{
    pred_which_sd[[i]] <- NA
  }
  
}

look_mat <- do.call(rbind,look)

pc <- princomp(x = prop.table(as.matrix(distA[,5:26]),1))

pdf('pred_biplot.pdf')
biplot(pc)

raw <- pc$scores[,1:2]
breaks <- c(seq(0,250,50))
color_num <- as.numeric(cut(distA$biomass, breaks=breaks,labels=1:(length(breaks)-1)))
colors <- terrain.colors(length(breaks)-1)

plot(raw[,1],raw[,2],bg=colors[color_num],pch=21,cex=2,main='colored by predicted biomass')
legend('topright',legend = breaks[1:(length(breaks)-1)],col=colors,pch=19)
dev.off()



pred_sd <- unlist(pred_which_sd)
sd_pred_pts <- which(pred_sd >= quantile(pred_sd,.75,na.rm=T))

unique(dataID$name)[order(pred_sd)]

breaks <- seq(0,25,5)
color_num <- as.numeric(cut(pred_sd, breaks=breaks,labels=1:(length(breaks)-1)))

linear_change <- c(8,9,13,57,65,77,79,70)
sparse_data <- c(3,39)
bimodal_change <- sd_pred_pts[-which(sd_pred_pts%in%c(linear_change,sparse_data,73))]

dev.off()

rownames(P)[calibration_underpredicted] <- paste0('U',calibration_underpredicted)
rownames(P)[calibration_bimodal] <- paste0('b',calibration_bimodal)

pc <- princomp(P)


pdf('predict.biplot.pdf')

biplot(pc)
raw_pred <- pc$scores[,1:2]
points(raw_pred[calibration_underpredicted,1], raw_pred[calibration_underpredicted,2],
       col='blue', pch=20,
       cex=2)
points(raw[calibration_bimodal,1], raw[calibration_bimodal,2],
       col='green', pch=20,
       cex=2)


plot(raw[,1], raw[,2],
     bg=colors[color_num], pch=21,
     cex=2,main='biplot colored by SD')


dev.off()


lapply(dist.all,function(x){
  c(sum(x[,1]),sum(x[,2]))
})

sum.keep <- list()
for(i in 1:80){
  x <- dist.all[[i]]
  if(!is.na(x))sum.keep[[i]] <- c(sum(x[,1]),sum(x[,2]))
}

sk <- do.call(rbind,sum.keep)
plot(sk[,1],abs(sk[,2]),pch=19)

plot(distA[,1],distA[,2],ylab = 'Difference in Biomass',
     xlab = 'Difference in Composition',pch=19,cex=.5)

pog <- distA[which(abs(distA[,2])>5),]

hist(pog[,5])

pog <- distA[-which(distA[,2] > 10|distA[,2]< -10),]
plot(pog[,1],abs(pog[,2]),ylab = 'Difference in Biomass',
     xlab = 'Difference in Composition',pch=19,cex=.5)

### Sensitivity of biomass prediction to pollen proportion

Y.all <- do.call(rbind,Y.keep)
prop.all <- prop.table(as.matrix(Y.all),margin = 1)

pdf('pollen.prediction.sensitivity.pdf')
par(mfrow=c(2,2))
for(t in 1:22){
  Y.use <- prop.table(as.matrix(Y.keep[[1]][,3:24]),margin = 1)
  plot(biomassCI[[1]][2,Y.keep[[1]][,1]],Y.use[,t],
       ylim=c(0,(max(prop.all[,t+2])+.1*max(prop.all[,t+2]))),xlim=c(0,150),pch=19,cex=.5,
       ylab='Pollen Prop',xlab='Biomass Prediction (Mg/ha)')
  title(colnames(Y.use)[t])
  for(i in 1:62){
    if(!is.null(biomassCI[[i]])){
      Y.use <- prop.table(as.matrix(Y.keep[[i]][,3:24]),margin = 1)
      points(biomassCI[[i]][2,Y.keep[[i]][,1]],Y.use[,t],pch=19,cex=.5)
    }else{
      print(paste('not doing',i))
    }
  }
}
dev.off()

### For Paleon MIP
biomass.mean.df <- data.frame(lat=unlist(lat)[-c(35,39)],lon=unlist(long)[-c(35,39)],biomassMean = unlist(lapply(biomassCI,function(x){mean(x[2,1:11])}))[-c(2,35)])
write.csv(biomass.mean.df,file='biomass.means.csv')


biomass.settle <- last.index <-numeric(62)
for(i in 1:62){
  if(!is.null(name.keep[[i]])&i!=35){
    biomass.settle[i] <- x.meta[x.meta$site.name==name.keep[[i]],'SettleBiomassMean'][1]
    if(!is.na(biomass.settle[i])){
      last.index[i] <- biomassCI[[i]][2,1]
    }else{
      last.index[i] <- NA
    }
  }
}


### Creating correlation plot between settlement PLS estimates and settlement ReFAB estiamtes
last.index <- last.index[-c(which(is.na(last.index)),which(last.index==0))]
biomass.settle <- biomass.settle[-c(which(is.na(biomass.settle)),which(biomass.settle==0))]

pdf('settle.biomass.corr.pdf')
par(mfrow=c(1,1))
plot(biomass.settle,last.index,pch=19,xlab='Settlement Biomass Mean',ylab='Last Prediction Biomass Mean',xlim=c(0,155),ylim=c(0,155))
abline(a=0,b=1)
legend('topleft',paste('cor =', cor(x=biomass.settle,y=last.index)),pch=NA)
dev.off()

### Mapping sites by the last data point to make sure there is no spatial correlation in end of data
pdf('last_age_in_samp_map.pdf')
map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(unlist(long),unlist(lat), pch=19,
       cex=1.1,lwd=.2,col=rainbow(20,start = 0,end = .8)[unlist(stop.spot)])
legend('topright',as.character(1:20*100),col=rainbow(20,start = 0,end = .8),pch=rep(19,20))
title('Map Colored by Last Sample Years BP')
dev.off()





### Old code for showing RW and genpareto outputs together
blue       <- col2rgb("blue")
alphablue  <- rgb(blue[1], blue[2], blue[3], 75, max = 255)
orange       <- col2rgb("orange")
alphaorange  <- rgb(orange[1], orange[2], orange[3], 85, max = 255)

pdf('compare.RW.v.genPareto.Sigma0.12.pdf')
par(mfrow=c(2,2))
plot.new()
legend('center',c('GenPareto Sigma = .12','RW'),pch=c(19,19),col=c('darkorange','blue'))
for(i in 1:62){
  locn <- names(how.many)[i]
  site_number = unique(x.meta[x.meta$site.name == locn,1])
  keep.dataset.id <- unique(x.meta[x.meta$site.id==site_number,4])
  
  ten_count_use = ten.count[which(x.meta$site.id == site_number), ]
  
  Y = as.matrix(ten_count_use)
  
  sample_ages <- x.meta[x.meta[,1] == site_number, ]$age_bacon
  age_bins <- seq(minAge, maxAge, ageInterval)
  age_index <- as.matrix(as.numeric(
    cut(sample_ages, breaks = age_bins, labels=seq(1:(length(age_bins)-1)))
  ))
  
  tmp <- data.frame(cbind(age_index, Y))
  names(tmp)[1] <- 'age_index'
  
  Y2 <- aggregate(tmp, by = list(tmp$age_index), FUN = sum)
  
  if(!is.null(group) | FALSE){
    Y2 <- Y2[-group.mat[group,],]
  }
  
  Y <- as.matrix(Y2[ , -c(1,2)])
  age_index <- Y2[,1]
  
  plot(biomassCI[[i]][2,],col='blue',typ='l',xlim=c(100,0),
       ylim=c(0,150),lwd=2,main=names(how.many)[i])
  ciEnvelope(x = 1:100, yhi = biomassCI[[i]][3,],ylo=biomassCI[[i]][1,],col=alphablue)
  lines(biomassCI.save[[i]][2,],col='darkorange',lwd=2)
  ciEnvelope(x = 1:100, yhi = biomassCI.save[[i]][3,],ylo=biomassCI.save[[i]][1,],col=alphaorange)
  rug(x.meta[x.meta[,1]==site_number,]$age_bacon/100,lwd=2)
  
  if(length(out.list[[i]])!=0){
    points(age_index,seq(5, bMax, by = 2)[apply(out.list[[i]],2,which.max)])
  }
  
  mean.keep <- x.meta[x.meta$site.name == locn,'SettleBiomassMean'][1]
  sd.keep <- x.meta[x.meta$site.name == locn,'SettleBiomassSD'][1]
  points(0,mean.keep,col='blue',cex=1.5,pch=19)
  segments(x0=0,y0=mean.keep-sd.keep,x1 = 0,y1=mean.keep+sd.keep)
  #rug(control.pts[which(control.pts[,2]%in%keep.dataset.id),]$geo_age/100,lwd=3,col="red")
}
dev.off()



