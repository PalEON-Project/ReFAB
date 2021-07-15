library(dtwclust)
library(tidyverse)
library(maps)
library(dendextend)
library(vioplot)
library(ggfortify)
library(gridExtra)
library(plotrix)
library(fields)

#####
##### load biomass and metadata
#####

colors_tri <- c('#1bceff','#f09719','#c41b15')

agb.mat <- read.csv('median_biomass_plus_meta.csv')
agb.list <- split(x = agb.mat[,1:100],f = agb.mat$name) %>%
            lapply(.,as.numeric)

lon_mat <- agb.mat$lon[order(agb.mat$name)]
lat_mat <- agb.mat$lat[order(agb.mat$name)]

#####
##### Make clusters
#####
clusters <- tsclust(agb.list,type = 'hierarchical',k=3)
save(clusters,file='clusters.Rdata')
clusters@cluster

dend <- as.dendrogram(clusters)

dend <- color_labels(dend,k=3, col = cluster_colors[c(3,1,2)])
dend1 <- color_branches(dend,k=3,groupLabels = c('E','W','C'), col = cluster_colors[c(3,1,2)])

plot(dend1)

#####
##### Make cluster dendrogram and map
#####

#pdf('cluster_3.pdf',height=20,width = 18)
layout(matrix(c(1,2,3,3,4,4,5,5),4,2,byrow=T))

plot(dend)

plot(unlist(lon),unlist(lat),col=colors_tri[clusters@cluster],pch=19,cex=2)
maps::map('state',add=T)
text(lon,lat+.2, labels = names(clusters@cluster),cex=.5)

#par(mfrow=c(3,1))
for(i in 1:3){
  plot_me <- do.call(cbind,agb.list[which(clusters@cluster==i)])
  matplot(plot_me,xlim=c(100,0),typ='l',col=colors_tri[i],ylab = 'Biomass (Mg/ha)',xlab='Time',lwd=2)
  points(rowMeans(plot_me),type = 'l',lwd=4,col='blue')
}
#dev.off()

#####
##### divide pollen into cluster groups
#####

pollen_all <- read.csv('refab_pollen_counts.csv')

cluster_all <- numeric(nrow(pollen_all))
for(ii in 1:nrow(pollen_all)){
  if(any(which(pollen_all$site.name[ii]==names(agb.list)))){
    cluster_all[ii] <- clusters@cluster[which(pollen_all$site.name[ii]==names(agb.list))]
  }   
}

pollen_clusters <- list()
for(i in 1:3){
  pollen_clusters[[i]] <- as.data.frame(prop.table(as.matrix(pollen_all[which(cluster_all==i),2:23]),margin = 1))
}


#### sensitivity analysis
load('info.Rdata')
now_all <- do.call(rbind,info) #info comes from the average.biomass figure

#### adding in biome type pollen
arboreal <- rowSums(now_all[,-c(2,5,6,23:27)])
savanna <- rowSums(now_all[,c('QUERCUS','CARYA','ULMUS','POPULUS')])
forest <- 1-rowSums(now_all[,c('prairie',"other_herbs",'CYPERACE','QUERCUS','CARYA','ULMUS','POPULUS')])

conifer <- rowSums(now_all[,c('PINUSX','TSUGAX','PICEAX','LARIXPSEU','CUPRESSA','ABIES')])
decidious <- rowSums(now_all[,c('QUERCUS','BETULA','ALNUSX','OSTRYCAR','ULMUS','ACERX','FRAXINUX','POPULUS','JUGLANSX')])

now_all_all <- cbind(conifer,decidious,arboreal,savanna,forest,now_all)

clusters_assign <- clusters@cluster[as.character(now_all$site)]
now_all_df <- data.frame(now_all_all,clusters = clusters_assign)


now_all_df$site <- as.numeric(as.factor(now_all_df$site))

#gut checks
plot(now_all_df[now_all_df$clusters==2,'lon'],now_all_df[now_all_df$clusters==2,'lat'])
maps::map('state',add=T)
plot(now_all_df[now_all_df$clusters==1,'age'],now_all_df[now_all_df$clusters==1,'max_est'],xlim=c(100,0))
plot(now_all_df[now_all_df$site==1,'age'],now_all_df[now_all_df$site==1,'max_est'],xlim=c(100,0))
plot(now_all_df[now_all_df$site==1,'age'],now_all_df[now_all_df$site==1,'TSUGAX'],xlim=c(100,0))


#sensitivity calculating
pollen_save <- biomass_save <- cluster_save <- site_dat_save <- site_lat_save <- site_lon_save <- list()
ages_set <- matrix(c(100,80,81,50,51,0),3,2,byrow = T)
for(ii in 1:77){
  site_dat <- now_all_df[now_all_df$site==ii,]
  pollen_save[[ii]] <- biomass_save[[ii]] <- cluster_save[[ii]] <- site_dat_save[[ii]] <- site_lat_save[[ii]] <- site_lon_save[[ii]] <- list()
  for(pp in 1:27){
    pollen_save[[ii]][[pp]] <- biomass_save[[ii]][[pp]] <- cluster_save[[ii]][[pp]] <-site_dat_save[[ii]][[pp]] <- site_lat_save[[ii]][[pp]] <- site_lon_save[[ii]][[pp]] <- list()
    ages <- site_dat[,'age']
    for(tt in 1:3){
      row_get <- which(ages <= ages_set[tt,1] &ages>=ages_set[tt,2])
      
      row_order <- order(ages[row_get],decreasing = T)
      
      biomasses <- site_dat[row_get[row_order],'max_est']
      pollens <- site_dat[row_get[row_order],pp]
  
      dist_mat <- outer(ages[row_get[row_order]],ages[row_get[row_order]],FUN='-')
      dist_mat_biomass <- outer(biomasses,biomasses,FUN='-')
      dist_mat_pollen <- outer(pollens,pollens,FUN='-')
    
      #dist_mat[lower.tri(dist_mat)] <- dist_mat_biomass[lower.tri(dist_mat_biomass)] <- dist_mat_pollen[lower.tri(dist_mat_pollen)] <- NA
      
      dist_mat_pollen[abs(dist_mat)>10 | abs(dist_mat)<1] <- dist_mat_biomass[abs(dist_mat)>10| abs(dist_mat)<1] <- NA
     
      remove <- c(which(lower.tri(dist_mat_pollen,diag = T)),which(dist_mat_pollen <= .03 & dist_mat_pollen >= -.03))#which(duplicated(c(dist_mat_pollen)))
      
      pollen_save[[ii]][[pp]][[tt]] <- dist_mat_pollen[-remove]
      biomass_save[[ii]][[pp]][[tt]] <- dist_mat_biomass[-remove]
      cluster_save[[ii]][[pp]][[tt]] <- rep(site_dat$cluster[1],length(dist_mat_biomass[-remove]))
      site_dat_save[[ii]][[pp]][[tt]] <- rep(site_dat$site[1],length(dist_mat_biomass[-remove]))
      site_lat_save[[ii]][[pp]][[tt]] <- rep(site_dat$lat[1],length(dist_mat_biomass[-remove]))
      site_lon_save[[ii]][[pp]][[tt]] <- rep(site_dat$lon[1],length(dist_mat_biomass[-remove]))
      #plot(dist_mat_pollen[-remove],dist_mat_biomass[-remove])
    }
  }
}
time_names <- c('early','mid','late')
colors <- colors_tri#c('black','darkred','darkgreen')
colors_tri <- cluster_colors

rlm_coeff <- matrix(NA,27,3)

#### Overall 
pdf('pol_biom_diff_scatters_limited.pdf',height = 20,width=20)
#par(mfrow=c(1,1))
#plot(allpol_df$lon,allpol_df$lat,col=allpol_df$clusters,pch=19,cex=3)
#maps::map('state',add=T)
par(mfrow=c(3,3),mar = c(3,3,2,0),oma = rep(4,4))
layout(matrix(1:9,3,3))
  for(ss in 3:5){
    all_pollen <- unlist(do.call(c, lapply(
      pollen_save,
      FUN = function(x) {
        x[[ss]]
      }
    )))
    all_biomass <- unlist(do.call(c, lapply(
      biomass_save,
      FUN = function(x) {
        x[[ss]]
      }
    )))
    
    all_cluster <- unlist(do.call(c, lapply(
      cluster_save,
      FUN = function(x) {
        x[[ss]]
      }
    )))
    for(cc in 1:3){
      Lab.palette <- colorRampPalette(c('white',colors[cc],'black'), space = "Lab")
      subset <- which(all_cluster==cc)#1:length(all_pollen)#runif(1000,1,length(all_pollen))
      
      
      sd_away <- 5*sd(all_pollen[subset][which(!is.na(all_pollen[subset]))])
      
      keep_idx <- which(all_pollen[subset] < sd_away & all_pollen[subset] > -sd_away)
      
      if(!any(all_pollen[subset][keep_idx])) next()
      
      plot(all_pollen[subset][keep_idx],
                    all_biomass[subset][keep_idx],
                    main = ,cex.axis=3,
                    #ylim = c(-250,250),
                    #xlim = c(-.5,.5),
                    xlab = 'Pollen Proportion Difference',
                    ylab = 'Biomass Difference',col = colors_tri[cc])
      title(paste(colnames(now_all_df)[ss]),cex.main=3)
      
      if(length(all_pollen[subset][keep_idx])<500) next()
      
      rlm_obj<- MASS::rlm(all_pollen[subset][keep_idx],
          all_biomass[subset][keep_idx])
      rlm_coeff[ss,cc] <- rlm_obj$coefficients
      abline(a=0,b = rlm_obj$coefficients,col='blue')
      mtext(text = paste('slope =',signif(rlm_obj$coefficients,digits = 3),side = 3),outer = F,line=-1)
      
      
    }
  }

mtext('Pollen Proportion Difference',outer = T,side = 1, line = 2)
mtext('Biomass Estimate Difference',outer = T,side = 2, line = 2)
dev.off()

### barplots
pdf('slope_barplots.pdf')
rownames(rlm_coeff) <- colnames(now_all_df)[1:25]
par(mfrow=c(3,1))
for(cc in 1:3){ 
  barplot(rlm_coeff[,cc],col=fields::tim.colors(25),beside=T,las=2,ylim=c(-150,150))
  title(paste('Cluster',cc))
}
dev.off()

cols2 <- colorRampPalette(c("blue","white","red"))(20)

rownames(rlm_coeff) <- c('Conifer','Decidious','Arboreal','Savanna',
                         'Forest', 'Pinus','Prairie','Quercus',
                         'Betula','Other Herb.','Cyperace',
                         'Alnus','Ostrycar','Ulmus',
                         'Tsuga','Picea','Acer',
                         'Fraxinus','Populus','Cupressa',
                         'Other Trees','Larix','Fagus',
                         'Carya','Tilia','Abies',
                         'Juglans')

###
### all together figure
###

settle <- read.csv('paleon_total_biomass_settlement.csv')

breaks <-  c(seq(0,50,10),seq(75,250,25),435)
colors <- rev(terrain.colors(length(breaks)-1))
data_binned <-  cut(settle$TotalBiomass_Mgha, breaks, include.lowest = TRUE, labels = FALSE)
legendName <- c("Biomass (Mg/ha)")#paste0("Biomass at Age = ",age_slice, " BP"
breaklabels <- apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1,  function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })

cols2 <-brewer.pal(7,'BrBG')#c('slateblue4',colorRampPalette(c("lightblue3","lemonchiffon","lightcoral"))(6),'red')# RColorBrewer::brewer.pal(11,'BrBG')
colors_tri <- c('navy','royalblue','magenta3')

dend <- as.dendrogram(clusters)
dend <- color_branches(dend,k=3, col = cluster_colors[c(3,2,1)],groupLabels = T)

# %>% 
#   dendextend::rotate(c(1:77))%>% #rotate to match labels new order
#   plot(main = "Rotated tree\n based on labels")
 dend2 <- click_rotate(rev(dend1))
 
 
#### pollen over time in a line calculation
 
clustline <- now_all_df %>%
  mutate(bins = cut(age,breaks=seq(0,100,10))) %>%
  group_by(clusters,bins) %>%
  summarise_all(list(mean))

clustmelt <- melt(as.data.frame(clustline[,1:29]),id.vars = c('clusters','bins'))

clustmelt1 <- clustmelt[clustmelt$clusters==1 & clustmelt$variable %in% c('forest'),]
clustmelt2 <- clustmelt[clustmelt$clusters==2 & clustmelt$variable %in% c('TSUGAX'),]
clustmelt3 <- clustmelt[clustmelt$clusters==3 & clustmelt$variable %in% c('PINUSX','conifer','FAGUS'),]

cm4 <- rbind(clustmelt1,clustmelt2,clustmelt3)


clustline_ages <- now_all_df %>%
  mutate(ages = cut(age,breaks=seq(0,100,1))) %>%
  group_by(clusters,ages) %>%
  summarise_all(list(mean))

clustmelt_ages <- melt(as.data.frame(clustline_ages[,1:29]),id.vars = c('clusters','ages'))

clustmelt1_ages <- clustmelt_ages[clustmelt_ages$clusters==1 & clustmelt_ages$variable %in% c('forest'),]
clustmelt2_ages <- clustmelt_ages[clustmelt_ages$clusters==2 & clustmelt_ages$variable %in% c('PINUSX','QUERCUS','TSUGAX'),]
clustmelt3_ages <- clustmelt_ages[clustmelt_ages$clusters==3 & clustmelt_ages$variable %in% c('PINUSX','conifer','QUERCUS','FAGUS'),]

cm5 <- rbind(clustmelt1_ages,clustmelt2_ages,clustmelt3_ages)


ggplot(data = cm4)+
  geom_tile(aes(x = bins,y = variable,fill = value)) +
  #scale_fill_brewer(palette = 1) + 
  facet_wrap(facets = vars(clusters),nrow=3)

load('~/Dropbox/ReFAB_outputs/split_calib_dat_v3.0.Rdata')

common_names <- c('Conifer','Decidious','Arboreal','Savanna','Forest','Pine','Prairie',
                  'Oak','Birch','Other Herb.','Sedge','Alder','Ironwood','Elm','Hemlock',
                  'Spruce','Maple','Ash','Poplar','Cypress','Other Trees','Larch','Beech',
                  'Hickory','Basswood','Fir','Walnut')

names_dat <- cbind(rownames(rlm_coeff),common_names)

rlm_order <- c(3,1,2,4,5,seq(6,27,1)[order(names_dat[6:27,2])])

pdf(paste0(Sys.Date(),'cluster_fig_refab_manuscript_2.pdf'),height=11,width = 9)
layout(matrix(c(1, 1, 2, 2,
                1, 1, 2, 2,
                3, 3, 3, 6,
                3, 3, 3, 6,
                4, 4, 4, 6,
                4, 4, 4, 6,
                5, 5, 5, 6,
                5, 5, 5, 6), 8, 4, byrow = T))
par(mar=c(1.5,5,1,3),oma=c(9,3,3,3))

#####
##### 2/3s to 1/3s fit R2 ### decided to remove 12/2
#####
# plot(biomass2_3, apply(samps2_3,2,FUN = quantile,.5),
#      xlim=c(0,bMax), ylim=c(0,bMax), pch=21,
#      xlab="True Biomass (Mg/ha)",
#      ylab="Predicted Mean Biomass (Mg/ha)")
# abline(a=0,b=1)
# lm.mod <- lm(apply(samps2_3,2,FUN = quantile,.5)~biomass2_3+0)
# abline(lm.mod,lty=2)
# 
# lm.mod.out <- lm(apply(samps1_3,2,FUN = quantile,.5)~biomass1_3+0)
# abline(lm.mod.out,lty=2,col='red')
# 
# points(biomass1_3,colMeans(samps1_3, na.rm = T),
#        col='red',pch=21)
#mtext(paste("r2-twothirds = ",signif(summary(lm.mod)$r.squared,digits=2),'r2-onethird = ',signif(summary(lm.mod.out)$r.squared,digits = 2)))

#R2training = signif(summary(lm.mod)$r.squared,digits=3)
#R2testing = signif(summary(lm.mod.out)$r.squared,digits = 3)
# 
# R2training <- signif(give_me_R2(preds = apply(samps2_3,2,FUN = quantile,.5),
#                                 actual =  biomass2_3),digits = 3)
# R2testing <- signif(give_me_R2(preds = apply(samps1_3,2,FUN = quantile,.5),
#                                actual =  biomass1_3),digits = 3)
# 
# mse_train <- signif(mse_func(preds = apply(samps2_3,2,FUN = quantile,.5),
#                              actual =  biomass2_3),digits = 3)
# mse_test <- signif(mse_func(preds = apply(samps1_3,2,FUN = quantile,.5),
#                             actual =  biomass1_3),digits = 3)
# 
# arrows(x0 = biomass2_3, y0 = apply(samps2_3,2,FUN = quantile,.05),
#        x1 = biomass2_3, y1 = apply(samps2_3,2,FUN = quantile,.975),
#        code = 0, lwd=1)
# arrows(x0 = biomass1_3, y0 = apply(samps1_3,2,FUN = quantile,.05),
#        x1 = biomass1_3, y1 = apply(samps1_3,2,FUN = quantile,.975),
#        code = 0, lwd=1, col = 'red')
# legend(
#   'bottomright',
#   c(paste0('Training (R2 = ', R2training,')'),#,', MSE = ',mse_train
#     paste0('Testing (R2 = ', R2testing,')')),#', MSE = ',mse_test,
#   pch=1,col=c('black','red'))



plot((dend2), leaflab = "none",lwd=3,las=2,cex.axis=1.5)
mtext('Height',side=2,line=5,cex=1.5)

plot(settle$lon,
     settle$lat,
     col = 'white',#colors[data_binned],
     ylim = c(41, 49),ylab='',xlab='',las=2,yaxt='n',xaxt='n',cex=.5)
#axis(side=1)

points(
  unlist(lon_mat),
  unlist(lat_mat),
  col = colors_tri[clusters@cluster],
  pch = 19,cex = 1.25
)
maps::map('state', add = T)

# plot.new()
# legend('left',breaklabels[1:7],pch=19,col=colors[1:7],cex=1.1,title = 'PLS Biomass (Mg/ha)')
# legend('right',breaklabels[8:14],pch=19,col=colors[8:14],cex=1.1,title = 'PLS Biomass (Mg/ha)')

# 
# colorbar.plot(
#   1,
#   50,
#   1:20,
#   strip.width = .5,
#   strip.length = 8.5,
#   col = colvals(length(breakvals) - 1),
#   adj.x = 0
# )
# axis(side=1,at=seq(1.25,7.75,.75),
#      labels = signif(seq(-1,1,length.out = length(seq(1.25,7.75,.75))),
#                      digits=2),las=1)

#par(mfrow=c(3,1))
for (i in c(1,2,3)) {
  plot_me <- do.call(cbind, agb.list[which(clusters@cluster == i)])
  matplot(
    plot_me[2:100,],
    xlim = c(100, -9),
    xaxt = 'n',
    yaxt = 'n',
    typ = 'l',
    col = adjustcolor(colors_tri[i],alpha.f = .5),
    ylab = '',
    xlab = '',
    lwd = 2,
    lty = 1,
    las =2,
    cex.axis = 1.5,
    ylim =c(0,320)
  )
  
  text(x = -7,y=125,labels= c('W','C','E')[i],cex=4,col = colors_tri[i])
  
  ciEnvelope(x=1:100,ylo = sd_clusts[i,],
             yhi = sd_clusts1[i,],
             col=adjustcolor('gray',alpha.f = .75))
  points(mean_clusts[i,],type = 'l',lty=2)
  points(sd_clusts[i,],type = 'l',lty=1)
  points(sd_clusts1[i,],type = 'l',lty=1)
  
  axis(2,at = seq(0,250,50),las=2,cex.axis = 1.5)
  
  abline(h = 250,col='black')
  abline(v=1,col='black')
  
  linesmat <- cast(cm4[cm4$clusters==i,],variable~bins)
  for(nnn in 1:nrow(linesmat)){
    linesvals <- (as.numeric(1*linesmat[nnn,-1]/max(linesmat[nnn,-1])))
    breakvals <- seq(-1,1,.01)
    cutvals <- cut(linesvals,breaks = breakvals,labels=F)
    colvals <- colorRampPalette(c('white','antiquewhite','peachpuff','#8db27d','#3b4c41'))
    
    if(nrow(linesmat)==1){ 
      y_do = 250+2*20 
      }else{
      y_do = 250+nnn*20 
    }
    
    segments(x0 = c(2,seq(10,90,10)),x1 = seq(10,100,10),y1 = y_do,y0=y_do,col=colvals(length(breakvals)-1)[cutvals],lwd = 10)
  }
  
  names_segs <- list(w = list('Forest'),c = list('Hemlock'),e = list(c('Conifer','Pine','Beech')))
  
  if(nrow(linesmat)==1){ 
    text(x=0,y = 250+2*20,unlist(names_segs[[i]]),cex=1.2,adj=c(0,.5))
  }else{
    text(x=0,y = 250+1:nrow(linesmat)*20,unlist(names_segs[[i]]),cex=1.2,adj=c(0,.5)) 
  }
  
  #if(i==1) colorbar.plot(6, 90, 1:100)
  
  # points(rowMeans(plot_me),
  #        type = 'l',
  #        lwd = 4,
  #        col = 'blue')
  if(i==2)     mtext('Biomass (Mg/ha)',side = 2,line = 4,cex = 1.5)
}
axis(side = 1, at = seq(0,100,10),labels = seq(0,10,1),cex.axis = 1.5)
mtext('Age (cal ka BP)',side = 1,line = 4,cex = 1.5)

image(
  t(rlm_coeff[rev(rlm_order),c(1,2,3)]),
  col = cols2,
  breaks = c(-200,seq(-100, 100, length.out = length(cols2) - 1),200),
  xaxt = 'n',
  yaxt = 'n'
)

axis(
  side = 2,
  at = seq(0, 1, length.out = 27),
  labels = rev(common_names[rlm_order]),
  las = 2,
  cex.axis = 1.5
)
for (i in 1:3)
  axis(
    side = 1,
    at = c(0, .5, 1)[i],
    labels = paste('Cluster', c('W','C','E'))[i],
    col.axis = colors_tri[c(1,2,3)][i],
    las = 2,
    cex.axis = 1.5
  )

image.plot(
  t(rlm_coeff[27:1,]),
  col = cols2,
  breaks = c(-200,seq(-100, 100, length.out = length(cols2) - 1),200),
  xaxt = 'n',
  yaxt = 'n',
  legend.width = 4,
  legend.mar = 2,
  #legend.lab = 'Slope',
  legend.cex = 1.5,
  legend.only = T
)

dev.off()

pdf('cluster_pollen_prop_legend.pdf')
nbins = 70
plot(
  breakvals[seq(1,length(breakvals),length.out = nbins)],
  rep(0, nbins),
  col = colvals(length(breakvals))[seq(1,length(breakvals),length.out = nbins)],
  pch = 15,
  cex = 4,
  ylim = c(-.01, .01),
  xlim = c(.03,1),
  bty="n",yaxt='n',ylab=NA
)
title('Cluster-normalized average pollen proportion')
dev.off()

# 
# lattice::levelplot(
#   t(rlm_coeff[25:1, ]),
#   at = seq(-200, 200, 25),
#   col.regions = cols2,
#   ylab = '',
#   xlab = '',
#   xaxt = 'n',
#   scales = list(
#     y = list(cex = 1.25),
#     x = list(
#       at = 1:3,
#       labels = paste('Cluster', 1:3),
#       col = colors_tri,
#       rot = 90,
#       cex = 1.5
#     )
#   )
# )













### maps

early <- now_all_df[which(now_all_df$age>80),]
mid <- now_all_df[which(now_all_df$age<=80&now_all_df$age>=50),]
late <- now_all_df[which(now_all_df$age<50),]

early_avg <- stats::aggregate(
  early[,1:25],
  by = list(site = early$site),
  data = early,
  FUN = "mean",
  na.rm = T
)

mid_avg <- stats::aggregate(
  mid[,1:25],
  by = list(site = mid$site),
  data = mid,
  FUN = "mean",
  na.rm = T
)

late_avg <- stats::aggregate(
  late[,1:25],
  by = list(site = late$site),
  data = late,
  FUN = "mean",
  na.rm = T
)

lat <- stats::aggregate(
  now_all_df$lat,
  by = list(site = now_all_df$site),
  data = now_all_df,
  FUN = "mean",
  na.rm = T
)
lon <- stats::aggregate(
  now_all_df$lon,
  by = list(site = now_all_df$site),
  data = now_all_df,
  FUN = "mean",
  na.rm = T
)

clust_ordered <- stats::aggregate(
  now_all_df$clusters,
  by = list(site = now_all_df$site),
  data = now_all_df,
  FUN = "mean",
  na.rm = T
)

color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

n_split <- 50
colors <- colorRampPalette(c('darkblue','blue', 'white', 'red','darkred'))(n_split)
my.colors <-colorRampPalette(c('darkblue','blue', 'white', 'red','darkred'))

pdf('pollen_diff_maps.pdf',width=10,height=10)
par(mfrow = c(4,4), mar = c(0, 0, 3, 0))
#mat_lay <- matrix(1:24,4,6)
#mat_lay[,5] <- 0
#mat_lay[,6] <- 17

#layout(mat_lay)

for (ss in 2:26) {
  diff_avg1 <-  mid_avg[, ss] - early_avg[, ss]
  diff_avg2 <-  late_avg[, ss] - mid_avg[, ss]
  
  range_use <- max(abs(range(diff_avg1,diff_avg2)))
  breaks <- seq(-range_use, range_use, length.out = n_split)
  cuts1 <- cut(diff_avg1, breaks = breaks, labels = F)
  plot(lon[, 2],
       lat[, 2],
       col = adjustcolor(viridis(3)[clust_ordered[, 2]],alpha.f = .5),
       bg = colors[cuts1],
       pch = 21,
       xaxt = 'n',
       yaxt = 'n',
       xlab = NA,
       ylab = NA)
  maps::map('state', add = T)
  title(paste('Mid - Early', colnames(now_all_df)[ss - 1]))
  mtext(paste('max pollen prop =',signif(range_use,digits=3)),side = 3,outer = F,line=0)
  
  cuts2 <- cut(diff_avg2, breaks = breaks, labels = F)
  plot(lon[, 2],
       lat[, 2],
       col = adjustcolor(viridis(3)[clust_ordered[, 2]],alpha.f = .5),
       bg = colors[cuts2],
       pch = 21,
       xaxt = 'n',
       yaxt = 'n',
       xlab = NA,
       ylab = NA)
  maps::map('state', add = T)
  title(paste('Late - Mid', colnames(now_all_df)[ss - 1]))
  mtext(paste('max pollen prop =',signif(range_use,digits=3)),side = 3,outer = F,line=0)
  
  if(FALSE) { #ss == 9 | ss == 17 | ss == 25 | 
    z=matrix(1:50,nrow=1)
    x=1
    y=breaks # supposing 3 and 2345 are the range of your data
    image(x,y,z,col=my.colors(50),axes=FALSE,xlab="",ylab="")
    axis(2)
  }
}

dev.off()

#### By Site
lm_obj <- array(NA,dim = c(77,25,3))

pdf('scatter_by_site.pdf')
for(ii in 1:77){
  
  par(mfrow=c(1,1))
  #plot(now_all_df$lon,now_all_df$lat,col=now_all_df$clusters,pch=1,cex=3)
  #maps::map('state',add=T)
  #text(now_all_df[now_all_df$site==ii,'lon'][1],now_all_df[now_all_df$site==ii,'lat'][1],pch=1,labels = now_all_df[now_all_df$site==ii,'site'][1])
  
  plot(allpol_df[allpol_df$site==ii,'age'],allpol_df[allpol_df$site==ii,'biomass'],
       xlab = 'Age',
       ylab = 'Max Lik Biomass', xlim=c(10000,0),ylim=c(0,250))
  
  par(mfrow=c(3,3))
  for(pp in 1:25){
    for(tt in 1:3){
      if(length(c(pollen_save[[ii]][[pp]][[tt]]))<5) next()
      plot(c(pollen_save[[ii]][[pp]][[tt]]),c(biomass_save[[ii]][[pp]][[tt]]),
           col=cluster_save[[ii]][[pp]][[tt]][1], pch = 1,
           xlab = 'Pollen Proportion Difference',ylab='Biomass Difference')
      if(max(c(pollen_save[[ii]][[pp]][[tt]]),na.rm=T)>.03) {
        keep_lm <- lm(c(biomass_save[[ii]][[pp]][[tt]])~c(pollen_save[[ii]][[pp]][[tt]])+0)
        sum_lm <- summary(keep_lm)
        if(!is.na(sum_lm$r.squared)){
          if(!is.null(sum_lm$r.squared) & sum_lm$r.squared>.25) lm_obj[ii,pp,tt] <- sum_lm$r.squared
          }
        }
      if(!is.na(lm_obj[ii,pp,tt])) abline(a=0,b= keep_lm$coefficients,col='blue')
      title(paste(ii,time_names[tt],colnames(now_all_df)[pp]),cex.main=1)
    }
  }
}
dev.off()


max_r2 <- list()
for(ee in 1:3){
  max_r2[[ee]] <- apply(lm_obj[,,3],1,FUN=function(x) {
    y <- which.max(x)
    if(!any(y)) {y <- NA}else {y<- which.max(x)}
  })
}

tab_max <- table(do.call(c,max_r2))
do_nums <- as.numeric(names(tab_max)[which(tab_max>3)])


par(mfrow=c(2,2))
for(cc in 1:3){
plot(now_all_df$lon,now_all_df$lat,col='white',pch=1,cex=3)
maps::map('state',add=T)

for(i in 1:77) {
  if(any(which.max(lm_obj[i,,cc]))){
    if(which.max(lm_obj[i,,cc])%in%do_nums){
  points(
  now_all_df[now_all_df$site==i,'lon'],
  now_all_df[now_all_df$site==i,'lat'],
  #col = now_all_df$clusters,
  pch = 19,
  cex = 1,
  col = fields::tim.colors(length(do_nums))[which.max(lm_obj[i,,cc])]
)
    }
  }
}

}
plot.new()
legend('center',
       colnames(now_all_df)[do_nums],
       pch = 19,
       col = fields::tim.colors(length(do_nums)),
       cex=1)



pdf('aggregated_slopes.pdf')
par(mfrow=c(3,3))
for(pp in 1:25){
  for(tt in 1:3){
    all_pollen <- do.call(c, lapply(
      pollen_save,
      FUN = function(x) {
        x[[pp]][[tt]]
      }
    ))
    all_biomass <- do.call(c, lapply(
      biomass_save,
      FUN = function(x) {
        x[[pp]][[tt]]
      }
    ))
    
    all_cluster <- do.call(c, lapply(
      cluster_save,
      FUN = function(x) {
        x[[pp]][[tt]]
      }
    ))
    plot(c(all_pollen),c(all_biomass),
         pch = 1,
         xlab = 'Pollen Proportion Difference',ylab='Biomass Difference')
    for(ii in 1:77){
      if(!is.na(lm_obj[ii,pp,tt])) abline(a=0,b=lm_obj[ii,pp,tt],col=now_all_df[now_all_df$site==ii,'clusters'][1])
    } 
    title(paste(time_names[tt],colnames(now_all_df)[pp]),cex.main=1)
  }
}
dev.off()

pdf('boxplot_slopes.pdf',height=15,width=15)
par(mfrow=c(3,3))
for(cc in 1:3){
clust_get <- unique(now_all_df[now_all_df$cluster==cc,'site'])
for(tt in 1:3){
  
  sort_by <- order(abs(apply(lm_obj[clust_get,,tt],2,quantile,.5,na.rm=T)))
  boxplot((lm_obj[clust_get,sort_by,tt]),xaxt='n',col=cc)
  abline(h=0,col='blue')
  axis(1,at = 1:25,labels = colnames(now_all_df)[sort_by],las=2)
}
}
dev.off()

#### averaged over site
library(nimble)
cluster_keep <- list()
cor_keep <- array(NA,dim = c(77,25,3))
for(ii in 1:77){
  use <- allpol_df[as.numeric(as.factor(allpol_df$site))==ii,]

  for(ss in 1:22){
    early <- which(use$age>8000 & use$age<10000 & use[,ss]!=0)
    mid <- which(use$age<8000 & use$age>5000 & use[,ss]!=0)
    late <- which(use$age<5000 & use[,ss]!=0)
    
    if(any(early) & length(early)>2) cor_keep[ii,ss,1] <- lm(use[early,'biomass']~logit(use[early,ss]))$coefficients[2]
    if(any(mid) & length(mid)>2) cor_keep[ii,ss,2] <- lm(use[mid,'biomass']~logit(use[mid,ss]))$coefficients[2]#cor(use[mid,ss],use[mid,'biomass'])
    if(any(late) & length(late)>2) cor_keep[ii,ss,3] <- lm(use[late,'biomass']~logit(use[late,ss]))$coefficients[2]#cor(use[late,ss],use[late,'biomass'])
    
  }
  
  early <- which(use$age>8000 & use$age<10000)
  mid <- which(use$age<8000 & use$age>5000)
  late <- which(use$age<5000)
  
  cor_keep[ii,23,1] <- cor(rowSums(use[early,-c(2,5,6,23:26)]),use[early,'biomass'])
  cor_keep[ii,23,2] <- cor(rowSums(use[mid,-c(2,5,6,23:26)]),use[mid,'biomass'])
  cor_keep[ii,23,3] <- cor(rowSums(use[late,-c(2,5,6,23:26)]),use[late,'biomass'])
  
  cor_keep[ii,24,1] <- cor(rowSums(use[early,c('QUERCUS','CARYA','ULMUS','POPULUS')]),use[early,'biomass'])
  cor_keep[ii,24,2] <- cor(rowSums(use[mid,c('QUERCUS','CARYA','ULMUS','POPULUS')]),use[mid,'biomass'])
  cor_keep[ii,24,3] <- cor(rowSums(use[late,c('QUERCUS','CARYA','ULMUS','POPULUS')]),use[late,'biomass'])
  
  cor_keep[ii,25,1] <- cor(1-rowSums(use[early,c('prairie',"other_herbs",'CYPERACE','QUERCUS','CARYA','ULMUS','POPULUS')]),use[early,'biomass'])
  cor_keep[ii,25,2] <- cor(1-rowSums(use[mid,c('prairie',"other_herbs",'CYPERACE','QUERCUS','CARYA','ULMUS','POPULUS')]),use[mid,'biomass'])
  cor_keep[ii,25,3] <- cor(1-rowSums(use[late,c('prairie',"other_herbs",'CYPERACE','QUERCUS','CARYA','ULMUS','POPULUS')]),use[late,'biomass'])
  
  cluster_keep[ii] <- use$cluster[1]
}

par(mfrow=c(3,3),mar=rep(2.5,4),oma=c(4,4,1,1))
for(cc in 1:3){
  for(tt in 1:3){
data <- (data.frame(cor_keep[which(unlist(cluster_keep)==cc),,tt]))

colnames(data) <- c('Pinus','Prairie','Quercus',
                    'Betula','Other Herb.','Cyperace',
                    'Alnus','Ostrycar','Ulmus',
                    'Tsuga','Picea','Acer',
                    'Fraxinus','Populus','Cupressa',
                    'Other Trees','Larix','Fagus',
                    'Carya','Tilia','Abies',
                    'Juglans','Arboreal', 'Savanna',
                    'Forest')

data <- data[,c('Pinus','Quercus',
                          'Betula','Other Herb.','Cyperace',
                          'Alnus','Ostrycar','Ulmus',
                          'Tsuga','Picea','Acer',
                          'Fraxinus','Populus','Cupressa',
                          'Other Trees','Larix','Fagus',
                          'Carya','Tilia','Abies',
                          'Juglans','Arboreal','Prairie',
                          'Savanna','Forest')]

for(ss in 1:ncol(data)){
  if(length(which(is.na(data[,ss])))>10){
    data[,ss] <- 0
  }
}

vioplot(as.data.frame(data),
        las = 2,
        col = cc,
        beside = T,
        cex.names = .75,
        ylim = c(-.2,.2))
abline(h = c(.2,0,-.2),lty=2)

if(cc==1){
  title(c('Early','Mid','Late')[tt])
}
if(cc==2&tt==1){
  mtext('Pollen Proportion Correlation with Biomass Mean Estimate',side=2,outer = T)
}
  }
}

#### not averaged
sub_allpol_df <- allpol_df[which(allpol_df$cluster==1),]
clust1 <- (apply(sub_allpol_df[,1:22],2,FUN = function(x) cor(x, sub_allpol_df$biomass)))

sub_allpol_df <- allpol_df[which(allpol_df$cluster==2),]
clust2 <- (apply(sub_allpol_df[,1:22],2,FUN = function(x) cor(x, sub_allpol_df$biomass)))

sub_allpol_df <- allpol_df[which(allpol_df$cluster==3),]
clust3 <- (apply(sub_allpol_df[,1:22],2,FUN = function(x) cor(x, sub_allpol_df$biomass)))

all_clust <- rbind(clust1, clust2, clust3)
barplot(all_clust[,order(colSums(all_clust))],las=2,col=1:3,beside=T, cex.names = .75)

colnames(allpol_df)[1:22] <- c('Pinus','Prairie','Quercus',
                         'Betula','Other Herb.','Cyperace',
                         'Alnus','Ostrycar','Ulmus',
                         'Tsuga','Picea','Acer',
                         'Fraxinus','Populus','Cupressa',
                         'Other Trees','Larix','Fagus',
                         'Carya','Tilia','Abies',
                         'Juglans')


#### correlation plots
pdf('corrs_over_time.pdf',height=10,width = 10)
layout(matrix(1:9,3,3))
par(mar=rep(1,4),oma=c(4,4,1,1))
#10k-8k
sub_allpol_df <- allpol_df[which(allpol_df$cluster==1 & allpol_df$age<10000 & allpol_df$age>8000),]
clust1 <- (apply(sub_allpol_df[,1:22],2,FUN = function(x) cor(x, sub_allpol_df$biomass)))

sub_allpol_df <- allpol_df[which(allpol_df$cluster==2 & allpol_df$age<10000 & allpol_df$age>8000),]
clust2 <- (apply(sub_allpol_df[,1:22],2,FUN = function(x) cor(x, sub_allpol_df$biomass)))

sub_allpol_df <- allpol_df[which(allpol_df$cluster==3 & allpol_df$age<10000 & allpol_df$age>8000),]
clust3 <- (apply(sub_allpol_df[,1:22],2,FUN = function(x) cor(x, sub_allpol_df$biomass)))

all_clust <- rbind(clust1, clust2, clust3)
barplot(
  clust1,
  las = 2,
  col = 1,
  beside = T,
  cex.names = .75,
  main = 'Early Holocene 10ka - 8ka',
  ylab = 'Cluster 1 Corrleation',
  ylim = c(-.25,.65)
)
abline(h = c(.3),lty=2)
barplot(
  clust2,
  las = 2,
  col = 2,
  beside = T,
  cex.names = .75,
  ylab = 'Cluster 2 Corrleation',
  ylim = c(-.25,.65)
)
abline(h = c(.3),lty=2)
barplot(
  clust3,
  las = 2,
  col = 3,
  beside = T,
  cex.names = .75,
  ylab = 'Cluster 3 Corrleation',
  ylim = c(-.25,.65)
)
abline(h = c(.3),lty=2)

#8k-5k
sub_allpol_df <- allpol_df[which(allpol_df$cluster==1 & allpol_df$age<8000 & allpol_df$age>5000),]
clust1 <- (apply(sub_allpol_df[,1:22],2,FUN = function(x) cor(x, sub_allpol_df$biomass)))

sub_allpol_df <- allpol_df[which(allpol_df$cluster==2 & allpol_df$age<8000 & allpol_df$age>5000),]
clust2 <- (apply(sub_allpol_df[,1:22],2,FUN = function(x) cor(x, sub_allpol_df$biomass)))

sub_allpol_df <- allpol_df[which(allpol_df$cluster==3 & allpol_df$age<8000 & allpol_df$age>5000),]
clust3 <- (apply(sub_allpol_df[,1:22],2,FUN = function(x) cor(x, sub_allpol_df$biomass)))

all_clust <- rbind(clust1, clust2, clust3)
barplot(
  clust1,
  las = 2,
  col = 1,
  beside = T,
  cex.names = .75,
  main = 'Mid Holocene 8ka - 5ka',
  ylim = c(-.25,.65)
)
abline(h = c(.3),lty=2)
barplot(
  clust2,
  las = 2,
  col = 2,
  beside = T,
  cex.names = .75,
  ylim = c(-.25,.65)
)
abline(h = c(.3),lty=2)
barplot(
  clust3,
  las = 2,
  col = 3,
  beside = T,
  cex.names = .75,
  ylim = c(-.25,.65)
)
abline(h = c(.3),lty=2)

#5k-now
sub_allpol_df <- allpol_df[which(allpol_df$cluster==1 & allpol_df$age<5000 & allpol_df$age>0),]
clust1 <- (apply(sub_allpol_df[,1:22],2,FUN = function(x) cor(x, sub_allpol_df$biomass)))

sub_allpol_df <- allpol_df[which(allpol_df$cluster==2 & allpol_df$age<5000 & allpol_df$age>0),]
clust2 <- (apply(sub_allpol_df[,1:22],2,FUN = function(x) cor(x, sub_allpol_df$biomass)))

sub_allpol_df <- allpol_df[which(allpol_df$cluster==3 & allpol_df$age<5000 & allpol_df$age>0),]
clust3 <- (apply(sub_allpol_df[,1:22],2,FUN = function(x) cor(x, sub_allpol_df$biomass)))

all_clust <- rbind(clust1, clust2, clust3)
barplot(
  clust1,
  las = 2,
  col = 1,
  beside = T,
  cex.names = .75,
  main = 'Late Holocene 5ka - Modern',
  ylim = c(-.25,.65)
)
abline(h = c(.3),lty=2)
barplot(
  clust2,
  las = 2,
  col = 2,
  beside = T,
  cex.names = .75,
  ylim = c(-.25,.65)
)
abline(h = c(.3),lty=2)
barplot(
  clust3,
  las = 2,
  col = 3,
  beside = T,
  cex.names = .75,
  ylim = c(-.25,.65)
)
abline(h = c(.3),lty=2)
dev.off()

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

#####
##### calculate correalations for different pollen groups
#####

pollen_sum <- rowSums(pollen_all[,-c(1,24:26)])


#### arboreal and nonarboreal
oftrees_df <- data.frame(arboreal = rowSums(pollen_all[,colnames(pollen_all)[-c(1,3,6,7,24:26)]])/pollen_sum,
                         nonarboreal = rowSums(pollen_all[,c('prairie','other_herbs','CYPERACE')])/pollen_sum,
                         age = pollen_all$age_bacon,
                         cluster = cluster_all,
                         site = pollen_all$site_index_pol,
                         biomass = rep(NA,length(pollen_all$site_index_pol)))

for(ii in 1:nrow(oftrees_df)){
  get_row <- round(oftrees_df$age[ii]/100)
  if(get_row >= 1 & get_row <= 100){
    oftrees_df$biomass[ii] <- agb.mat[which(agb.mat$site_index %in% oftrees_df$site[ii]), seq(2,101,1)[get_row]]
  }
}

oftrees_df <- oftrees_df[-which(is.na(oftrees_df$biomass)),]

sub_oftrees_df <- oftrees_df[which(oftrees_df$cluster==1),]
cor(sub_oftrees_df$arboreal,sub_oftrees_df$biomass)

#### all pollen types
allpol_df <- data.frame( prop.table(as.matrix(pollen_all[,2:23]),1),
                         age = pollen_all$age_bacon,
                         cluster = cluster_all,
                         site = pollen_all$site.name,
                         biomass = rep(NA,length(pollen_all$site_index_pol)))

### add biomass
for(ii in 1:nrow(allpol_df)){
  get_row <- round(allpol_df$age[ii]/100)
  if(get_row >= 1 & get_row <= 100){
    allpol_df$biomass[ii] <- agb.mat[which(agb.mat$name %in% allpol_df$site[ii]), seq(2,101,1)[get_row]]
  }
}

### remove samples without biomass estimate
allpol_df <- allpol_df[-which(is.na(allpol_df$biomass)),]

### gut checks
plot(allpol_df[which(allpol_df$cluster==3),'age'],
     allpol_df[which(allpol_df$cluster==3),'biomass'])

par(mfrow=c(4,5))
for(ii in 1:22){
  plot(allpol_df[which(allpol_df$cluster==2),'age'],
       allpol_df[which(allpol_df$cluster==2),ii])
}

#### adding in biome type pollen
arboreal <- rowSums(allpol_df[,-c(2,5,6,23:26)])
savanna <- rowSums(allpol_df[,c('QUERCUS','CARYA','ULMUS','POPULUS')])
forest <- rowSums(allpol_df[,c('prairie',"other_herbs",'CYPERACE','QUERCUS','CARYA','ULMUS','POPULUS')])

allpol_df <- cbind(arboreal,savanna,forest,allpol_df)

