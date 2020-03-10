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

plot(unlist(long),unlist(lat),col=clusters@cluster,pch=19,cex=2)
maps::map('state',add=T)
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
  if(any(which(pollen_all$site.name[ii]==names(agb.list)))){
    cluster_all[ii] <- clusters@cluster[which(pollen_all$site.name[ii]==names(agb.list))]
  }   
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
plot(allpol_df[which(allpol_df$cluster==2),'age'],
     allpol_df[which(allpol_df$cluster==2),'biomass'])

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

#### sensitivity analysis

now_all <- do.call(rbind,info) #info comes from the average.biomass figure

#### adding in biome type pollen
arboreal <- rowSums(now_all[,-c(2,5,6,23:27)])
savanna <- rowSums(now_all[,c('QUERCUS','CARYA','ULMUS','POPULUS')])
forest <- 1-rowSums(now_all[,c('prairie',"other_herbs",'CYPERACE','QUERCUS','CARYA','ULMUS','POPULUS')])

now_all_all <- cbind(arboreal,savanna,forest,now_all)

clusters_assign <- clusters@cluster[as.character(now_all$site)]
now_all_df <- data.frame(now_all_all,clusters = clusters_assign)

#gut checks
plot(now_all_df[now_all_df$clusters==3,'lon'],now_all_df[now_all_df$clusters==3,'lat'])
maps::map('usa',add=T)
plot(now_all_df[now_all_df$clusters==1,'age'],now_all_df[now_all_df$clusters==1,'max_est'],xlim=c(100,0))
plot(now_all_df[now_all_df$site==1,'age'],now_all_df[now_all_df$site==1,'max_est'],xlim=c(100,0))
plot(now_all_df[now_all_df$site==1,'age'],now_all_df[now_all_df$site==1,'TSUGAX'],xlim=c(100,0))


#sensitivity calculating
now_all_df$site <- as.numeric(as.factor(now_all_df$site))
pollen_save <- biomass_save <- cluster_save <- site_dat_save <- site_lat_save <- site_lon_save <- list()
ages_set <- matrix(c(100,80,81,50,51,0),3,2,byrow = T)
for(ii in 1:77){
  site_dat <- now_all_df[now_all_df$site==ii,]
  pollen_save[[ii]] <- biomass_save[[ii]] <- cluster_save[[ii]] <- site_dat_save[[ii]] <- site_lat_save[[ii]] <- site_lon_save[[ii]] <- list()
  for(pp in 1:25){
    pollen_save[[ii]][[pp]] <- biomass_save[[ii]][[pp]] <- cluster_save[[ii]][[pp]] <-site_dat_save[[ii]][[pp]] <- site_lat_save[[ii]][[pp]] <- site_lon_save[[ii]][[pp]] <- list()
    ages <- site_dat[,'age']
    for(tt in 1:3){
      row_get <- which(ages <= ages_set[tt,1] &ages>=ages_set[tt,2])
      
      row_order <- order(ages[row_get],decreasing = T)
      
      biomasses <- site_dat[row_get,'max_est']
      pollens <- site_dat[row_get,pp]
      dist_mat <- outer(ages[row_get],ages[row_get],FUN='-')
      dist_mat_biomass <- outer(biomasses,biomasses,FUN='-')
      dist_mat_pollen <- outer(pollens,pollens,FUN='-')
    
      #dist_mat[lower.tri(dist_mat)] <- dist_mat_biomass[lower.tri(dist_mat_biomass)] <- dist_mat_pollen[lower.tri(dist_mat_pollen)] <- NA
      
      dist_mat_pollen[abs(dist_mat)>10] <- dist_mat_biomass[abs(dist_mat)>10] <- NA
     
      remove <- which(duplicated(c(dist_mat_pollen)))
      
      pollen_save[[ii]][[pp]][[tt]] <- dist_mat_pollen#[-remove]
      biomass_save[[ii]][[pp]][[tt]] <- dist_mat_biomass#[-remove]
      cluster_save[[ii]][[pp]][[tt]] <- rep(site_dat$cluster[1],length(dist_mat_biomass))
      site_dat_save[[ii]][[pp]][[tt]] <- rep(site_dat$site[1],length(dist_mat_biomass))
      site_lat_save[[ii]][[pp]][[tt]] <- rep(site_dat$lat[1],length(dist_mat_biomass))
      site_lon_save[[ii]][[pp]][[tt]] <- rep(site_dat$lon[1],length(dist_mat_biomass))
      #plot(dist_mat_pollen[-remove],dist_mat_biomass[-remove])
    }
  }
}
time_names <- c('early','mid','late')
colors <- c('black','darkred','darkgreen')

#### Overall 
pdf('pol_biom_diff_scatters_max_liks_time_periods_same_scale.pdf',height = 20,width=20)
par(mfrow=c(1,1))
plot(now_all_df$lon,now_all_df$lat,col=now_all_df$clusters,pch=19,cex=3)
maps::map('state',add=T)
par(mfrow=c(3,3),mar = c(3,3,2,0),oma = rep(4,4))
layout(matrix(1:9,3,3))
  for(ss in 1:25){
    for(tt in 1:3){
    all_pollen <- do.call(c, lapply(
      pollen_save,
      FUN = function(x) {
        x[[ss]][[tt]]
      }
    ))
    all_biomass <- do.call(c, lapply(
      biomass_save,
      FUN = function(x) {
        x[[ss]][[tt]]
      }
    ))
    
    all_cluster <- do.call(c, lapply(
      cluster_save,
      FUN = function(x) {
        x[[ss]][[tt]]
      }
    ))
    for(cc in 1:3){
      Lab.palette <- colorRampPalette(c('white',colors[cc],'black'), space = "Lab")
      subset <- which(all_cluster==cc)#1:length(all_pollen)#runif(1000,1,length(all_pollen))
      smoothScatter((all_pollen[subset]),
                    (all_biomass[subset]),
                    main = ,cex.axis=3,
                    ylim = c(-250,250),
                    xlim = c(-.5,.5),
                    xlab = 'Pollen Proportion Difference',
                    ylab = 'Biomass Difference',colramp = Lab.palette)
      title(paste(time_names[tt],colnames(now_all_df)[ss]),cex.main=3)
      abline(h=0,col='blue')
    }
  }
}

mtext('Pollen Proportion Difference',outer = T,side = 1, line = 2)
mtext('Biomass Estimate Difference',outer = T,side = 2, line = 2)
dev.off()

#### By Site
lm_obj <- array(NA,dim = c(77,25,3))

pdf('scatter_by_site.pdf')
for(ii in 1:77){
  
  par(mfrow=c(1,1))
  plot(now_all_df$lon,now_all_df$lat,col=now_all_df$clusters,pch=1,cex=3)
  maps::map('state',add=T)
  text(now_all_df[now_all_df$site==ii,'lon'][1],now_all_df[now_all_df$site==ii,'lat'][1],pch=1,labels = now_all_df[now_all_df$site==ii,'site'][1])
  
  plot(now_all_df[now_all_df$site==ii,'age'],now_all_df[now_all_df$site==ii,'max_est'],
       xlab = 'Age',
       ylab = 'Max Lik Biomass', xlim=c(100,0),ylim=c(0,250))
  
  par(mfrow=c(3,3))
  for(pp in 1:25){
    for(tt in 1:3){
      plot(c(pollen_save[[ii]][[pp]][[tt]]),c(biomass_save[[ii]][[pp]][[tt]]),
           col=cluster_save[[ii]][[pp]][[tt]][1], pch = 1,
           xlab = 'Pollen Proportion Difference',ylab='Biomass Difference')
      if(max(c(pollen_save[[ii]][[pp]][[tt]]),na.rm=T)>.03 & length(!is.na(pollen_save[[ii]][[pp]][[tt]]))>10) {
        keep_lm <- rlm(c(biomass_save[[ii]][[pp]][[tt]])~c(pollen_save[[ii]][[pp]][[tt]])+0)
        sum_lm <- summary(keep_lm)
        lm_obj[ii,pp,tt] <- keep_lm$coefficients
      }
      if(!is.na(lm_obj[ii,pp,tt])) abline(a=0,b= lm_obj[ii,pp,tt],col='blue')
      title(paste(ii,time_names[tt],colnames(now_all_df)[pp]),cex.main=1)
    }
  }
}
dev.off()

pdf('aggergated_slopes.pdf')
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
clust_get <- unique(now_all_df[now_all_df$clusters==cc,'site'])
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
