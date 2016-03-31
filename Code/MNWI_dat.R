data.dir = c("/Users/paleolab/babySTEPPS/Data/")
fig.dir = c('/Users/paleolab/babySTEPPS/Figures/')

#####
##### Prediction Data -- BIGWOODS 8/13/14 #####
#####

##### Install #####
library(gtools)
library(neotoma)
library(reshape)
library(ggplot2)
require(grid)
require(plyr)
require(maps)
library(gtools)
library(sp)
library(rgdal)
library(fields)
library(maptools)
library(mapplots)
gpclibPermit()



#####
##### Download Data for MN and WI #####
#####
gpids <- get_table(table.name='GeoPoliticalUnits')
gpid <- gpids[which(gpids$GeoPoliticalName == 'Minnesota'),1]
gpid1 <- gpids[which(gpids$GeoPoliticalName == 'Wisconsin'),1]
gpid2 <- gpids[which(gpids$GeoPoliticalName == 'Michigan'),1]
meta <- get_dataset(datasettype='pollen', gpid=c(gpid), ageyoung=0) #Pollen data for all of MN&WI
meta1 <- get_dataset(datasettype='pollen', gpid=c(gpid1), ageyoung=0) #Pollen data for all of MN&WI
meta2 <- get_dataset(datasettype='pollen', gpid=c(gpid2), ageyoung=0) #Pollen data for all of MN&WI

loc.PI = rep(0,length(c(meta)))
for(i in 1:length(loc.PI)) loc.PI[i] = as.character(meta[[i]]$pi.data$ContactName)

hk_counts = read.csv(paste0(data.dir,"hotchkiss_lynch_calcote_counts_v0.csv"))
hk_meta = read.csv(paste0(data.dir,"hotchkiss_lynch_calcote_meta_v0.csv"))

site.locs <- ldply(meta, function(x) c(x$site.data$long, x$site.data$lat))
site.locs1 <- ldply(meta1, function(x) c(x$site.data$long, x$site.data$lat))
site.locs2 <- ldply(meta2, function(x) c(x$site.data$long, x$site.data$lat))
site.locs<-rbind(site.locs,site.locs1,site.locs2)

pdf(paste0(fig.dir,"all.sites.neotoma.hotchkiss.pdf"))
map('state', xlim=range(site.locs$V1)+c(-2, 2), ylim=range(site.locs$V2)+c(-1, 1))
points(site.locs$V1, site.locs$V2, pch=19, cex=1,col="black")
points(hk_meta$long, hk_meta$lat, pch=19, cex=1,col="black")
title("All Pollen Sites")
dev.off()

mnwi1 <- rep(0,length(c(meta)))
mnwi2 <- rep(0,length(c(meta1)))
mnwi3 <- rep(0,length(c(meta2)))
for(i in 1:length(meta)) mnwi1[i] <- meta[[i]]$dataset.meta$dataset.id
for(i in 1:length(meta1)) mnwi2[i] <- meta1[[i]]$dataset.meta$dataset.id
for(i in 1:length(meta2)) mnwi3[i] <- meta2[[i]]$dataset.meta$dataset.id
datasets <- as.vector(c(mnwi1,mnwi2,mnwi3))
dat.mnwi <- lapply(datasets, function(x)try(get_download(x)))#get_download(x = as.vector(c(mnwi1,mnwi2,mnwi3)))
save(dat.mnwi,file="mnwi1.rdata")
load(file="mnwi.rdata") # start here unless you think there might be new data in neotoma

##### Get only metadata
pol.cal.data=matrix(0,length(c(meta,meta1,meta2)),5)
for(i in 1:length(c(meta,meta1,meta2))){	
  pol.cal.data[i,1] <- dat.mnwi[[i]][[1]]$dataset$site.data$site.id
  pol.cal.data[i,2] <- dat.mnwi[[i]][[1]]$dataset$site.data$lat
  pol.cal.data[i,3] <- dat.mnwi[[i]][[1]]$dataset$site.data$long
  pol.cal.data[i,4] <- dat.mnwi[[i]][[1]]$dataset$site.data$site.name
}

pol.cal.data<-as.data.frame(pol.cal.data)

##### prepare data frame for adding counts
get.count <- function(x,i){	
	compile_taxa(x[[i]],'WhitmoreSmall')
}

x=dat.mnwi

pol.cal.count=data.frame(cbind(rep(pol.cal.data[1,1],
              length(get.count(dat.mnwi,1)$sample.meta$Age)),
              rep(pol.cal.data[1,2],length(get.count(dat.mnwi
              ,1)$sample.meta$Age)),rep(pol.cal.data[1,3],
              length(get.count(dat.mnwi,1)$sample.meta$Age)),
              rep(pol.cal.data[1,4],length(get.count(dat.mnwi,
              1)$sample.meta$Age)),rep(pol.cal.data[1,5],
              length(get.count(dat.mnwi,1)$sample.meta$Age)),
              get.count(dat.mnwi,1)$sample.meta$Age,
              as.matrix(get.count(dat.mnwi,1)$count)))
head(pol.cal.count)

##### add counts to data frame using function "smartbind" because 
##### rows are of varying lengths i.e. some have more types of pollen than others
##### expect it to take a while and several warnings

for(i in 2:nrow(pol.cal.data)){
		save1<-data.frame(cbind(rep(pol.cal.data[i,1],
    length(get.count(dat.mnwi,i)$sample.meta$Age)),
    rep(pol.cal.data[i,2],length(get.count(dat.mnwi,
    i)$sample.meta$Age)),rep(pol.cal.data[i,3],
    length(get.count(dat.mnwi,i)$sample.meta$Age)),
    rep(pol.cal.data[i,4],length(get.count(dat.mnwi,
    i)$sample.meta$Age)),rep(pol.cal.data[i,5],
    length(get.count(dat.mnwi,i)$sample.meta$Age)),
    get.count(dat.mnwi,i)$sample.meta$Age,
    as.matrix(get.count(dat.mnwi,i)$count)))
	pol.cal.count=smartbind(pol.cal.count,save1)	
	}	

##### Fix problems with matrix
rownames(pol.cal.count)<-seq(1,nrow(pol.cal.count),1)
pol.cal.count[is.na(pol.cal.count)]<-0
colnames(pol.cal.count)<-c("SiteID","LatitudeNorth","LongitudeWest","dataset.id","ContactName","Age",colnames(pol.cal.count[,7:ncol(pol.cal.count)]))
#save(pol.cal.count,file="pol.cal.count.mnwi1.csv")

load(paste0(data.dir,"pol.cal.count.mnwi1.csv"))
load("hk_counts3.csv")

to_pol_mat = which(colnames(hk_counts3[,18:ncol(hk_counts3)])%in%colnames(pol.cal.count))

colnames(hk_counts3)
colnames(pol.cal.count)

add_on = matrix(0,ncol=(ncol(hk_counts3) + 6 - 17),nrow=nrow(hk_counts3))
add_on[,1] = hk_counts3$name
add_on[,2] = hk_counts3$lat
add_on[,3] = hk_counts3$long
add_on[,6] = hk_counts3$AGE
add_on[,7:ncol(add_on)] = as.matrix(hk_counts3[,18:ncol(hk_counts3)])

colnames(add_on) = c(colnames(pol.cal.count[,1:6]),colnames(hk_counts3[,18:ncol(hk_counts3)]))

pol.cal.count = smartbind(pol.cal.count,add_on)

##### Fix problems with matrix
rownames(pol.cal.count)<-seq(1,nrow(pol.cal.count),1)
pol.cal.count[is.na(pol.cal.count)]<-0
colnames(pol.cal.count)<-c("SiteID","LatitudeNorth","LongitudeWest","dataset.id","ContactName","Age",colnames(pol.cal.count[,7:ncol(pol.cal.count)]))
#save(pol.cal.count,file="pol.cal.count.mnwi1.csv")

head(pol.cal.count)
ponds = pol.cal.count[pol.cal.count$Age>0,]
#ponds = ponds[ponds$Age<=5000,]
ponds = ponds[order(ponds$LatitudeNorth),]
site.factor = factor(ponds[,1],labels = seq(1,length(unique(ponds[,1])),1))
ponds1 = cbind(site.factor,ponds)

ponds_rm = which(ponds1$LatitudeNorth<44.5&ponds1$LongitudeWest>c(-86))
ponds1 = ponds1[-ponds_rm,]

pdf(paste0(fig.dir,"chrono.10k.pdf"))
par(mfrow=c(1,1))
#plot(ponds1$Age,ponds1$LatitudeNorth,pch=19,cex=.25,main="All Records",ylab="Latitude",xlab="Age BP")
#abline(v=2000)
plot(ponds1$Age,ponds1$LatitudeNorth,pch=19,cex=.5,xlim=c(0,10000),main="All Pollen Records",ylab="Latitude",xlab="Age BP")
dev.off()

plot.seq = seq(0,20000,500)

#quartz()
par(mfrow = c(2,2))
for(i in 2:length(plot.seq)){
  map('state', xlim=range(ponds1$LongitudeWest)+c(-2, 2), ylim=range(ponds1$LatitudeNorth)+c(-1, 1))
  points(ponds1[ponds1$Age>plot.seq[i-1]&ponds1$Age<plot.seq[i],]$LongitudeWest, 
         ponds1[ponds1$Age>plot.seq[i-1]&ponds1$Age<plot.seq[i],]$LatitudeNorth, 
         pch=19, cex=.5)
  title(c(plot.seq[i-1],"-",plot.seq[i]))
}
dev.off()


