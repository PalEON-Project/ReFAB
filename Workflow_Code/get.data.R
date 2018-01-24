library(neotoma)
library(maptools)
library(ggplot2)
library(reshape)
library(fields)
library(oce)
library(splines)
library(mgcv)

data.dir = c("Data/")
fig.dir = c('Figures/')
DRAW = TRUE

usShp <- readShapeLines(file.path(data.dir, '/us_alb.shp'), proj4string=CRS('+init=epsg:3175'))
usShp@data$id <- rownames(usShp@data)
usFortified <- fortify(usShp, region='id')

#####
##### Download Data for MN and WI and MI #####
#####

meta <- get_dataset(datasettype='pollen', gpid=c("Wisconsin", 
                                                 "Michigan",
                                                 "Minnesota"), ageyoung=0)
meta_dl <- get_download(meta)

comp.tax <- compile_taxa(meta_dl, 'WhitmoreSmall')
pol_cal_count <- compile_downloads(comp.tax)

save(pol_cal_count,file=paste0('nimble_pull',Sys.Date(),'.Rdata'))

#####
##### Add Bacon Pollen Dates #####
#####

bacon <- read.csv(paste0(data.dir,'/sediment_ages_v1.0_varves.csv'))
cal <- read.csv(paste0(data.dir,'/cal_data_mid_depth_2015-06-10.csv'))

cal.new.pol <- merge(x = pol_cal_count, y = cal, 
      by.x = c('dataset','depth'), 
      by.y = c('id','depth'))
length(unique(cal.new.pol$dataset))

# site.find <- cal[-na.omit(match(cal$id,cal.new.pol$dataset)),'site']
# 
# d.id <- list()
# for(i in 1:length(site.find)){
#   cal[cal$site==site.find[i],'id']<-pol_cal_count[grep(site.find[i],pol_cal_count$site.name),'dataset'][1]
# }


#bacon <- read.csv(paste0(data.dir,'/pollen_ts_bacon_v8.csv'))
new.pollen <- merge(x = pol_cal_count, y = bacon, 
                    by.x = c('.id','age','depth','PINUSX','lat','long'), 
                    by.y = c('id','age_default','depth','PINE','lat','long'))
pol_list <- list()
pol_ids <- unique(pol_cal_count$dataset)
for(i in 1:length(pol_ids)){
  pol_list[[i]] <- pol_cal_count[pol_cal_count$dataset==pol_ids[i],]
}
names(pol_list) <- pol_ids


### Test
sum(new.pollen$BIRCH - new.pollen$BETULA) #should equal zero

#####
###### Format settlement pollen data #####
#####

#original
x <- new.pollen[new.pollen$age_bacon >= 100, ]
x <- x[x$age_bacon<=250,]

#just settlement horizon
colnames(cal.new.pol)[9] <- c('lat')
colnames(cal.new.pol)[10] <- c('long')
x <- cal.new.pol

#settlement horizon + previous 100 years
pol_settle <- list()
for(i in 1:length(pol_list)){
  #cal_list <- merge(x = pol_list[[i]], y = cal, 
   #                 by.x = c('dataset','depth'), 
    #                by.y = c('id','depth'))
  cal_stop <- cal[which(cal$id==names(pol_list)[i]),]
  cal_list <- pol_list[[i]][which(pol_list[[i]]$depth==cal_stop$depth),]
  # this is not using bacon ages I think... Ask Andria if cal has bacon ages. 
  #We need the bacon ages to aggregate 100 years back from settlement horizon.
  i_cal <- which(pol_list[[i]]$age >= cal_list$age & pol_list[[i]]$age < (cal_list$age + 100))
  pol_settle[[i]] <- pol_list[[i]][i_cal,]
  
}
## might not be getting all of the sites
x <- pol_settle_mat <- do.call(rbind,pol_settle)
length(unique(x$dataset))

all.pollen.taxa.names <- colnames(pol_cal_count)[11:length(colnames(pol_cal_count))]
melt.x <- melt(x, id.vars=c('.id','age','lat','long','site.name'),
               measure.vars = all.pollen.taxa.names)
cast.x <- cast(melt.x, lat + long + .id + site.name ~ variable, sum) #summing over pollen taxa for years around settlement
cast.x <- cast.x.full <- as.data.frame(cast.x)


#####
##### Remove sites and take subset of points
#####

if(DRAW == TRUE) pdf(paste0(fig.dir,paste0("all_sites",Sys.Date(),".pdf")))
map('state', ylim=range(cast.x$lat)+c(-2, 2), xlim=range(cast.x$long)+c(-1, 1),main=NA)
points(cast.x$long, cast.x$lat, pch=19, cex=1,col="gray")
title(main="all sites")

set.seed(4)
sites_rm = sample(1:nrow(cast.x),round(nrow(cast.x)/3))

points(cast.x$long[-sites_rm], cast.x$lat[-sites_rm], pch=19, cex=1,col="blue")
if(DRAW == TRUE) dev.off()

#plot_biomass_pollen = plot_biomass_pollen[-c(which(props$Other>.5),which(props$POACEAE>.8)),]

#####
##### Adding settlement biomass data #####
#####

biomass_dat_est <- read.csv(paste0(data.dir,"biomass_prediction_v0.9-10_bam.csv"))


##### Changing pollen coordinates so that we can find the right rows when we find the biomass for each pond.
lat.long.reg <- cbind(as.numeric(as.character(cast.x$long)),as.numeric(as.character(cast.x$lat)))
lat.long.reg.df = data.frame(lat.long.reg)
colnames(lat.long.reg.df) = c('x', 'y')

coordinates(lat.long.reg.df) <- ~ x + y
proj4string(lat.long.reg.df) <- CRS('+proj=longlat +ellps=WGS84')

albers <- spTransform(lat.long.reg.df, CRS('+init=epsg:3175'))
albers <- as.matrix(data.frame(albers))

save(albers,file=paste0(Sys.Date(),'calibration.albers.Rdata'))

centers_biomass = cbind(biomass_dat_est$x,biomass_dat_est$y)
idx_cores = vector(length=nrow(cast.x))

for(i in 1:nrow(cast.x)){   
  core_site = albers[i,]
  d = rdist(matrix(core_site, ncol=2), as.matrix(centers_biomass))
  if(min(d)<8000){
  	idx_cores[i] = which.min(d) 
  }else{
  	idx_cores[i] = NA 
  }
  
}

### Map checking for biomass grid cell that pollen core is contained within
if(DRAW == TRUE) pdf(paste0(fig.dir,"check_points.pdf")) else (quartz())
plot(albers[,1], albers[,2])
points(centers_biomass[idx_cores,1], centers_biomass[idx_cores,2], col='blue', pch=8)
plot(usShp, add=T, lwd=2)
if(DRAW == TRUE) dev.off()

biomass <- list()
for(i in 1:length(idx_cores)){ 
  biomass[[i]] = biomass_dat_est[idx_cores[i],'Total']
 }
 
cast.x <- cbind(cast.x,unlist(biomass))

#### Remove places where we do not have settlement biomass
cast.x <- cast.x[-which(is.na(cast.x[,ncol(cast.x)])),]
biomass <- biomass[-which(is.na(biomass))]

biomass <- unlist(biomass)

breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)-1))
data_binned <-  cut(cast.x[,ncol(cast.x)], c(breaks), include.lowest = FALSE, labels = FALSE)

pdf(paste0(fig.dir,paste0("biomass.pts.settlement",Sys.Date(),".pdf")))
map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(cast.x$long,cast.x$lat, pch=21,
		cex=1.1, bg=colors[data_binned],lwd=.2)
title("Biomass Point Estimates at Settlement")
plotInset(-90,47,-82.5,50,
          expr={
          	hist(data_binned,col=colors,xaxt="n",xlab=NA,
          	ylab=NA,main=NA,cex.lab=.5,cex.axis=.5)
          	axis(side=1,breaks,at = seq(1,12,1),cex.axis = .5,las=2,line=0)
          	mtext(side = 1, "Biomass (Mg/ha)", line = 1.5,cex=.5)
         mtext(side = 2, "Frequency", line = 1.7,cex=.5)
          	})
dev.off()

#####
##### Creating a calibration dataset with the species we want to use
#####

### I do this after indexing because some of the places where we have settlement pollen we don't have settlement biomass.

trees <- c('PINUSX',"ALNUSX","JUGLANSX","ACERX","CUPRESSA","FRAXINUX","FAGUS","LARIXPSEU","TSUGAX","QUERCUS","TILIA",
           "BETULA","PICEAX","OSTRYCAR","ULMUS","ABIES","POPULUS","CARYA","CYPERACE")
other.trees <- c("TAXUS","NYSSA","CASTANEA","PLATANUS","SALIX","LIQUIDAM")
ten.count = matrix(0,nrow(cast.x),length(trees)+3)
prairie <- c("ARTEMISIA","ASTERX","POACEAE","AMBROSIA","CHENOAMX","CORYLUS")
ten.count[,1] <- unlist(rowSums(cast.x[,prairie],na.rm = T))
ten.count[,2] <- unlist(rowSums(cast.x[,other.trees],na.rm = T))
ten.count[,3:(length(trees)+2)] <- as.matrix(cast.x[,trees])
ten.count[,(length(trees)+3)] <- rowSums(cast.x[,all.pollen.taxa.names],na.rm = T) - rowSums(ten.count,na.rm = T)
colnames(ten.count)<-c("prairie","other_trees",trees,"other_herbs")

total_counts = round(rowSums(ten.count,na.rm = TRUE))
counts = round(ten.count)
props = counts/rowSums(counts,na.rm = TRUE)

total_counts_spp = colSums(counts)

props = props[,order(total_counts_spp,decreasing=TRUE)]

if(DRAW==TRUE) pdf(paste0(fig.dir,"scatter.newdata",Sys.Date(),".pdf"))
#quartz()
par(mfrow=c(4,4))
for(i in 1:ncol(props)){
  if(length(unique(props[-sites_rm,i]))>=9){
    plot(biomass[-sites_rm],props[-sites_rm,i],main=colnames(props)[i],xlab="biomass",ylab="pollen prop",pch = 19, cex = .5)
  }
}  
if(DRAW==TRUE) dev.off()

props = as.data.frame(props)

#####
##### Drawing Splines (with all data) #####
#####

Z = bs(unlist(biomass),intercept=TRUE,df=5)
betas = matrix(0,ncol(Z),ncol(counts)); betas.save = betas


counts[is.na(counts)] <- 0 ### important. Could put above in merge or ten.count creation.

if(DRAW == TRUE) pdf(paste0(fig.dir,"splines.new.pdf"))
#quartz()
par(mfrow=c(3,3))
for(i in 1:ncol(counts)){
  gam_mod = gam(cbind(counts[,i],total_counts-counts[,i]) ~ s(biomass),family=binomial(link="logit"))
  plot(biomass,counts[,i]/total_counts,pch=19,cex=.4,col='grey',ylab="Pollen Prop",main=colnames(counts)[i])
  points(biomass,predict(gam_mod,type="response"),pch=19,col="green")
  
  glm_mod = glm(cbind(counts[,i],total_counts-counts[,i]) ~ Z - 1,family=binomial(link="logit"))   
  points(biomass,counts[,i]/total_counts,pch=19,cex=.4,col='grey')
  new.biomass = seq(1,200,1)
  Z.new = bs(new.biomass,intercept=TRUE,df = ncol(Z))
  lines(new.biomass, predict(glm_mod,newdata=list(Z=Z.new),type="response"),col="blue")  
  
  betas[,i] = glm_mod$coefficients
  betas.save[,i] = glm_mod$coefficients
  
}

if(DRAW == TRUE) dev.off()

Z.knots<- Z

u<-c(rep(attr(Z.knots,"Boundary.knots")[1],1),attr(Z.knots,"knots"),rep(attr(Z.knots,"Boundary.knots")[2],1))

library(nimble)
source('Workflow_Code/utils/bs_nimble.R')

Z.knots.check = matrix(0,nrow=u[length(u)],ncol=(length(u)+2));

for(i in 1:u[length(u)]){
	u_given <-i
	Z.knots.check[i,] = bs_nimble(u_given, u=u, N0 = rep(0, (length(u)-1)),N1 = rep(0, (length(u))), N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))
}

for(i in 1:length(biomass)){
    u_given <- biomass[i]
	Z.knots[i,] = bs_nimble(u_given, u=u, N0 = rep(0, (length(u)-1)),N1 = rep(0, (length(u))), N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))
}

#### Plot basis functions ####
if(DRAW==TRUE) pdf(file.path(fig.dir,'basis_function_check.pdf'))
plot(Z.knots.check[,1],xlim=c(0,u[length(u)]),pch=19,ylim=c(0,1),xlab="Biomass")
for(i in 2:ncol(Z.knots.check)){
  points(Z.knots.check[,i],col=i,pch=19)
}
abline(v=u,lwd=2)
title("Basis Functions")
if(DRAW==TRUE) dev.off()

#####
##### Create final calibration datasets #####
#####

Y = counts[-sites_rm,] #remove sites "sites_rm" defined above
Y <- Y[,rev(order(colMeans(Y)))]
biomass = biomass[-sites_rm]
counts = counts[-sites_rm,]
counts <- Y[,rev(order(colMeans(Y)))]
total_counts = rowSums(counts)

dim(Y)[1]-length(biomass) # should be zero

save(Y,biomass,file=paste0(Sys.Date(),'calibration.data.Rdata'))

#save.image(file="add.bacon2.Rdata")

u<-c(rep(attr(Z,"Boundary.knots")[1],1),attr(Z,"knots"),rep(attr(Z,"Boundary.knots")[2],1))

#####
##### Create final prediction datasets #####
#####

x = new.pollen[new.pollen$age_bacon>=200,]
x = x[x$age_bacon<=10000,]

x.meta = x[,c('.id','lat',"long","dataset","site.name","age_bacon")]
colnames(x.meta)[1] <- c('site.id')

ten.count = matrix(0,nrow(x),length(trees)+3)
ten.count[,1] <- unlist(rowSums(x[,prairie],na.rm = TRUE))
ten.count[,2] <- unlist(rowSums(x[,other.trees],na.rm = TRUE))
ten.count[,3:(length(trees)+2)] <- as.matrix(x[,trees])
ten.count[,(length(trees)+3)] <- rowSums(x[,all.pollen.taxa.names],na.rm = TRUE) - rowSums(ten.count,na.rm = TRUE)
colnames(ten.count)<-c("prairie","other_trees",trees,"other_herbs")

ten.count.save = ten.count
ten.count = round(ten.count.save)
ten.count <- ten.count[,colnames(counts)]

ten.count[is.na(ten.count)] <- 0

save(x.meta, ten.count, Z, u, file = 'prediction.data.Rdata')

#####
##### Create paleon mip datasets #####
#####

x = new.pollen[new.pollen$age_bacon>=200,]
x = x[x$age_bacon<=2000,]

x.meta = x[,c('.id','lat',"long","dataset","site.name","age_bacon")]
colnames(x.meta)[1] <- c('site.id')

ten.count = matrix(0,nrow(x),length(trees)+3)
ten.count[,1] <- unlist(rowSums(x[,prairie],na.rm = TRUE))
ten.count[,2] <- unlist(rowSums(x[,other.trees],na.rm = TRUE))
ten.count[,3:(length(trees)+2)] <- as.matrix(x[,trees])
ten.count[,(length(trees)+3)] <- rowSums(x[,all.pollen.taxa.names],na.rm = TRUE) - rowSums(ten.count,na.rm = TRUE)
colnames(ten.count)<-c("prairie","other_trees",trees,"other_herbs")

ten.count.save = ten.count
ten.count = round(ten.count.save)
ten.count <- ten.count[,colnames(counts)]


save(u,Z,x.meta,ten.count,file = 'paleon.data.Rdata')


#####
##### Plots #####
#####

if(DRAW==TRUE) pdf(paste0(fig.dir,"all.sites.neotoma",Sys.Date(),".pdf"))
map('state', xlim=range(x.meta$long)+c(-2, 2), ylim=range(x.meta$lat)+c(-1, 1))
points(x.meta$long, x.meta$lat, pch=19, cex=1,col="black")
title("All Pollen Sites")
if(DRAW==TRUE) dev.off()

if(DRAW==TRUE) pdf(paste0(fig.dir,"chrono.10k",Sys.Date(),".pdf"))
par(mfrow=c(1,1))
plot(x.meta$age_bacon,x.meta$lat,pch=19,cex=.5,xlim=c(0,10000),main="All Pollen Records",ylab="Latitude",xlab="Age BP")

plot.seq = seq(0,10000,500)

par(mfrow = c(2,2))
for(i in 2:length(plot.seq)){
  map('state', xlim=range(as.numeric(as.character(x.meta$long)))+c(-2, 2), ylim=range(as.numeric(as.character(x.meta$lat)))+c(-1, 1))
  points(x.meta[x.meta$age>plot.seq[i-1]&x.meta$age<plot.seq[i],]$long, 
         x.meta[x.meta$age>plot.seq[i-1]&x.meta$age<plot.seq[i],]$lat, 
         pch=19, cex=.5)
  title(c(plot.seq[i-1],"-",plot.seq[i]))
}
if(DRAW==TRUE) dev.off()






