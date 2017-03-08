#Comes after future.prediction.R
load("~/Downloads/future.samples.pred.Rdata")
samples.pred1 <- samples.pred
load("~/Downloads/future.samples.pred1.Rdata")
samples.pred1<-cbind(samples.pred1,samples.pred)
load("~/Downloads/future.samples.pred2.Rdata")
samples.pred1<-cbind(samples.pred1,samples.pred)
load("~/Downloads/future.samples.pred3.Rdata")
samples.pred<-cbind(samples.pred1,samples.pred)

all.samps.dat <- cbind(x.meta,colMeans(samples.pred))
hist(samples.pred,main="Age All Samples",col=colors,breaks=breaks)


par(mfrow=c(5,5))
for(i in 1:length(unique(all.samps.dat$SiteID))){
	plot(all.samps.dat[all.samps.dat[,1]==unique(all.samps.dat$SiteID)[i],7],main=unique(all.samps.dat$SiteID)[i])
 }


head(all.samps.dat)
colnames(all.samps.dat)<-c("SiteID","y","x","dat.id","contact",
                        "age","mean")
                        
new.site.locs <- cbind(all.samps.dat$x,all.samps.dat$y)
centers_pol = data.frame(new.site.locs)
colnames(centers_pol) = c('x', 'y')

coordinates(centers_pol) <- ~ x + y
proj4string(centers_pol) <- CRS('+proj=longlat +ellps=WGS84')

centers_polA <- spTransform(centers_pol, CRS('+init=epsg:3175'))
centers_polA <- as.matrix(data.frame(centers_polA))

all.samps.dat$x<-centers_polA[,1] 
all.samps.dat$y<-centers_polA[,2]                       

b <- gam(log(mean) ~ te(y,x,age, d = c(2,1),bs = c("tp","cr"),k=c(25,10)), data = as.data.frame(all.samps.dat))
summary(b)
vis.gam(b)



biomass_dat_est <- read.csv("~/babySTEPPS/Data/biomass_prediction_v0.2.csv")
d.save = list()
age.seq = seq(500,3000,500)
breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)-1))

par(mfrow=c(2,3))
for(i in 1:length(age.seq)){
pts.inc <- all.samps.dat[all.samps.dat$age>(age.seq[i]-100)&all.samps.dat$age<(age.seq[i]+100),]
data_binned1 <-  cut(pts.inc$mean, breaks, include.lowest = TRUE, labels = FALSE)

plot(pts.inc$x,pts.inc$y,ylim=range(all.samps.dat$y)+c(-100000,100000),xlim=range(all.samps.dat$x)+c(-100000,100000),pch=19,cex=2,col= colors[data_binned1],xlab=NA,ylab=NA,main=age.seq[i])
plot(usShp,add=TRUE)
}


for(i in 1:length(age.seq)){
	d.save[[i]]<-plot_age(age.seq[i])
}
grid.arrange(d.save[[1]],d.save[[2]],d.save[[3]],d.save[[4]],d.save[[5]],d.save[[6]],nrow=3)

usShp <- readShapeLines(file.path("/Users/paleolab/Documents/babySTEPPS/", 'us_alb.shp'), proj4string=CRS('+init=epsg:3175'))
usShp@data$id <- rownames(usShp@data)
usFortified <- fortify(usShp, region='id')

##### add state lines function
add_map_albers <- function(plot_obj, map_data = usShp, dat){
  p <- plot_obj + geom_path(data = map_data, aes(x = long, y = lat, group = group), size = 1) + 
    scale_x_continuous(limits = c(min(dat$x, na.rm = TRUE), max(dat$x, na.rm = TRUE))) +
    scale_y_continuous(limits = c(min(dat$y, na.rm = TRUE), max(dat$y, na.rm = TRUE)))
  return(p)
}


##### make a heat plot of biomass
theme_clean <- function(plot_obj){
  plot_obj <- plot_obj + theme(axis.ticks = element_blank(),
                               axis.text.y = element_blank(),
                               axis.text.x = element_blank(),
                               axis.title.x = element_blank(),
                               axis.title.y = element_blank())
  
  return(plot_obj)
}

plot_age <- function(age){
newdata = as.data.frame(cbind(biomass_dat_est[,1:2],rep(age,nrow(biomass_dat_est[,1:2]))))
colnames(newdata)<-c('x','y','age')

pred_biomass_gam = exp(predict(b,newdata = newdata,type="response"))


for(n in 1:length(pred_biomass_gam)){
	if(as.vector(pred_biomass_gam[n])>200){
		pred_biomass_gam[n]=200
	} 
}

#hist(pred_biomass_gam)

full.mat <- cbind(biomass_dat_est[,1:2], pred_biomass_gam)
colnames(full.mat) <- c("x","y","Predicted Biomass")
y = as.data.frame(full.mat)

#Regular breaks and colors
breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)-1))

legendName <- "Biomass (Mg/ha)"

data_binned <-  cut(y[,3], breaks, include.lowest = TRUE, labels = FALSE)

breaklabels <- apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1,  function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })

inputData <- data.frame(X = y[,1], Y = y[,2], Preds = cbind(data_binned,data_binned))
inputData_long <- melt(inputData, c('X', 'Y'))

input_points <- data.frame(all.samps.dat$x,all.samps.dat$y)
colnames(input_points) <- c('lon','lat')
input_points<- input_points[all.samps.dat$age>(age-100)&all.samps.dat$age<(age+100),]

d <- ggplot() + geom_raster(data = inputData_long, aes(x = X, y = Y, fill = factor(value))) + scale_fill_manual(labels = breaklabels, name = legendName, drop = FALSE, values = colors, guide = FALSE) + 
  theme(strip.text.x = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16)) +
  geom_point(data = input_points, aes(x=lon,y=lat), pch=16, size=1,colour="black") +
  ggtitle(paste("Biomass",age))

add_map_albers <- function(plot_obj, map_data = usFortified, dat){
  p <- plot_obj + geom_path(data = map_data, aes(x = long, y = lat, group = group), size = 0.1) +
    scale_x_continuous(limits = c(min(inputData$X, na.rm = TRUE), max(inputData$X, na.rm = TRUE))) +
    scale_y_continuous(limits = c(min(inputData$Y, na.rm = TRUE), max(inputData$Y, na.rm = TRUE)))
  return(p)
}

d <- add_map_albers(plot_obj = d, map_data = usFortified, dat = inputData_long)

return(d=d)  
}

