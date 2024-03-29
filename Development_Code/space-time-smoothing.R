load("all.preds.min.list.RData")
all.preds = summary(csamp.real.pred)
#follows from Run_Model.R ##all.preds = summary(csamp.real.pred)

site.factor = factor(x[,1],labels = seq(1,142,1))
all.preds1 = cbind(site.factor,x[,1:6],all.preds$quantiles[,1],all.preds$quantiles[,3],all.preds$quantiles[,5])
all.preds1 = cbind(all.preds1,ten.count)
all.preds1 = all.preds1[order(all.preds1$LatitudeNorth),]
all.preds1[,1] = factor(all.preds1[,2],labels = seq(1,142,1),
                        levels=unique(all.preds1[,2]),ordered=FALSE)

#ten.count = ten.count[order(all.preds1$LatitudeNorth),]

count_mat = count(df = as.data.frame(all.preds1),vars = "SiteID")

add1 = matrix(0,nrow(all.preds1),1)

for(i in 1:nrow(all.preds1)){
  add1[i,] = count_mat[count_mat[,1]==all.preds1$SiteID[i],2]
}

rbPal <- colorRampPalette(c("yellow",'red','blue',"green","black"))

#This adds a column of color values
# based on the y values
color_dat <- rbPal(5)[as.numeric(cut(add1,breaks = 5))]

breaks <-  c(0,5, 10,15, 20,25,30,35, 40,45, 50)
colors <- rbPal(length(breaks)-1)
data_binned <- as.character(cut(add1, breaks,
                                include.lowest = TRUE,
                                labels=colors))

map('state', xlim=range(all.preds1[,4])+c(-2, 2), ylim=range(all.preds1[,3])+c(-1, 1))
points(all.preds1[,4],all.preds1[,3], pch=19, cex=1,col=data_binned)
title("Colored by Number of Samples")

hist(add1,col=colors,breaks=10,xlab="Number of samples between Ages 200 - 1000")

breaks <-  seq(200,1000,100)
colors <- rbPal(length(breaks)-1)
data_binned1 <- as.character(cut(all.preds1$Age, breaks,
                                include.lowest = TRUE,labels=colors))

map('state', xlim=range(all.preds1[,4])+c(-2, 2), ylim=range(all.preds1[,3])+c(-1, 1))
points(all.preds1[,4],all.preds1[,3], pch=19, cex=1,col=data_binned1)
title("Colored by Max Age")

hist(all.preds1$Age,col=colors,breaks=8,xlab="Sample Ages")


#pdf("time.series1_with_pie.pdf")
quartz()
par(mfrow=c(2,2))
for(i in 1:142){
  if(length(all.preds1[all.preds1[,1]==i,7])>5){
    plot(all.preds1[all.preds1[,1]==i,7],all.preds1[all.preds1[,1]==i,9],xlab="Age",ylab="Biomass",
         main = c("site",as.character(unique(all.preds1[all.preds1[,1]==i,2])[1])),
         ylim=c(0,400),xlim=c(150,1000),pch=19,cex=1)
    map('state', xlim=range(all.preds1[,4])+c(-2, 2), ylim=range(all.preds1[,3])+c(-1, 1))
    points(unique(all.preds1[all.preds1[,1]==i,4]),unique(all.preds1[all.preds1[,1]==i,3]), pch=19, cex=1)
   
    min.calc = min(all.preds1[all.preds1[,1]==i,7])
    min.plot = as.numeric(all.preds1[which(all.preds1[,1]==i&all.preds1[,7]==min.calc),
                                     11:ncol(all.preds1)])
    pie(min.plot,col=rainbow(ncol(ten.count)),main=c("Age BP",min.calc),
        labels = colnames(all.preds1[,11:ncol(all.preds1)]))
    max.calc = max(all.preds1[all.preds1[,1]==i,7])
    max.plot = as.numeric(all.preds1[which(all.preds1[,1]==i&all.preds1[,7]==max.calc),
                                     11:ncol(all.preds1)])
    pie(max.plot,col=rainbow(ncol(ten.count)), main=c("Age BP",max.calc),
        labels = colnames(all.preds1[,11:ncol(all.preds1)]))
  }
}
#dev.off()

plot.pie.time.series = function(site){
  plot(all.preds1[all.preds1[,2]==site,7],all.preds1[all.preds1[,2]==site,9],xlab="Age",ylab="Biomass",
       main = c("site",as.character(unique(all.preds1[all.preds1[,2]==site,2])[1])),
       ylim=c(0,400),pch=19,cex=1)
  text(all.preds1[all.preds1[,2]==site,7],all.preds1[all.preds1[,2]==site,9]+40,labels = all.preds1[all.preds1[,2]==site,7],cex=.6)
  map('state', xlim=range(all.preds1[,4])+c(-2, 2), ylim=range(all.preds1[,3])+c(-1, 1))
  points(unique(all.preds1[all.preds1[,2]==site,4]),unique(all.preds1[all.preds1[,2]==site,3]), pch=19, cex=1)
  
  which.pie = all.preds1[all.preds1[,2]==site,]
  #min.plot = as.numeric(all.preds1[which(all.preds1[,1]==i&all.preds1[,7]==min.calc),
  #                                 11:ncol(all.preds1)])
  for(i in 1:nrow(which.pie)){
    pie(as.numeric(which.pie[i,11:ncol(which.pie)]),col=rainbow(ncol(ten.count)),main=c("Age BP",which.pie$Age[i]),
        labels = colnames(all.preds1[,11:ncol(all.preds1)])) 
  }
}

quartz()
#pdf("site_1988.pdf")
par(mfrow=c(3,3))
plot.pie.time.series(site = 1988)
#dev.off()

new.site.locs <- cbind(all.preds1$LongitudeWest,all.preds1$LatitudeNorth)
centers_pol = data.frame(new.site.locs)
colnames(centers_pol) = c('x', 'y')

coordinates(centers_pol) <- ~ x + y
proj4string(centers_pol) <- CRS('+proj=longlat +ellps=WGS84')

centers_polA <- spTransform(centers_pol, CRS('+init=epsg:3175'))
centers_polA <- as.matrix(data.frame(centers_polA))

all.preds1 = cbind(all.preds1[,1:10],centers_polA)

#colnames(all.preds1[,8:10]) <- c("2.5%","50%","97.5")
head(all.preds1)
colnames(all.preds1)<-c("site_factor","SiteID","Lat","Long","dataset.id",
                        "ContactName","Age","low_bound","Median","high_bound","x","y")

b <- gam(log(Median) ~ te(x,y,Age, d = c(2,1),bs = c("tp","cr"),k=30), data = as.data.frame(all.preds1))
summary(b)
vis.gam(b)  

load("/Users/paleolab/Documents/babySTEPPS/biomass_dat5.Rdata")
age_slice = 300

diff_map_plor <- function(age_slice){

  pred_data = cbind(biomass_dat5[,1:2],rep(age_slice,nrow(biomass_dat5)))
  colnames(pred_data)<- c("x","y","Age")
  pred_biomass_gam = exp(predict(b,newdata = as.data.frame(pred_data)))
  #hist(pred_biomass_gam)

for(i in 1:length(pred_biomass_gam)){
  if(pred_biomass_gam[i]>400) pred_biomass_gam[i] = 400
}

full.mat <- cbind(biomass_dat5[,1:2],as.vector(pred_biomass_gam))
colnames(full.mat) <- c("x","y","pred biomass")
y = as.data.frame(full.mat) #rowSums(biomass_dat_est) to make xiaopings

breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)-1))
legendName <- c("Biomass (Mg/ha)")#paste0("Biomass at Age = ",age_slice, " BP")

biomass_dat_est <- read.csv(paste0(data.dir,"biomass_prediction_v0.2.csv"))
xiao_ests <- rowSums(biomass_dat_est[,4:23])
breaks <-  c(-300, -200, -100, 100, 200, 300, 400, 500, 600)
colors <- c("blue", "light blue", "white", "pink","darkred", "red","orange","yellow")
legendName <- paste0("Xiao - Pred at ", age_slice)

data_binned <-  cut(xiao_ests, breaks, include.lowest = TRUE, labels = FALSE)

breaklabels <- apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1,  function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })

inputData <- data.frame(X = y[,1], Y = y[,2], Preds = cbind(data_binned,data_binned))
inputData_long <- melt(inputData, c('X', 'Y'))

colnames(centers_polA) <- c('lat','lon')
input_points <- data.frame(centers_polA[all.preds1$Age >(age_slice-50)&all.preds1$Age <(age_slice+50),])

d <- ggplot() + geom_raster(data = inputData_long, aes(x = X, y = Y, fill = factor(value))) +
  scale_fill_manual(labels = breaklabels, name = legendName, drop = FALSE, values = colors, guide = "legend") + 
  theme(strip.text.x = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16)) + #legend.position="none" removes legend
  geom_point(data = input_points, aes(x=lat,y=lon), pch=16, size=2,colour="black") +
  ggtitle("Xiaoping's Biomass Estimates")#paste0("Pred Map ","at Age = ",age_slice, " BP")

add_map_albers <- function(plot_obj, map_data = usFortified, dat){
  p <- plot_obj + geom_path(data = map_data, aes(x = long, y = lat, group = group), size = 0.1) +
    scale_x_continuous(limits = c(min(inputData$X, na.rm = TRUE), max(inputData$X, na.rm = TRUE))) +
    scale_y_continuous(limits = c(min(inputData$Y, na.rm = TRUE), max(inputData$Y, na.rm = TRUE)))
  return(p)
}

d <- add_map_albers(plot_obj = d, map_data = usFortified, dat = inputData_long)
return(d)
}
quartz()
print(d)


d_200 = diff_map_plor(200)
d_300 = diff_map_plor(300)
d_400 = diff_map_plor(400)
d_500 = diff_map_plor(500)

d_600 = diff_map_plor(600)
d_700 = diff_map_plor(700)
d_800 = diff_map_plor(800)
d_900 = diff_map_plor(900)

quartz()
pdf("initial.pred.diff.maps.pdf")
grid.arrange(d_200,d_300)
grid.arrange(d_400,d_500)
grid.arrange(d_600,d_700)
grid.arrange(d_800,d_900)
dev.off()

## 4/4/2015 TO DO
#figure out how to add covariates
#run model for longer time series