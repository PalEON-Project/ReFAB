load("all.preds.min.list.RData")
#follows from Run_Model.R ##all.preds = summary(csamp.real.pred)

site.factor = factor(x[,1],labels = seq(1,142,1))
all.preds1 = cbind(site.factor,x[,1:6],all.preds$quantiles[,1],all.preds$quantiles[,3],all.preds$quantiles[,5])
all.preds1 = all.preds1[order(all.preds1$LatitudeNorth),]
all.preds1[,1] = factor(all.preds1[,2],labels = seq(1,142,1),
                        levels=unique(all.preds1[,2]),ordered=FALSE)

pdf("time.series1.pdf")
par(mfrow=c(4,4))
for(i in 1:142){
  if(length(all.preds1[all.preds1[,1]==i,7])>5){
    plot(all.preds1[all.preds1[,1]==i,7],all.preds1[all.preds1[,1]==i,9],xlab="Age",ylab="Biomass",
         main = c("site",as.character(unique(all.preds1[all.preds1[,1]==i,2])[1])),
         ylim=c(0,400),xlim=c(150,1000),pch=19,cex=1)
    map('state', xlim=range(all.preds1[,4])+c(-2, 2), ylim=range(all.preds1[,3])+c(-1, 1))
    points(unique(all.preds1[all.preds1[,1]==i,4]),unique(all.preds1[all.preds1[,1]==i,3]), pch=19, cex=1)
  }
}
dev.off()

new.site.locs <- cbind(all.preds1$LongitudeWest,all.preds1$LatitudeNorth)
centers_pol = data.frame(new.site.locs)
colnames(centers_pol) = c('x', 'y')

coordinates(centers_pol) <- ~ x + y
proj4string(centers_pol) <- CRS('+proj=longlat +ellps=WGS84')

centers_polA <- spTransform(centers_pol, CRS('+init=epsg:3175'))
centers_polA <- as.matrix(data.frame(centers_polA))

all.preds1 = cbind(all.preds1,centers_polA)

#colnames(all.preds1[,8:10]) <- c("2.5%","50%","97.5")
head(all.preds1)
colnames(all.preds1)<-c("site_factor","SiteID","Lat","Long","dataset.id",
                        "ContactName","Age","low_bound","Median","high_bound","x","y")

b <- gam(log(Median) ~ te(x,y,Age), data = as.data.frame(all.preds1), d = c(2,1),bs = c("tp","cr"))
summary(b)
vis.gam(b)  

load("/Users/paleolab/Documents/babySTEPPS/biomass_dat5.Rdata")
age_slice = 300
pred_data = cbind(biomass_dat5[,1:2],rep(age_slice,nrow(biomass_dat5)))
colnames(pred_data)<- c("x","y","Age")
pred_biomass_gam = exp(predict(b,newdata = as.data.frame(pred_data)))

full.mat <- cbind(biomass_dat5[,1:2],as.vector(pred_biomass_gam))
colnames(full.mat) <- c("x","y","pred biomass")
y = as.data.frame(full.mat) #rowSums(biomass_dat_est) to make xiaopings

breaks <-  c(0, 25, 50, 75, 100, 125, 150, 200, 250, 300, 400)
colors <- rev(terrain.colors(length(breaks)-1))

legendName <- "Biomass at Age = 300"

data_binned <-  cut(y[,3], breaks, include.lowest = TRUE, labels = FALSE)

breaklabels <- apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1,  function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })

inputData <- data.frame(X = y[,1], Y = y[,2], Preds = cbind(data_binned,data_binned))
inputData_long <- melt(inputData, c('X', 'Y'))

colnames(centers_polA) <- c('lat','lon')
input_points <- data.frame(centers_polA[all.preds1$Age == age_slice,])

d <- ggplot() + geom_raster(data = inputData_long, aes(x = X, y = Y, fill = factor(value))) + scale_fill_manual(labels = breaklabels, name = legendName, drop = FALSE, values = colors, guide = "legend") + 
  theme(strip.text.x = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16)) +
  geom_point(data = input_points, aes(x=lat,y=lon), pch=16, size=2,colour="black") +
  ggtitle("Prediction Map")

add_map_albers <- function(plot_obj, map_data = usFortified, dat){
  p <- plot_obj + geom_path(data = map_data, aes(x = long, y = lat, group = group), size = 0.1) +
    scale_x_continuous(limits = c(min(inputData$X, na.rm = TRUE), max(inputData$X, na.rm = TRUE))) +
    scale_y_continuous(limits = c(min(inputData$Y, na.rm = TRUE), max(inputData$Y, na.rm = TRUE)))
  return(p)
}

d <- add_map_albers(plot_obj = d, map_data = usFortified, dat = inputData_long)

quartz()
print(d)












