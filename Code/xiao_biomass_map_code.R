library(reshape)
library(ggplot2)
library(sp)
library(rgdal)
library(fields)
library(maptools)
library(neotoma)
require(grid)
require(plyr)
require(maps)
library(gridExtra)
gpclibPermit()
library(mgcv)
library(splines)
library(boot)
library(gtools)
library(rjags)

biomass_dat_est <- read.csv(paste0("~/Downloads/","biomass_prediction_v0.9-10_bam.csv"))
xiao_ests <- biomass_dat_est$Total#rowSums(biomass_dat_est[,4:23])

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

full.mat <- cbind(biomass_dat_est[,c('x','y')],xiao_ests)
colnames(full.mat) <- c("x","y","Xiao Total Biomass")
y = as.data.frame(full.mat)

breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)-1))

legendName <- "Biomass (Mg/ha)"

data_binned <-  cut(y[,3], breaks, include.lowest = TRUE, labels = FALSE)

breaklabels <- apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1,  function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })

inputData <- data.frame(X = y[,1], Y = y[,2], Preds = cbind(data_binned,data_binned))
inputData_long <- melt(inputData, c('X', 'Y'))

input_points <- data.frame(albers) # how to add points part A
colnames(input_points) <- c('lat','lon')

d <- ggplot() + geom_raster(data = inputData_long, aes(x = X, y = Y, fill = factor(value))) + scale_fill_manual(labels = breaklabels, name = legendName, drop = FALSE, values = colors, guide = "legend") + 
  theme(strip.text.x = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16)) +
  geom_point(data = input_points, aes(x=lat,y=lon), pch=16, size=2,colour="black") + # how to add points part B
  ggtitle("Xiaoping Estimates")

add_map_albers <- function(plot_obj, map_data = usFortified, dat){
  p <- plot_obj + geom_path(data = map_data, aes(x = long, y = lat, group = group), size = 0.1) +
    scale_x_continuous(limits = c(min(inputData$X, na.rm = TRUE), max(inputData$X, na.rm = TRUE))) +
    scale_y_continuous(limits = c(min(inputData$Y, na.rm = TRUE), max(inputData$Y, na.rm = TRUE)))
  return(p)
}

d <- add_map_albers(plot_obj = d, map_data = usFortified, dat = inputData_long)

quartz()
pdf('biomass_est_1.pdf')
print(d)
dev.off()
