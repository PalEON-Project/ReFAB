
###copied from gam_spatial.R

diff_do <-
  coors_dat_esa[which(!is.na(ove)), 3] - coors_dat_pls[which(!is.na(ove_pls)), 3]
diff_breaks <- seq(-400, 400, 10)
diff_cuts <-
  cut(
    diff_do,
    breaks = diff_breaks,
    include.lowest = T,
    labels = F
  )
diff_colors <-
  colorRampPalette(c('purple', 'blue', 'white', 'red', 'maroon'))(length(diff_breaks) -
                                                                    1)

plot(
  coors_dat_esa[which(!is.na(ove)), 1],
  coors_dat_esa[which(!is.na(ove)), 2],
  pch = 15,
  col = diff_colors[diff_cuts],
  main = 'ESA minus PLS',
  ylab = NA,
  xlab = NA,
  cex=.75
)
maps::map('state',add=T)
legend('topright',
       as.character(pretty(diff_breaks, n = 8)),
       col = diff_colors[c(which(diff_breaks %in% pretty(diff_breaks, n = 8))[-9], length(diff_colors))],
       pch = 15)

hist(diff_do, main = 'Histogram of ESA minius PLS', xlab = 'Biomass Difference',breaks=1000)

require(sf)
lu <- raster('~/Downloads/NLCD_2011_Land_Cover_L48_20190424/NLCD_2011_Land_Cover_L48_20190424.img')
### transform pts instead

library(raster)
# create spatial points data frame
spg <- data.frame(coors_dat_esa[which(!is.na(ove)),1:2])
colnames(spg) <- c('x','y')
coordinates(spg) <- ~ x + y
crs(spg) <- "+proj=longlat +ellps=WGS84 +no_defs"

spg_aea <- spTransform(spg,CRSobj = CRS('+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,-0,-0,-0,0 +units=m +no_defs'))

ncld_vals <- raster::extract(lu,spg)



diffdat <- data.frame(diff_do,ncld_vals)

#unique vals legend names
class_nums <- unique(diffdat$ncld_vals)

ncld_class <- c('open water','woody wetlands','cultivated crops','decidous forest','pasture/hay',
                'emergent herbaceous wetlands','evergreen forest','developed, open space',
                'developed, medium intensity','grassland/herbaceous','developed,low intensity',
                'shrub/scrub','developed high intensity','mixed forest','barren land')

pdf('diff_hists_by_nlcd.pdf',height=12,width=12)
par(mfrow=c(4,4),oma=c(1,7,1,1))
barplot(table(factor(x = diffdat$ncld_vals ,levels = class_nums,labels = ncld_class)),las=2,horiz=T)
plot.new()
for(ii in 1:length(ncld_class)){
  val_do <- unique(diffdat$ncld_vals)[ii]
  ncld_cols <- rep('gray',length(ncld_vals))
  ncld_cols[ncld_vals==val_do] <- 'darkgreen'
  plot(spg,col=ncld_cols,main=ncld_class[ii])
  
  hist(diffdat[which(diffdat$ncld_vals==val_do),'diff_do'],xlim=c(-400,400),breaks=diff_breaks,col=diff_colors,main=ncld_class[ii],border='gray',xlab = 'PLS to ESA difference')
}
dev.off()


diffacdat <- data.frame(fact = factor(x = diffdat$ncld_vals ,levels = class_nums,labels = ncld_class),diff_do = diffdat$diff_do)


sumstatsdiff <- diffacdat %>% group_by(fact) %>% summarise(sum = sum(diff_do),mean = mean(diff_do),median = median(diff_do),sd=sd(diff_do),pct_inc = length(which(diff_do>0))/length(diff_do),pct_dec = length(which(diff_do<0))/length(diff_do))

sumstatsdiff[order(sumstatsdiff$sum),]

# coerce to SpatialPixelsDataFrame
#gridded(spg) <- TRUE
# coerce to raster
rasterDF <- raster(spg,vals=NA)
cells <- cellFromXY(rasterDF, coors_dat_esa[,1:2])

rasterDF <- writeRaster(rasterDF, filename='test.grd',overwrite=T)

rasterDF <- update(rasterDF, cell=cells, v=as.numeric(coors_dat_esa[, 3] - coors_dat_pls[, 3]))

crs(rasterDF) <- CRS("+proj=longlat +ellps=WGS84 +no_defs")#('+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,-0,-0,-0,0 +units=m +no_defs')
r.new <- projectRaster(rasterDF,crs='+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,-0,-0,-0,0 +units=m +no_defs')
lu.new <- crop(lu,r.new)
x <- resample(lu.new,r.new,method='bilinear')

