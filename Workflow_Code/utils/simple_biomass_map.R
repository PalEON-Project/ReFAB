
library(ncdf4)
library(raster)

nc <- nc_open(file.path('~','Downloads','PLS_biomass_agb_western_point_v1.0rc2.nc'))
#nc_fia <- nc_open('~/Downloads/FIA_biomass_point_v0.999.nc')

x <- nc$dim$x$vals
y <- nc$dim$y$vals
data <- ncvar_get(nc, varid = c('Total'))

rownames(data) <- x
colnames(data) <- y

r1 <- raster(list(x = x, y = y, z = data))

#can do
#plot(r1)
breaks <-  c(0,25, 50,75,100,125,150,175,200,225,250,300,800)#c(0, seq(25, 100, 25), seq(150, 300, 50), max(y$Tot_Biomass))

colors <- rev(terrain.colors(length(breaks))[-4])

biomass_dat_est <- as.data.frame(rasterToPoints(r1))
colnames(biomass_dat_est) <- c('x', 'y', 'Total')
xiao_ests <- biomass_dat_est$Total


### Biomass
full.mat <- cbind(biomass_dat_est[, c('x', 'y')], xiao_ests)
colnames(full.mat) <- c("x", "y", "Tot_Biomass")
tot.biom.df <- y <- as.data.frame(full.mat)

data_binned <-
  cut(full.mat$Tot_Biomass,
      breaks,
      include.lowest = TRUE,
      labels = FALSE)

inputData <- data.frame(
  X = full.mat$x,
  Y = full.mat$y,
  Preds = data_binned,
  tot_biomass = full.mat$Tot_Biomass
)

inputData<-inputData[inputData$Y>500000,]
inputData<-inputData[inputData$X<1200000,]

lat.long.reg.df = data.frame(inputData[,c('X','Y')])
colnames(lat.long.reg.df) = c('x', 'y')

coordinates(lat.long.reg.df) <- ~ x + y
proj4string(lat.long.reg.df) <-  CRS('+init=epsg:3175')
inputData_latlon <- spTransform(lat.long.reg.df, CRS('+proj=longlat +ellps=WGS84'))
inputData_latlon <- as.matrix(data.frame(inputData_latlon))
inputData_latlon <- cbind(inputData_latlon,inputData[,c('Preds','tot_biomass')])

colnames(inputData_latlon) <- c('X','Y','Preds','tot_biomass')

inputData = inputData_latlon


#' rev(
#'   c(
#'     '#d73027',
#'     '#f46d43',
#'     '#fdae61',
#'     '#fee090',
#'     '#ffffbf',
#'     '#e0f3f8',
#'     '#abd9e9',
#'     '#74add1',
#'     #'#4575b4',
#'     cols[3]
#'   )
#' )#)) #

legendName <- "Biomass (Mg/ha)"

data_binned <-
  cut(y$Tot_Biomass,
      breaks,
      include.lowest = TRUE,
      labels = FALSE)

breaklabels <- apply(cbind(breaks[1:(length(breaks) - 1)],
                           breaks[2:length(breaks)]), 1,
                     function(r) {
                       sprintf("%0.2f - %0.2f",
                               r[1], r[2])
                     })

breaklabels <- signif(breaks[2:length(breaks)],digits = 2)

YBP100 <- read.csv('refab_final_datasets/data_products/ReFAB_spatial_reconstruction_means_v1.0.csv')
data_binned <-
  cut(YBP100$AgeYBP1000,
      breaks,
      include.lowest = TRUE,
      labels = FALSE)

inputData <- data.frame(
  X = YBP100$x,
  Y = YBP100$y,
  Preds = data_binned,
  tot_biomass = YBP100$AgeYBP1000
)

data_binned <-
  cut(as.numeric(inputData$tot_biomass),
      breaks,
      include.lowest = TRUE,
      labels = FALSE)

calib_meta <- read.csv('refab_final_datasets/data_products/ReFAB_calibration_data_v1.2.csv') 

input_points <- data.frame(calib_meta[,c('lat','long','in.calib.sample','calib_pred')]) # how to add points part A
colnames(input_points) <- c('lat','lon')

pdf(paste0(Sys.Date(), 'new_simple_biomass_map.pdf'))
plot(
  inputData$X,
  inputData$Y,
  col = colors[data_binned],
  pch = 19,
  xlab = NA,
  ylab = NA,ylim=c(41.5,49.5)
)
maps::map('state', add = T)
points(
  input_points$lon,
  input_points$lat,
  col = c('magenta3', 'darkblue')[input_points[, 3] + 1],
  bg = colors[cut(input_points[,4],breaks,
                  include.lowest = TRUE,
                  labels = FALSE)],
  #bg = c('magenta3', 'darkblue')[input_points[, 3] + 1],
  pch = c(25,21)[input_points[, 3] + 1],
  cex = 1
)
dev.off()

plot.new()
legend('center',breaklabels,col=colors,pch=19)

# 
# 
# 
#  d <-  ggplot() +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.background = element_blank(),
#     legend.text=element_text(size=18),
#     legend.title=element_text(size=20),
#     title = element_text(size=0),
#     aspect.ratio = 1,
#    # axis.ticks = element_blank(),
#     #axis.text = element_blank()
#   ) +
#   geom_point(data = inputData_latlon, aes(
#     x = X,
#      y = Y,
#      color = factor(data_binned)
#   )) +
#   scale_color_manual(
#     labels = breaklabels,
#     name = legendName,
#     drop = FALSE,
#     values = colors,
#     guide = "legend"
#   ) +
#   geom_point(
#     data = input_points,
#     aes(x = lat, y = lon),
#     pch = 16,
#     size = 1.25,
#     alpha = .9 ,
#     colour = "black"
#   )
# 
# add_map_albers <- function(plot_obj, map_data = usFortified, dat) {
#   p <-
#     plot_obj + geom_path(data = map_data,
#                          aes(x = long, y = lat, group = group),
#                          size = .5) +
#     scale_x_continuous(limits = c(min(inputData$X, na.rm = TRUE), max(inputData$X, na.rm = TRUE))) +
#     scale_y_continuous(limits = c(min(inputData$Y, na.rm = TRUE), max(inputData$Y, na.rm = TRUE)))
#   return(p)
# }
# 
# d <-
#   add_map_albers(plot_obj = d,
#                  map_data = usFortified,
#                  dat = inputData_long)
# #quartz()
# print(d)
# 
# ggsave(d, filename = paste0(Sys.Date(), 'simple_biomass_map.pdf'))
usShp <-
  rgdal::readOGR(
    file.path("~/Downloads/", 'us_alb.shp')#,
    # proj4string = CRS('+init=epsg:3175')
  )
usShp@data$id <- rownames(usShp@data)
usFortified <- fortify(usShp, region = 'id')

##### add state lines function
add_map_albers <- function(plot_obj, map_data = usShp, dat) {
  p <-
    plot_obj + geom_path(data = map_data,
                         aes(x = long, y = lat, group = group),
                         size = 1) +
    scale_x_continuous(limits = c(min(dat$x, na.rm = TRUE), max(dat$x, na.rm = TRUE))) +
    scale_y_continuous(limits = c(min(dat$y, na.rm = TRUE), max(dat$y, na.rm = TRUE)))
  return(p)
}

##### make a heat plot of biomass
theme_clean <- function(plot_obj) {
  plot_obj <- plot_obj + theme(
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
  
  return(plot_obj)
}

