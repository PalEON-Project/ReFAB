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

library(ncdf4)
library(raster)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 3
cols = gg_color_hue(n)



#biomass_dat_est <- read.csv(paste0("~/Downloads/","biomass_prediction_v0.9-10_bam.csv"))

nc <-
  nc_open(file.path('Data', 'PLS_biomass_western_point_v0.999.nc'))
nc_fia <- nc_open('~/Downloads/FIA_biomass_point_v0.999.nc')


x <- nc$dim$x$vals
y <- nc$dim$y$vals
data <- ncvar_get(nc, varid = c('Total'))

biomass_hemlock <- ncvar_get(nc, varid = c('Hemlock'))
biomass_beech <- ncvar_get(nc, varid = c('Beech'))

rownames(data) <- x
colnames(data) <- y

r1 <- raster(list(x = x, y = y, z = data))

r_hem <- raster(list(x = x, y = y, z = biomass_hemlock))
r_beech <- raster(list(x = x, y = y, z = biomass_beech))

#can do
#plot(r1)

biomass_dat_est <- as.data.frame(rasterToPoints(r1))
colnames(biomass_dat_est) <- c('x', 'y', 'Total')
xiao_ests <- biomass_dat_est$Total

#### FIA
x_fia <- nc_fia$dim$x$vals
y_fia <- nc_fia$dim$y$vals
data_fia <- ncvar_get(nc_fia, varid = c('Total'))

rownames(data_fia) <- x_fia
colnames(data_fia) <- y_fia

r1_fia <- raster(list(x = x_fia, y = y_fia, z = data_fia))

biomass_dat_est_fia <- as.data.frame(rasterToPoints(r1_fia))
colnames(biomass_dat_est_fia) <- c('x', 'y', 'FIA')


all <- merge(biomass_dat_est,biomass_dat_est_fia, by = c('x','y'))

plot(all$Total,all$FIA)
all<-all[all$y>500000,]
all<-all[all$x<1200000,]
hist((all$Total-all$FIA))

hist((all$FIA)/(all$Total),xlim=c(0,1))

plot(all$x[all$Total>75],all$y[all$Total>75])

1-sum(all$FIA[all$Total>75])/sum(all$Total[all$Total>75])

hem_beech_mat <-
  cbind(rasterToPoints(r_hem), rasterToPoints(r_beech)) #
hem.beech <- hem_beech_mat[,3]#rowSums(hem_beech_mat[, c(3, 6)])

usShp <-
  readShapeLines(
    file.path("/Users/paleolab/Documents/babySTEPPS/", 'us_alb.shp'),
    proj4string = CRS('+init=epsg:3175')
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

### Biomass
full.mat <- cbind(biomass_dat_est[, c('x', 'y')], xiao_ests, rowSums(hem_beech_mat[, c(3, 6)]))
colnames(full.mat) <- c("x", "y", "Tot_Biomass", 'LS')
tot.biom.df = as.data.frame(full.mat)

### Density

nc_d <- nc_open('~/Downloads/PLS_density_western_point_v0.999.nc')
nc_d <- nc_open('~/Downloads/FIA_density_point_v0.999.nc')


x <- nc_d$dim$x$vals
y <- nc_d$dim$y$vals
data <- ncvar_get(nc_d, varid = c('Total'))
r2 <- raster(list(x = x, y = y, z = data))

density.df = as.data.frame(rasterToPoints(r2))
colnames(density.df) <- c('x', 'y', 'density')

y <- merge(tot.biom.df, density.df, by = c('x', 'y'))

breaks <-  c(0, seq(25, 100, 25), seq(150, 300, 50), max(y$Tot_Biomass))
colors <- rev(
  c(
    '#d73027',
    '#f46d43',
    '#fdae61',
    '#fee090',
    '#ffffbf',
    '#e0f3f8',
    '#abd9e9',
    '#74add1',
    #'#4575b4',
    cols[3]
  )
)#)) #rev(topo.colors(length(breaks)-1))

density.breaks <- c(0, 1, 48, max(y$density))

legendName <- "Biomass (Mg/ha)"

data_binned <-
  cut(y$Tot_Biomass,
      breaks,
      include.lowest = TRUE,
      labels = FALSE)
data_binned_density <-
  cut(y$density,
      density.breaks,
      include.lowest = TRUE,
      labels = FALSE)

breaklabels <- apply(cbind(breaks[1:(length(breaks) - 1)],
                           breaks[2:length(breaks)]), 1,
                     function(r) {
                       sprintf("%0.2f - %0.2f",
                               r[1], r[2])
                     })

inputData <- data.frame(
  X = y$x,
  Y = y$y,
  Preds = data_binned,
  tot_biomass = y$Tot_Biomass,
  density = y$density,
  LS = y$LS
  #beech = y$beech
)

inputData<-inputData[inputData$Y>500000,]
inputData<-inputData[inputData$X<1200000,]

data_binned <-
  cut(as.numeric(inputData$tot_biomass),
      breaks,
      include.lowest = TRUE,
      labels = FALSE)
data_binned_density <-
  cut(inputData$density,
      density.breaks,
      include.lowest = TRUE,
      labels = FALSE)

inputData_long <- melt(inputData, c('X', 'Y'))

load('2018-11-12calibration.albers.Rdata')
input_points <- data.frame(albers) # how to add points part A
colnames(input_points) <- c('lat', 'lon')

inputdata1 <-
  inputData[inputData$density > 47, ] #ask kelly about this number
colnames(inputdata1)[1:2] <- c('a', 'b')
col.plot <- rep('Forest', nrow(inputdata1))
data_binned1 <-
  cut(inputdata1$tot_biomass,
      breaks,
      include.lowest = TRUE,
      labels = FALSE)

inputdata2 <- inputData[inputData$LS > 1, ] #what should the cutoff be?
colnames(inputdata2)[1:2] <- c('c', 'd')
col.plot2 <- rep('LS', nrow(inputdata1))
data_binned2 <-
  cut(inputdata2$tot_biomass,
      breaks,
      include.lowest = TRUE,
      labels = FALSE)

#inputdata4 <- inputData[inputData$beech > 1, ] #what should the cutoff be?
#colnames(inputdata4)[1:2] <- c('g', 'h')
#col.plot4 <- rep('Beech', nrow(inputdata1))
#data_binned4 <-
#  cut(inputdata4$tot_biomass,
 #     breaks,
  #    include.lowest = TRUE,
   #   labels = FALSE)


inputdata3 <-
  inputData[inputData$X > 430000 &
              inputData$X < 557000 & inputData$Y > 678000 & inputData$Y < 1200000, ]
colnames(inputdata3)[1:2] <- c('e', 'f')

d <-   ggplot() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.text=element_text(size=18),
    legend.title=element_text(size=20),
    title = element_text(size=0),
    aspect.ratio = 1,
    axis.text = element_text(size=15)
    
    
  ) +
  geom_tile(data = inputData, aes(
    x = X,
    y = Y,
    fill = factor(data_binned)
  )) +
  
  scale_fill_manual(
    labels = breaklabels,
    name = legendName,
    drop = FALSE,
    values = colors,
    guide = "legend"
  ) +
  geom_tile(
    data = inputdata1,
    aes(x = a, y = b, colour = 'Forest'),
    color = cols[1],
    alpha = 0,
    size = 4
  ) +
  geom_tile(
    data = inputdata1,
    aes(x = a, y = b, colour = 'Forest'),
    color = 'white',
    alpha = 0,
    size = .5
  ) +
  
  #scale_color_manual(values='black', name=element_blank())
  geom_tile(data = inputdata1, aes(
    x = a,
    y = b,
    fill = factor(data_binned1)
  )) + #puts biomass colors on
  
  
 #puts biomass colors on
  geom_tile(
    data = inputdata2,
    aes(x = c, y = d, colour = col.plot2),
    color = cols[2],
    alpha = 0,
    size = 4
  ) +
  geom_tile(
    data = inputdata2,
    aes(x = c, y = d, colour = col.plot2),
    color = 'gray',
    alpha = 0,
    size = .5
  ) +
  geom_tile(data = inputdata2, aes(
    x = c,
    y = d,
    fill = factor(data_binned2)
  )) +

  geom_tile(
    data = inputdata3,
    aes(x = e, y = f, color = 'Transect'),
    fill = 'gray',
    color = NA,
    size = .5,
    alpha = .5
  ) +
  geom_point(
    data = input_points,
    aes(x = lat, y = lon),
    pch = 16,
    size = 1.25,
    alpha = .9 ,
    colour = "black"
  )


#+ # how to add points part B
# geom_text(data = input_points,aes(x=lat,y=lon,label=unlist(name.keep))) +
# ggtitle("Xiaoping Estimates")

add_map_albers <- function(plot_obj, map_data = usFortified, dat) {
  p <-
    plot_obj + geom_path(data = map_data,
                         aes(x = long, y = lat, group = group),
                         size = 0.1) +
    scale_x_continuous(limits = c(min(inputData$X, na.rm = TRUE), max(inputData$X, na.rm = TRUE))) +
    scale_y_continuous(limits = c(min(inputData$Y, na.rm = TRUE), max(inputData$Y, na.rm = TRUE)))
  return(p)
}

d <-
  add_map_albers(plot_obj = d,
                 map_data = usFortified,
                 dat = inputData_long)
#quartz()
print(d)

ggsave(d, filename = paste0(Sys.Date(), 'fia_biomass_map.pdf'))

temp <-
  read.csv('~/Downloads/reclimatedata/tmean_yr_Prism_1900_full.csv')
precip <-
  read.csv('~/Downloads/reclimatedata/pr_monthly_Prism_1900_full.csv')

temp.precip <- merge(x = temp, y = precip, by = c('x', 'y'))
colnames(inputdata3)[1:2] <- c('x', 'y')
temp.precip.bio <- merge(temp.precip, inputdata3, by = c('x', 'y'))

#temp.precip.bio$Mean*temp.precip.bio$total #smooth metric

LS_binned <- rep('Prairie/Savanna', length(temp.precip.bio$LS))
LS_binned[temp.precip.bio$density > 47] <- 'Forested (No Hemlock)'
LS_binned[temp.precip.bio$LS > 1] <- 'Hemlock'

colors_use <- c(cols[1:2], colors[1])

A <-
  ggplot() + theme_bw() + theme(legend.text=element_text(size=15),
                                legend.title=element_text(size=17),
                                title = element_text(size=18),
                                axis.text = element_text(size = 12)) +
  labs(x = '', y = "Biomass (Mg/ha)", color = "Ecosystem Type") +
  geom_jitter(data = temp.precip.bio,
              aes(x = y, y = tot_biomass, color = LS_binned),
              size = 2)

B <-
  ggplot() + theme_bw() + theme(legend.text=element_text(size=15),
                              legend.title=element_text(size=17),
                              title = element_text(size=18),
                              axis.text = element_text(size = 12)) +
  labs(x = '', y = "Mean Temperature (ºC)", color = "Ecosystem Type") +
  geom_jitter(data = temp.precip.bio,
              aes(x = y, y = Mean, color = LS_binned),
              size = 2)#+ stat_smooth(data = temp.precip.bio, aes(x = y, y = Mean), method = "lm",
# formula = y ~ poly(x, 10), se = FALSE) #mean temp
C <-
  ggplot() + theme_bw() + theme(legend.text=element_text(size=15),
                                legend.title=element_text(size=17),
                                title = element_text(size=18),
                                axis.text = element_text(size = 12)) +
  labs(x = "Latitude (Albers Projection)", y = "Total Precip (cm)", color = "Ecosystem Type") +
  geom_jitter(data = temp.precip.bio,
              aes(x = y, y = total, color = LS_binned),
              size = 2)#+ stat_smooth(data = temp.precip.bio, aes(x = y, y = total), method = "lm",
# formula = y ~ poly(x, 10), se = FALSE) #total precip

library(egg)
F <- ggarrange(A, B, C, ncol = 1)
#print(F)
ggsave(F, file = paste0(Sys.Date(), 'transect_output_scatters.pdf'))

#### Mapping to potential vegetation

rast <- raster::raster('~/Downloads/potentialveg/potveg/w001001.adf')

if(FALSE){
pv <- read.asciigrid('~/Downloads/ISLSCP_II_POT_VEG_961/data/potential_veg_hd.asc')
rast <- raster::raster(pv)
}

df <- as.data.frame(raster::rasterToPoints(rast))

if(FALSE){
  plot(df[,1],df[,2],col=df[,3])
}

lat.long.reg <- cbind(as.numeric(as.character(df$x)),as.numeric(as.character(df$y)))
lat.long.reg.df = data.frame(lat.long.reg)
colnames(lat.long.reg.df) = c('x', 'y')

coordinates(lat.long.reg.df) <- ~ x + y
proj4string(lat.long.reg.df) <- CRS('+proj=longlat +ellps=WGS84')

albers <- spTransform(lat.long.reg.df, CRS('+init=epsg:3175'))
albers <- as.matrix(data.frame(albers))

idx_centers = vector(length=nrow(inputdata3))

for(i in 1:length(idx_centers)){   
  center = inputdata3[i,1:2]
  d = rdist(matrix(center, ncol=2), as.matrix(albers))
  idx_centers[i] = which.min(d) 
}

table(df[idx_centers,3])
pv_names <- rep(NA,15)

pv_names[c(4,5,8,9,10)] <- c(
  'Temp NE EG F/W',
  'Bore EV F/W',
  'Everg/Decid Mix',
  'Savanna',
  'Grassland'
)

if(FALSE){
  pv_names[c(1, 5, 9, 10, 11)] <- c(
    'No Data',
    'Temp. NL EG F/W',
    'Mixed Forest',
    'Savanna',
    'Grassland/Steppe'
  )
}

if(FALSE){
  pv_keep <- data.frame(x=albers[idx_centers, 1],
                        y=albers[idx_centers, 2],
                        biomass = inputdata3$tot_biomass,
                        pv_use=pv_names[df[idx_centers,3]+1])
}

pv_keep <- data.frame(x=albers[idx_centers, 1],
                      y=albers[idx_centers, 2],
                      biomass = inputdata3$tot_biomass,
                      pv_use=pv_names[df[idx_centers,3]])


pdf('RF_map.pdf')
if(FALSE) pv_colors <- c('white',tim.colors(16))
pv_colors <- c(tim.colors(15))
plot(albers[,1],albers[,2],col=pv_colors[(df[,3])],
     cex=.5, pch=19,
     xlim=range(albers[idx_centers, 1])+c(-400000,400000),
     ylim=range(albers[idx_centers, 2])+c(-400000,400000))
plot(usShp, add=T, lwd=2)
points(albers[idx_centers, 1],
     albers[idx_centers, 2],
     cex = .1,
     col = 'black')
legend('topright',legend=pv_names,pch=19,col=pv_colors)
dev.off()

PV <-
  ggplot() + theme_bw() + theme(legend.text=element_text(size=15),
                                legend.title=element_text(size=17),
                                title = element_text(size=15)) +
  labs(x = "Latitude (Albers Projection)", y = "Biomass (Mg/ha)", color = "Potential Vegetation") +
  geom_jitter(data = pv_keep,
              aes(x = y, y = biomass, color = pv_use),
              size = 4)+
  scale_color_manual(values = pv_colors[c(4,5,8,9,10)])#+ stat_smooth(data = temp.precip.bio, aes(x = y, y = total), method = "lm",
print(PV)
ggsave(PV, file = 'RF_v_biomass.pdf')


library(egg)
A_PV <- ggarrange(A,PV, ncol = 1)
ggsave(A_PV, file = 'RF_v_biomass_together.pdf')

