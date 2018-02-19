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
hem.beech <- rowSums(biomass_dat_est[,c('Hemlock','Beech')])

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

### Biomass
full.mat <- cbind(biomass_dat_est[,c('x','y')],xiao_ests,hem.beech)
colnames(full.mat) <- c("x","y","Tot_Biomass",'late_succ')
tot.biom.df = as.data.frame(full.mat)

### Density
density <- read.csv('~/Downloads/plss_density_alb_v0.9-10.csv')
density.df = data.frame(x = density$x, y = density$y, 
           density = rowSums(density[,4:ncol(density)]))

y <- merge(tot.biom.df,density.df,by=c('x','y'))

breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)-1))

density.breaks <- c(0,1,48,1300)

legendName <- "Biomass (Mg/ha)"

data_binned <-  cut(y$Tot_Biomass, breaks, include.lowest = TRUE, labels = FALSE)
data_binned_density <-  cut(y$density, density.breaks, include.lowest = TRUE, labels = FALSE)

breaklabels <- apply(cbind(breaks[1:(length(breaks)-1)],
                           breaks[2:length(breaks)]), 1,
                     function(r) { sprintf("%0.2f - %0.2f", 
                                           r[1], r[2]) })

inputData <- data.frame(X = y$x, Y = y$y,
                        Preds = data_binned,tot_biomass = y$Tot_Biomass,
                        density = y$density, LS = y$late_succ)
inputData_long <- melt(inputData, c('X', 'Y'))

load('calibration.albers.rdata')
input_points <- data.frame(albers) # how to add points part A
colnames(input_points) <- c('lat','lon')

inputdata1 <- inputData[inputData$density>47,] #ask kelly about this number
colnames(inputdata1)[1:2] <- c('a','b')
col.plot <- rep('Forest',nrow(inputdata1))
data_binned1 <-  cut(inputdata1$tot_biomass, breaks, include.lowest = TRUE, labels = FALSE)

inputdata2 <- inputData[inputData$LS>1,] #what should the cutoff be?
colnames(inputdata2)[1:2] <- c('c','d')
col.plot2 <- rep('Late Successional',nrow(inputdata1))
data_binned2 <-  cut(inputdata2$tot_biomass, breaks, include.lowest = TRUE, labels = FALSE)



inputdata3 <- inputData[inputData$X>430000&inputData$X<557000&inputData$Y>678000&inputData$Y<1200000,]
colnames(inputdata3)[1:2] <- c('e','f')

d <-   ggplot() +
  geom_tile(data = inputData, aes(x = X, y = Y,fill = factor(data_binned))) +
  scale_fill_manual(labels = breaklabels, name = legendName, drop = FALSE, values = colors, guide = "legend") +
  geom_tile(data = inputdata1, aes(x = a, y = b, colour = 'Forest'), color='black',alpha=0, size=1) +
    #scale_color_manual(values='black', name=element_blank())
  geom_tile(data = inputdata1, aes(x = a, y = b, fill = factor(data_binned1))) + #puts biomass colors on
  geom_tile(data = inputdata2, aes(x = c, y = d, colour = col.plot2), color = 'blue', alpha=0, size=1) +
    #scale_color_manual(values='blue', name=element_blank())
  geom_tile(data = inputdata2, aes(x = c, y = d, fill = factor(data_binned2))) + #puts biomass colors on
  geom_tile(data = inputdata3, aes(x = e, y = f, color = 'Transect'),
              fill = 'gray', color = NA, size = .25, alpha = .5) +
  geom_point(data = input_points, aes(x=lat,y=lon), pch=16, size=1, alpha = .9 ,colour="brown") 
  
  
  #+ # how to add points part B
 # geom_text(data = input_points,aes(x=lat,y=lon,label=unlist(name.keep))) +
 # ggtitle("Xiaoping Estimates")

add_map_albers <- function(plot_obj, map_data = usFortified, dat){
  p <- plot_obj + geom_path(data = map_data, aes(x = long, y = lat, group = group), size = 0.1) +
    scale_x_continuous(limits = c(min(inputData$X, na.rm = TRUE), max(inputData$X, na.rm = TRUE))) +
    scale_y_continuous(limits = c(min(inputData$Y, na.rm = TRUE), max(inputData$Y, na.rm = TRUE)))
  return(p)
}

d <- add_map_albers(plot_obj = d, map_data = usFortified, dat = inputData_long)

quartz()
print(d)

ggsave(d,filename = paste0(Sys.Date(),'biomass_map.pdf'))



temp <- read.csv('~/Downloads/reclimatedata/tmean_yr_Prism_1900_full.csv')
precip <- read.csv('~/Downloads/reclimatedata/pr_monthly_Prism_1900_full.csv')

temp.precip <- merge(x = temp, y = precip, by = c('x','y'))
colnames(inputdata3)[1:2] <- c('x','y')
temp.precip.bio <- merge(temp.precip,inputdata3,by=c('x','y'))

#temp.precip.bio$Mean*temp.precip.bio$total #smooth metric

LS_binned <- rep('Prairie/Savanna',length(temp.precip.bio$LS))
LS_binned[temp.precip.bio$density>47] <- 'Forested (No Hemlock)'
LS_binned[temp.precip.bio$LS>1] <- 'Hemlock'

A <- ggplot()+ labs(x = '',y = "Biomass (Mg/ha)", color = "Ecosystem Type")+geom_jitter(data = temp.precip.bio,aes(x=y,y=tot_biomass,color=LS_binned))

B <- ggplot()+labs(x = '', y = "Mean Temperature [ºC]", color = "Ecosystem Type")+ geom_jitter(data = temp.precip.bio, aes(x=y,y=Mean,color=LS_binned))#+ stat_smooth(data = temp.precip.bio, aes(x = y, y = Mean), method = "lm",
                                                                                              # formula = y ~ poly(x, 10), se = FALSE) #mean temp
C <- ggplot()+ labs(x = "Latitude (Albers Projection)", y = "Total Precip [cm]", color = "Ecosystem Type") +geom_jitter(data = temp.precip.bio, aes(x=y,y=total,color=LS_binned))#+ stat_smooth(data = temp.precip.bio, aes(x = y, y = total), method = "lm",
                                                                              # formula = y ~ poly(x, 10), se = FALSE) #total precip 

library(egg)
F<-ggarrange(A,B,C, ncol = 1)
ggsave(F,file='transect_output_scatters.pdf')



