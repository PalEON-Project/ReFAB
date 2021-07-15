
library(ncdf4)

pls_biomass <- readRDS('~/Downloads/PLS_aboveground_biomass_by_taxa_and_climate.rds')

breaks <- seq(0,15,1)
colors <- rainbow(length(breaks)-1)
data_binned <- cut(pls_biomass$MAT,breaks)
plot(pls_biomass$x,pls_biomass$y,col=colors[data_binned])
data_binned_fia <- cut(fia_biomass$MAT,breaks)
plot(fia_biomass$x,fia_biomass$y,col=colors[data_binned_fia])

fia_biomass <- readRDS('~/Downloads/FIA_aboveground_biomass_by_taxa_and_climate.rds')

just_temp <- pls_biomass[,c('x','y','MAT')]

pca_dens <- read.csv('~/Downloads/density_full_unc_v1.0.csv')
head(pca_dens)
nc <- nc_open(file.path('~','Downloads','PLS_biomass_agb_western_point_v1.0rc1.nc'))

x <- nc$dim$x$vals
y <- nc$dim$y$vals
data <- ncvar_get(nc,varid = c('Total'))

rownames(data) <- x
colnames(data) <- y

r1 <- raster(list(x=x,y=y,z=data))

#can do
#plot(r1)

biomass_dat_est <- as.data.frame(rasterToPoints(r1))
colnames(biomass_dat_est) <- c('x','y','pls')
head(biomass_dat_est)


nc <- nc_open(file.path('~','Downloads','FIA_biomass_agb_point_v1.0rc1.nc'))
x <- nc$dim$x$vals
y <- nc$dim$y$vals
data <- ncvar_get(nc,varid = c('Total'))

rownames(data) <- x
colnames(data) <- y

r2 <- raster(list(x=x,y=y,z=data))

biomass_dat_est_fia <- as.data.frame(rasterToPoints(r2))
colnames(biomass_dat_est_fia) <- c('x','y','fia')

r3 <-
  merge(biomass_dat_est,
       biomass_dat_est_fia, by = c('x','y'))
plot(r3$pls, r3$fia)

r4_1 <- 
  merge(r3,pca_dens,by = c('x','y'))

r5 <- 
  merge(r4_1,just_temp,by = c('x','y'))


r4 <- r5[sample(x = 1:nrow(r5),size = 5000,replace = F),]

pdf('pls_v_fia_plot.pdf')
layout(matrix(c(1,2,3,1,2,3,4,4,4),3,3))
par(mar=c(4,4,1,0),cex.lab=1.25)
#Temperature
plot(r4$pls, r4$MAT, ylab = 'Mean Annual Temp (C)',
     xlab = NA,cex=.1)
points(r4$fia, r4$modtmean, col = 'red',cex=.1)

#Precip
plot(r4$pls,
     r4$MAP1910,
     ylim = range(c(r4$MAP1910, r4$MAP2011), na.rm = T),
     ylab = 'Mean Annual Precip (cm)',
     xlab = NA,cex=.1,pch=19)
points(r4$fia, r4$MAP2011, col = 'red',cex=.1,pch=19)

#PCA
plot(r4$pls, r4$PC1,
     ylab = 'PCA 1',
     xlab = NA,cex=.1,pch=19)
points(r4$fia, r4$PC1fia, col = 'red',cex=.1,pch=19)
mtext(text = 'Biomass (Mg/ha)',side = 1,outer = F,line=2.5)

plot.new()
legend('center',c('PLS','FIA'),pch=19,cex=2,col=c('black','red'))
dev.off()
