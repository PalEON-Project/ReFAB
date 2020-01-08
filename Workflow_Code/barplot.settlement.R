

#### This script is for plotting figure 4 A and B in the supplements

library(ncdf4)
library(raster)
library(reshape)
library(ggplot2)

nc <- nc_open(file.path('~','Downloads','PLS_biomass_agb_western_point_v1.0rc1.nc'))

x <- nc$dim$x$vals
y <- nc$dim$y$vals
vars <- names(nc$var)
data <- list()
biomass_mat <- matrix(NA,11784,length(vars))
for(i in 1:length(vars)){
  data[[i]] <- ncvar_get(nc,varid = vars[[i]])
  r1 <- raster(list(x=x,y=y,z=data[[i]]))
  biomass_mat[,i] <- rasterToPoints(r1)[,3]
}
colnames(biomass_mat) <- vars

b = biomass_mat#biomass_dat_est
tot_b = b[,'Total']

#### biomass_dat_get from get_data.R
load('biomass_dat_est_agb.Rdata')

breaks <-  c(seq(0, 100, 25), seq(150, 250, 50))
data_binned <-
  cut(biomass_dat_est$layer,
      c(breaks),
      include.lowest = FALSE,
      labels = FALSE)
b <- as.data.frame(cbind(biomass_mat, data_binned))
b.melt <- melt(data = b)

c <- matrix(NA, max(data_binned, na.rm = T), 22)

for (i in 1:max(data_binned, na.rm = T)) {
  c[i, ] <- colSums(b[b$data_binned == i, -1], na.rm = T)
}
colnames(c) <- colnames(b[, -1])

c <- c[, -22]


c1 <- prop.table(c, 1)
c2 <- c1[, order(colMeans(c1))]
c1.melt <- melt(c2)

c2.melt <-
  within(c1.melt, X2 <-
           factor(X2, names(sort(
             colMeans(c1, na.rm = T), decreasing = FALSE
           )))) #has to be reshape not reshape2

pdf('[5-A]settlement.biomass.tiles.pdf')
 ppp <- ggplot() + geom_tile(data = c2.melt,
                     aes(x = X1, y = X2, fill = value),
                     colour = 'grey') +
  
  scale_fill_gradient(name = 'Proportion of \n Total Biomass', low = "white", high = "black") +
  ylab('Tree Taxa') +
  xlab('Biomass Categories (Mg/ha)') +
  scale_x_continuous(breaks = seq(.5, (max(data_binned, na.rm = T) + .5), 1), 
                     labels = breaks[1:(max(data_binned, na.rm = T) + 1)])
 ppp <- ppp + theme(legend.text=element_text(size=12),
             legend.title=element_text(size=10),
             title = element_text(size=16),
             axis.text = element_text(size = 10))
 ppp
dev.off()



load('threethirds_v3.0.Rdata')
breaks <-  c(seq(0, 100, 25), seq(150, 450, 50))
data_binned <-  cut(biomass, c(breaks), include.lowest = FALSE, labels = FALSE)
Y.prop <- prop.table(as.matrix(Y),1)
colnames(Y.prop) <- c('Pine','Prairie','Oak','Birch','Other Herb.','Sedge',
                      'Alder','Ironwood','Elm','Hemlock',
                      'Spruce','Maple',
                      'Ash','Poplar','Cypress','Other Arboreal',
                      'Tamarack','Beech','Hickory','Fir','Basswood',
                      'Walnut')
b<-as.data.frame(cbind(Y.prop,data_binned))

c <- matrix(NA,max(data_binned,na.rm = T),(ncol(b)-1))

for(i in 1:max(data_binned,na.rm = T)){
  c[i,] <- colSums(b[b$data_binned==i,1:(ncol(b)-1)],na.rm=T)
}

colnames(c) <- colnames(b[,1:(ncol(b)-1)])

c1 <- prop.table(c,1)

#c <- c[-which(rowSums(c)==0),]
c1[is.na(c1)] <- 0
c2 <- c1[,order(colMeans(c1))]
c1.melt<-melt(c2)

c2.melt <- within(c1.melt, X2 <- factor(X2,names(sort(colMeans(c1),decreasing = FALSE))))

pdf('[5-B]settlement.pollen.tiles.pdf')
ggplot() + geom_tile(data = c2.melt, aes(x = X1, y = X2, fill = value), colour = 'grey') +
  scale_fill_gradient(name = 'Pollen \n Proportion',low="white",high="black") +
  ylab('Tree Taxa') +
  xlab('Biomass Categories Mg/ha') +
  scale_x_continuous(breaks=seq(.5,( max(data_binned,na.rm = T) +.5), 1), labels = breaks[1:(max(data_binned,na.rm = T)+1)])+ 
  theme(legend.text=element_text(size=12),
        legend.title=element_text(size=10),
        title = element_text(size=16),
        axis.text = element_text(size = 10))
dev.off()



################################ End figures begin exploration notes


#biomass_dat_est <- read.csv(paste0('~/Downloads/',"biomass_prediction_v0.9-10_bam.csv"))

nc <- nc_open(file.path('~','Downloads','PLS_biomass_agb_western_point_v1.0rc1.nc'))

x <- nc$dim$x$vals
y <- nc$dim$y$vals
vars <- names(nc$var)
data <- list()
biomass_mat <- matrix(NA,11784,length(vars))
for(i in 1:length(vars)){
  data[[i]] <- ncvar_get(nc,varid = vars[[i]])
  r1 <- raster(list(x=x,y=y,z=data[[i]]))
  biomass_mat[,i] <- rasterToPoints(r1)[,3]
}
colnames(biomass_mat) <- vars

b = biomass_mat#biomass_dat_est
tot_b = b[,'Total']

par(mfrow=c(1,1))
plot(rowSums(b[,-which(colnames(b)%in%'Total')]), tot_b,
     ylim=c(0,400),xlim=c(0,400),col='darkgray')
abline(b=1,a=0,lwd=2)

table <- rbind(colMeans(prop.table(as.matrix(b[tot_b>=0 & tot_b<10,-1]),margin = 1)),
               colMeans(prop.table(as.matrix(b[tot_b>=10 & tot_b<20,-1]),margin = 1)),
               colMeans(prop.table(as.matrix(b[tot_b>=20 & tot_b<30,-1]),margin = 1)),
               colMeans(prop.table(as.matrix(b[tot_b>=30 & tot_b<40,-1]),margin = 1)),
               colMeans(prop.table(as.matrix(b[tot_b>=40 & tot_b<50,-1]),margin = 1)),
               colMeans(prop.table(as.matrix(b[tot_b>=50 & tot_b<75,-1]),margin = 1)),
               colMeans(prop.table(as.matrix(b[tot_b>=75 & tot_b<100,-1]),margin = 1)),
               colMeans(prop.table(as.matrix(b[tot_b>=100 & tot_b<125,-1]),margin = 1)),
               colMeans(prop.table(as.matrix(b[tot_b>=125 & tot_b<150,-1]),margin = 1)),
               colMeans(prop.table(as.matrix(b[tot_b>=150,-1]),margin = 1)))

bluefunc <- colorRampPalette(c('red','orange','yellow','green','blue','purple'))
bluefuncs <- bluefunc(21)

quartz()
pdf('barplot.settlement.comp.pdf')
barplot(t(table[,rev(order(colMeans(table)))[1:10]]),col=bluefuncs,names.arg = c('0-10','10-20','20-30','30-40',
                                             '40-50','50-75','75-100','100-125',
                                             '125-150','>150'),cex.names = .7,
        ylab = 'species prop',xlab = 'biomass categories (Mg/ha)')
plot.new()
legend('left',colnames(b[,5:15]),col=bluefuncs[1:11],pch=rep(19,20),cex=2)
legend('right',colnames(b[,16:25]),col=bluefuncs[12:21],pch=rep(19,20),cex=2)
dev.off()

carrots <- data.frame(length = rnorm(100000, 6, 2))
cukes <- data.frame(length = rnorm(50000, 7, 2.5))

#Now, combine your two dataframes into one.  First make a new column in each.
carrots$veg <- 'carrot'
cukes$veg <- 'cuke'

#and combine into your new data frame vegLengths
vegLengths <- rbind(carrots, cukes)


bio.df<-melt(biomass_dat_est[,-1])
#now make your lovely plot
ggplot(bio.df, aes(length, fill = variable)) + geom_density(alpha = 0.2)

ggplot(vegLengths, aes(length, fill = veg)) + geom_density(alpha = 0.2)
#You can say something about biomass because you can see they are different
#Is it gradual? or sudden changes in biomass?
#pca or squared cord dist #metric of changing composition over biomass
#hist of each spp over biomass


plot.which <- rev(names(sort(colMeans(b[,-1]))))

pdf('pls.explor.scatters.pdf')
par(mfrow=c(3,3))
plot(tot_b,bio.df[bio.df$variable==plot.which[1],'value'],
     xlim = c(0,200),col=rainbow(length(plot.which),alpha=.5)[1],pch=19,
     main=plot.which[1],ylab='spp biomass prop',xlab='total biomass')
n=1
for(i in plot.which[2:length(plot.which)]){
  n = n + 1
  plot(tot_b,bio.df[bio.df$variable==i,'value'],
         col=rainbow(length(plot.which),alpha = .5)[n],pch=19,main=plot.which[n],
       ylab='spp biomass prop',xlab='total biomass')
}
#legend('right',as.character(plot.which),col=rainbow(length(plot.which)),pch=rep(19,22),cex=1)
dev.off()






