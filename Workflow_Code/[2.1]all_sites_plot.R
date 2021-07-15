
library(maps)

load('~/ReFAB/nimble_pull2018-10-31.Rdata') #all sites in neotoma
load('cast.x.Rdata') #our calibration dataset
load('prediction.data_v5.Rdata') #our prediction dataset
dataID <- read.csv('dataID_v5.csv')


x.meta.plot <- x.meta[which(x.meta$site.name %in% dataID$name),]
x.meta.plot <- x.meta.plot[-which(x.meta.plot$site.name == 'Lily Lake'),]
x.meta.plot <- x.meta.plot[-which(x.meta.plot$site.name == 'Mud Lake'),]

pdf('[1]all-sites.pdf')
map('state', ylim=range(pol_cal_count$lat)+c(4, 2), xlim=range(pol_cal_count$long)+c(-1, 1),main=NA)
points(pol_cal_count$long, pol_cal_count$lat, pch=19, cex=1,col="gray")
points(cast.x$long, cast.x$lat, pch=4, cex=1,col="red")
points(x.meta.plot$long, x.meta.plot$lat, pch=1, cex=1,col="blue")

legend('topright',c('Available in Neotoma on 10-31-2018','Calibration','Prediction'),
       col = c('gray','red','blue'),pch = c(19,4,1))
dev.off()



