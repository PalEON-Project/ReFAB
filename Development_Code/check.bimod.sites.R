library(analogue)

get.sites <- read.csv(file='~/babystepps/site.num.bimodal.csv')
lat.long.keep <- matrix(NA,nrow = nrow(get.sites),ncol=2)
sample.bimod <- matrix(NA,nrow = nrow(get.sites),ncol=21)
row.nums <- list()

age.bins <- seq(0,10000,100)

x.meta$age.index <- as.matrix(as.numeric(cut(x.meta$age_bacon,breaks = age.bins,labels=seq(1:(length(age.bins)-1)))))

for(i in 1:nrow(get.sites)){
  lat.long.keep[i,] <- as.numeric(x.meta[x.meta$site.name==as.character(get.sites[i,1]),c('lat','long')][1,])
  sample.bimod[i,] <- ten.count[x.meta$site.name==as.character(get.sites[i,1])&x.meta$age.index==get.sites[i,2],]
  }

colnames(sample.bimod) <- colnames(ten.count)

x <- get.sites[,2]
breaks <-  seq(0,100,10)
colors <- rev(rainbow(length(breaks)-1,start=0,end=.7))

data_binned <-  cut(x, c(breaks), include.lowest = FALSE, labels = FALSE)

map('state', xlim=range(lat.long.keep[,2])+c(-2, 2), ylim=range(lat.long.keep[,1])+c(-1, 1))
points(lat.long.keep[,2], lat.long.keep[,1], pch=19, cex=1, col = colors[data_binned])
legend('bottomright',legend=signif(breaks[1:(length(breaks)-1)],digits=2),pch=rep(19,length(breaks)-1),col=colors)
title('Sites and Times with Bimod Liks')

load("~/ReFAB/Data/calibration.data.Rdata")
calib.data <- prop.table(Y,margin = 1)
testing.data <- prop.table(sample.bimod,margin=1) #prop.table(ten.count,margin=1)#
analog1 <- analog(x=calib.data,y=testing.data,'euclidean')
sum.analog1 <- summary(analog1)
sum.analog1

pdf('Distribution_dissim_all.pdf')
plot(dissim(analog1))
dev.off()

plot(minDC(analog1),depth=1:nrow(testing.data),type="p")
title('Dissimilarity of Bimodal Samples')

x <- as.numeric(unlist(minDC(analog1),use.names = FALSE))[1:21]
breaks <- c(0,as.numeric(unlist(minDC(analog1),use.names = FALSE))[23:26],1)
colors <- rev(rainbow(length(breaks)-1,start=0,end=.7))

data_binned <-  cut(x, c(breaks), include.lowest = FALSE, labels = FALSE)

map('state', xlim=range(lat.long.keep[,2])+c(-2, 2), ylim=range(lat.long.keep[,1])+c(-1, 1))
points(lat.long.keep[,2], lat.long.keep[,1], pch=19, cex=1, col = colors[data_binned])
legend('bottomright',legend=signif(breaks[1:(length(breaks)-1)],digits=2),pch=rep(19,length(breaks)-1),col=colors)
title('Euclidean distance from calibration')


pdf('Biplots.pdf')
par(mfrow=c(2,2))
all.rda <- princomp(prop.table(ten.count,margin = 1))
biplot(all.rda)
title('All Samples')

calib.rda <- princomp(calib.data[,])
biplot(calib.rda)
title('Calibration Samples')

test.rda <- princomp(testing.data[1:20,1:20])
biplot(test.rda)
title('Bimodal')
dev.off()

quartz()
pdf('boxplots.bimodal.pdf')
par(mfrow=c(2,1))
boxplot(testing.data,main='Bimodal Samples')
boxplot(calib.data,main='Calibration Dataset')
dev.off()
