library(nimble)
library(splines)
library(maps)
library(analogue)

setwd("/Users/paleolab/babySTEPPS/")

load("min.list.june24.Rdata")
load(file="/Users/paleolab/babySTEPPS/Data/pol.cal.count.mnwi1.csv") 

training = prop.table(as.matrix(ten.count),1)
rownames(training)<-cast.x$SiteID

x = pol.cal.count[pol.cal.count$Age>=200,]
x = x[x$Age<=10000,]

x.meta = x[,1:6]
x = x[,7:ncol(x)]

x = x[,-which(colnames(x)==c("PINUSX"))]
trees <- c("ALNUSX","JUGLANSX","ACERX","CUPRESSA","FRAXINUX","FAGUS","CYPERACE","LARIXPSEU","TSUGAX","QUERCUS","TILIA","BETULA","PICEAX","OSTRYCAR","ULMUS","ABIES","POPULUS")
other.trees <- c("TAXUS","NYSSA","CASTANEA","PLATANUS","SALIX","LIQUIDAM")
ten.count = matrix(0,nrow(x),length(trees)+3)
prairie <- c("ARTEMISIA","ASTERX","POACEAE","AMBROSIA","CHENOAMX","CORYLUS")
ten.count[,1] <- as.numeric(rowSums(x[,prairie]))
ten.count[,2] <- as.numeric(rowSums(x[,other.trees]))
ten.count[,3:(length(trees)+2)] <- as.matrix(x[,trees])
ten.count[,(length(trees)+3)] <- as.numeric(rowSums(x)) - as.numeric(rowSums(ten.count))
colnames(ten.count)<-c("prairie","other trees",trees,"other herbs")

ten.count.save = ten.count
ten.count = round(ten.count.save)

testing.save = prop.table(as.matrix(ten.count),1)
rownames(testing.save)<-x.meta$Age
testing.meta <- x.meta
testing <- testing.save

analog1 <- analog(x=training, y=testing,'SQchord')
sum.analog1 <- summary(analog1)
sum.analog1

plot(dissim(analog1))

plot(minDC(analog1),depth=as.numeric(rownames(testing)),type="p")

cma.analog1<-cma(analog1,cutoff=.3)
plot(cma.analog1)

pdf("dissim.maps.pdf")
par(mfrow=c(2,2))
age.vec <- seq(1000,10000,500)
plot.leg <- rep(c(FALSE,FALSE,TRUE),length.out=length(age.vec))
for(i in 2:length(age.vec)){
map('state', xlim=c(-98,-81), ylim=c(41,50))
plot.age.low <- age.vec[i-1]
plot.age.high <- age.vec[i]
plot.pts <- x.meta[x.meta$Age>=plot.age.low&x.meta$Age<=plot.age.high,]
dissim.values <- minDC(analog1)$minDC[which(x.meta$Age>=plot.age.low&x.meta$Age<=plot.age.high)]
	breaks <-  seq(0,1,.1)
	colors <- rainbow(length(breaks),start = 0, end = .75)
	data_binned <-  cut(dissim.values, breaks,
	 include.lowest = FALSE, labels = FALSE)
points(plot.pts[,3], plot.pts[,2], pch=19,
		cex=1,col = colors[data_binned])
title(paste(plot.age.low,"-",plot.age.high,"Before Present"))
cuts<-levels(cut(dissim.values
,breaks = breaks))
cuts<-gsub(","," - ",cuts)
cuts<-gsub("\\(","[",cuts)
cuts
if(plot.leg[i] == TRUE){
	plot.new()
	legend("center",legend = cuts,col=colors,pch=19)
}	
}
dev.off()
















