library(nimble)
library(splines)
library(maps)
library(analogue)
library(fields)

load("threethirds_v2.0.Rdata")

training = prop.table(as.matrix(Y),1)
rownames(training)<-1:nrow(Y)

load('prediction.data_v4.Rdata')
dataID <- read.csv('dataID_bacon_v4.csv')

ten.count.used <- ten.count[x.meta$site.name%in%unique(dataID$name),]
x.meta.used <- x.meta[x.meta$site.name%in%unique(dataID$name),]


testing.save = prop.table(as.matrix(ten.count.used),1)
rownames(testing.save)<-x.meta.used$age_bacon
testing.meta <- x.meta.used
testing <- testing.save

analog1 <- analog(x=training, y=testing,'SQchord')
sum.analog1 <- summary(analog1)
sum.analog1

plot(dissim(analog1))

plot(minDC(analog1),depth=as.numeric(rownames(testing)),type="p")

cma.analog1<-cma(analog1,cutoff=.3)
plot(cma.analog1)


age.vec <- seq(1000,10000,100)
plot.leg <- rep(c(FALSE,FALSE,TRUE),length.out=length(age.vec))

plot.leg <- rep(FALSE,length(age.vec))
plot.leg[2] <- TRUE

pdf("dissim.maps.pdf")

par(mfrow=c(4,5),mar=c(0,0,0,0),oma=c(0,0,0,0))
for(i in rev(2:length(age.vec))) {
  map('state', xlim = c(-98, -81), ylim = c(41, 50))
  plot.age.low <- age.vec[i - 1]
  plot.age.high <- age.vec[i]
  pick.pts <- which(x.meta$age_bacon >= plot.age.low &
                      x.meta$age_bacon <= plot.age.high)
  plot.pts <- x.meta[pick.pts, ]
  dissim.values <- minDC(analog1)$minDC[pick.pts]
  breaks <-  seq(0, 1, .1)
  colors <- rainbow(length(breaks), start = 0, end = .75)
  data_binned <-  cut(dissim.values,
                      breaks,
                      include.lowest = FALSE,
                      labels = FALSE)
  points(plot.pts[, 'long.x'],
         plot.pts[, 'lat.x'],
         pch = 19,
         cex = 1,
         col = colors[data_binned])
  title(paste(plot.age.low, "-", plot.age.high, "Before Present"))
  cuts <- levels(cut(dissim.values, breaks = breaks))
  cuts <- gsub(",", " - ", cuts)
  cuts <- gsub("\\(", "[", cuts)
  cuts
  if (plot.leg[i] == TRUE) {
    plot.new()
    legend("center",
           legend = cuts,
           col = colors,
           pch = 19)
  }
}

dev.off()

##### bigwoods
if(FALSE){
  bigwoods_calib <- c(2,10,11,13,14,18,20,26,29,32)
  
  testing.save_bw = prop.table(as.matrix(ten.count.used),1)
  rownames(testing.save_bw)<-x.meta.used$age_bacon
  testing.meta_bw <- x.meta.used
  testing_bw <- testing.save_bw
  analog1_bw <- analog(x=training_bw, y=testing_bw,'SQchord')
  mins_bw <- minDC(analog1_bw)$minDC
  breaks_bw <-  c(0, .05, 1)
  colors_bw <- c(adjustcolor('darkgray',alpha.f = .4))#adjustcolor(rev(tim.colors(length(breaks_bw))),alpha.f = .75)#rainbow(length(breaks), start = 0, end = .75)
  data_binned_bw <-  cut(mins_bw,
                         breaks_bw,
                         include.lowest = FALSE,
                         labels = FALSE)
}


dissim.values <- minDC(analog1)$minDC
breaks <-  seq(0, 1, .1)
colors <- adjustcolor(rev(tim.colors(length(breaks)-1)),alpha.f = .75)#rainbow(length(breaks), start = 0, end = .75)
data_binned <-  cut(dissim.values,
                    breaks,
                    include.lowest = FALSE,
                    labels = FALSE)
pdf('no_analog_assessment.pdf')
par(mfrow=c(1,1))
plot(
  x.meta.used$age_bacon,
  x.meta.used$lat.x,
  bg = colors[data_binned],
  pch = 21,
  col = c(adjustcolor('darkgray',alpha.f = .4)),
  cex = 2,
  xlab = 'Sample Age (years before present)',
  ylab = 'Latitude',
  main = 'Pollen Dissimilarity'
)
abline(v=10000,lwd=2)
legend("bottomright",
       legend = c(cuts),
       pt.bg = c(colors),
       col = c(rep(colors_bw[2],length(breaks)-1)),
       pch = 21,
       cex = 1,
       title = 'minDC value')
dev.off()


#####
##### Bigwoods
#####

load("twothirds_v2.0.Rdata")

breaks <-  c(seq(0,50,10),seq(75,250,25))
colors <- rev(terrain.colors(length(breaks)-1))
data_binned_biomass <-  cut(biomass, c(breaks), include.lowest = FALSE, labels = FALSE)


bigwoods_calib <- c(2,10,11,13,14,18,20,26,29,32)

map('state', xlim=c(-100,-80), ylim=c(39.5,49.5))
load('two.thirds.cast.x.Rdata')
points(ag.two.thirds.cast.x[bigwoods_calib,'long'],
       ag.two.thirds.cast.x[bigwoods_calib,'lat'])
text(ag.two.thirds.cast.x[,'long'],
     ag.two.thirds.cast.x[,'lat'],
     labels=1:nrow(Y),cex=.3,
     col=colors[data_binned_biomass])

training_bw = prop.table(as.matrix(Y[bigwoods_calib,]),1)
rownames(training_bw)<-1:nrow(training_bw)

load('prediction.data_v4.Rdata')

dataID <- read.csv('dataID_bacon_v4.csv')

ten.count.used <- ten.count[x.meta$site.name%in%unique(dataID$name),]
x.meta.used <- x.meta[x.meta$site.name%in%unique(dataID$name),]

testing.save_bw = prop.table(as.matrix(ten.count.used),1)
rownames(testing.save_bw)<-x.meta.used$age_bacon
testing.meta_bw <- x.meta.used
testing_bw <- testing.save_bw
analog1_bw <- analog(x=training_bw, y=testing_bw,'SQchord')


mins_bw <- minDC(analog1_bw)$minDC

order_mins_bw <- mins_bw[order(mins_bw)]

breaks_bw <-  seq(0, .25, .05)
colors_bw <- adjustcolor(rev(tim.colors(length(breaks_bw))),alpha.f = .75)#rainbow(length(breaks), start = 0, end = .75)
data_binned_bw <-  cut(order_mins_bw,
                    breaks,
                    include.lowest = FALSE,
                    labels = FALSE)

map('state', xlim=c(-100,-80), ylim=c(39.5,49.5))
points(x.meta.used[order(mins_bw),c('long.x','lat.x')],
       col=colors_bw[data_binned_bw],pch=19)

length(which(mins_bw<.25))/length(mins_bw)

sort(table(x.meta.used[which(mins_bw<.1),'site.name'])/table(x.meta.used[,'site.name']))

dist2bigwoods <- mins
save(dist2bigwoods,file='dist2bigwoods.Rdata')
