library(nimble)
library(splines)
library(maps)
library(analogue)

load("threethirds_v2.0.Rdata")

training = prop.table(as.matrix(Y),1)
rownames(training)<-1:nrow(Y)

load('prediction.data_v4.Rdata')

testing.save = prop.table(as.matrix(ten.count),1)
rownames(testing.save)<-x.meta$age_bacon
testing.meta <- x.meta
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




dissim.values <- minDC(analog1)$minDC
breaks <-  seq(0, 1, .1)
colors <- rev(tim.colors(length(breaks)))#rainbow(length(breaks), start = 0, end = .75)
data_binned <-  cut(dissim.values,
                    breaks,
                    include.lowest = FALSE,
                    labels = FALSE)
pdf('no_analog_assessment.pdf')
par(mfrow=c(1,1))
plot(
  x.meta$age_bacon,
  x.meta$lat.x,
  bg = colors[data_binned],
  pch = 21,
  cex = 1.5,
  xlab = 'Sample Age',
  ylab = 'Latitude',
  main = 'Pollen Dissimilarity'
)
abline(v=10000,lwd=2)
legend("bottomright",
       legend = cuts,
       col = colors,
       pch = 19,
       title = 'minDC value')
dev.off()








