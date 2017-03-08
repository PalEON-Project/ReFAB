biomass_dat_est <- read.csv(paste0('~/Downloads/',"biomass_prediction_v0.9-10_bam.csv"))

b = biomass_dat_est
tot_b = b[,ncol(b)]

hist(rowSums(b[,4:25]) - b$Total)

table <- rbind(colMeans(prop.table(as.matrix(b[tot_b>=0 & tot_b<10,5:25]),margin = 1)),
               colMeans(prop.table(as.matrix(b[tot_b>=10 & tot_b<20,5:25]),margin = 1)),
               colMeans(prop.table(as.matrix(b[tot_b>=20 & tot_b<30,5:25]),margin = 1)),
               colMeans(prop.table(as.matrix(b[tot_b>=30 & tot_b<40,5:25]),margin = 1)),
               colMeans(prop.table(as.matrix(b[tot_b>=40 & tot_b<50,5:25]),margin = 1)),
               colMeans(prop.table(as.matrix(b[tot_b>=50 & tot_b<75,5:25]),margin = 1)),
               colMeans(prop.table(as.matrix(b[tot_b>=75 & tot_b<100,5:25]),margin = 1)),
               colMeans(prop.table(as.matrix(b[tot_b>=100 & tot_b<125,5:25]),margin = 1)),
               colMeans(prop.table(as.matrix(b[tot_b>=125 & tot_b<150,5:25]),margin = 1)),
               colMeans(prop.table(as.matrix(b[tot_b>=150,5:25]),margin = 1)))

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


bio.df<-melt(biomass_dat_est[,4:25])
#now make your lovely plot
ggplot(bio.df, aes(length, fill = variable)) + geom_density(alpha = 0.2)

ggplot(vegLengths, aes(length, fill = veg)) + geom_density(alpha = 0.2)
#You can say something about biomass because you can see they are different
#Is it gradual? or sudden changes in biomass?
#pca or squared cord dist #metric of changing composition over biomass
#hist of each spp over biomass


plot.which <- rev(names(sort(colMeans(b[,4:25]))))

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
breaks <-  seq(0,200,10)#c(seq(0,50,10),seq(75,200,25))
data_binned <-  cut(biomass_dat_est$Total, c(breaks), include.lowest = FALSE, labels = FALSE)
b<-cbind(biomass_dat_est,data_binned)
b.melt <- melt(data = b)

c <- matrix(NA,max(data_binned),21)

for(i in 1:max(data_binned)){
  c[i,] <- colSums(b[b$data_binned==i,c(5:25)])
}

colnames(c) <- colnames(b[,c(5:25)])

c1 <- prop.table(c,1)
c2 <- c1[,order(colMeans(c1))]
c1.melt<-melt(c2)


c2.melt <- within(c1.melt, X2 <- factor(X2,names(sort(colMeans(c1),decreasing = FALSE))))

pdf('settlement.biomass.tiles.pdf')
ggplot() + geom_tile(data = c2.melt, aes(x = X1, y = X2, fill = value), colour = 'grey') +
  scale_fill_gradient(name = 'Total Bio. Prop.',low="white",high="black") +
  ylab('Tree Taxa') +
  xlab('Biomass Categories Mg/ha') +
  scale_x_continuous(breaks=seq(.5,( max(data_binned) +.5), 1), labels = breaks[1:(max(data_binned)+1)])
dev.off()


load('add.bacon2.Rdata')
data_binned <-  cut(biomass, c(breaks), include.lowest = FALSE, labels = FALSE)
b<-as.data.frame(cbind(Y,data_binned))

c <- matrix(NA,max(data_binned),(ncol(Y)-1))

for(i in 1:max(data_binned)){
  c[i,] <- colSums(b[b$data_binned==i,1:(ncol(Y)-1)])
}

colnames(c) <- colnames(b[,1:(ncol(Y)-1)])

c1 <- prop.table(c,1)
c2 <- c1[,order(colMeans(c1))]
c1.melt<-melt(c2)


c2.melt <- within(c1.melt, X2 <- factor(X2,names(sort(colMeans(c1),decreasing = FALSE))))

pdf('settlement.pollen.tiles.pdf')
ggplot() + geom_tile(data = c2.melt, aes(x = X1, y = X2, fill = value), colour = 'grey') +
  scale_fill_gradient(name = 'Count Prop.',low="white",high="black") +
  ylab('Tree Taxa') +
  xlab('Biomass Categories Mg/ha') +
  scale_x_continuous(breaks=seq(.5,( max(data_binned) +.5), 1), labels = breaks[1:(max(data_binned)+1)])
dev.off()






