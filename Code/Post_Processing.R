sum.real.pred = summary(csamp.real.pred)$statistics[,1]

round(rbind(colMeans(counts[biomass<160&sum.real.pred>200,]/rowSums(counts[biomass<160&sum.real.pred>200,])),
            colMeans(counts[biomass<100&biomass<160&sum.real.pred<100,]/rowSums(counts[biomass<100&biomass<160&sum.real.pred<100,]))),digits = 2)

par(mfrow=c(1,1))
map('state', xlim=range(plot_biomass_pollen[,2])+c(-2, 2), ylim=range(plot_biomass_pollen[,3])+c(-1, 1))
pine.prop = counts[,4]/rowSums(counts)
points(plot_biomass_pollen[biomass<150&biomass>50&pine.prop<.2,2], plot_biomass_pollen[biomass<150&biomass>50&pine.prop<.2,3], pch=19, cex=1)

par(mfrow=c(3,3))
for(i in 1:ncol(counts)){
  gam_mod = gam(cbind(counts[,i],total_counts-counts[,i]) ~ s(biomass),family=binomial(link="logit"))
  plot(biomass,counts[,i]/total_counts,pch=19,cex=.4,col='grey',ylab="Pollen Prop",main=colnames(counts)[i])
  points(biomass,predict(gam_mod,type="response"),pch=19,col="green")
  
  glm_mod = glm(cbind(counts[,i],total_counts-counts[,i]) ~ Z.knots - 1,family=binomial(link="logit"))   
  points(biomass,counts[,i]/total_counts,pch=19,cex=.4,col='grey')
  new.biomass = seq(0,400,1)
  #Z.new = bs(new.biomass,intercept=TRUE)
  #Z.new = bs(new.biomass,intercept=TRUE,knots = attr(Z,"knots"),Boundary.knots = attr(Z,"Boundary.knots"))
  #lines(new.biomass, predict(glm_mod,newdata=list(Z=Z.new),type="response"),col="blue")  
  points(biomass,exp(Z.knots%*%beta.est.real)[,i]/rowSums(exp(Z.knots%*%beta.est.real)),col="red")
}

#post processing
new.site.locs = new.site.locs[-c(61,67),]
biomass.preds = cbind(new.site.locs[-sites_rm,],summary(csamp.real.pred)$statistics[,1])
colnames(biomass.preds)<-c("x","y","biomass") #important step #gam has problems with formatting otherwise

biomass_gam_mod = gam(log(biomass) ~ s(x,y,k=80),data = as.data.frame(biomass.preds))
load("/Users/paleolab/Documents/babySTEPPS/biomass_dat5.Rdata")
biomass_dat_est <- read.csv(paste0(data.dir,"biomass_estimate_v1.9.csv"))
xiao_ests <- rowSums(biomass_dat_est)
pred_biomass_gam = exp(predict(biomass_gam_mod,newdata = as.data.frame(biomass_dat5[,1:2])))

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

thresh = 400
values <- c(seq(0,100,25),seq(100,200,50),300,400)

#for(i in 2:n.iter){
colnames(biomass.preds)<-c("x","y","biomass")
xiao_ests1 = numeric(length(xiao_ests))
for(i in 1:length(xiao_ests)){
  if(xiao_ests[i] > 400) xiao_ests1[i] = 400 else xiao_ests1[i] = xiao_ests[i]
}

for(i in 1:length(pred_biomass_gam)){
  if(pred_biomass_gam[i] > 400) pred_biomass_gam[i] = 400 else pred_biomass_gam[i] = pred_biomass_gam[i]
}
full.mat <- cbind(biomass_dat5[,1:2],xiao_ests1,as.vector(pred_biomass_gam))
colnames(full.mat) <- c("x","y","Xiao biomass","pred biomass")
y = as.data.frame(full.mat) #rowSums(biomass_dat_est) to make xiaopings
#id = which(y[,3] > thresh,arr.ind=T)
#y[id,3] = thresh

biomass_dat6 = melt(y, c('x','y'))

d <- ggplot() + geom_raster(data=biomass_dat6, aes(x=x,y=y,fill=value)) +
  scale_fill_gradientn(colours=sort(terrain.colors(length(values)),decreasing=TRUE),values=values, 
                       rescaler = function(x, ...) x, oob = identity) + 
  coord_fixed()

d <- d + facet_wrap(~ variable, ncol=1)

d <- theme_clean(d) + theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5))) +
  ggtitle("Biomass")

d <- add_map_albers(plot_obj = d,  map_data = usFortified, dat = biomass_dat6)

quartz()
print(d)

# chris's code for plotting discrete values
# data would be a vector of biomas values

breaks <-  c(0, 25, 50, 75, 100, 125, 150, 200, 250, 300, 400)

colors <- rev(terrain.colors(length(breaks)-1))

breaks <-  c( -400, -300, -200, -100, 100, 200, 300, 400)

colors <- c("dark blue", "blue", "light blue", "white", "pink", "red")
legendName <- "Difference"

data_binned <-  cut(y[,3] - y[,4], breaks, include.lowest = TRUE, labels = FALSE)

breaklabels <- apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1,  function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })

inputData <- data.frame(X = y[,1], Y = y[,2], Preds = cbind(data_binned,data_binned))
inputData_long <- melt(inputData, c('X', 'Y'))

colnames(new.site.locs) <- c('lat','lon')
input_points <- data.frame(new.site.locs[-sites_rm,])

d <- ggplot() + geom_raster(data = inputData_long, aes(x = X, y = Y, fill = factor(value))) + scale_fill_manual(labels = breaklabels, name = legendName, drop = FALSE, values = colors, guide = "legend") + 
  theme(strip.text.x = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16)) +
  geom_point(data = input_points, aes(x=lat,y=lon), pch=16, size=2,colour="black") +
  ggtitle("Xiao Est - Max List Pred")

add_map_albers <- function(plot_obj, map_data = usFortified, dat){
  p <- plot_obj + geom_path(data = map_data, aes(x = long, y = lat, group = group), size = 0.1) +
    scale_x_continuous(limits = c(min(inputData$X, na.rm = TRUE), max(inputData$X, na.rm = TRUE))) +
    scale_y_continuous(limits = c(min(inputData$Y, na.rm = TRUE), max(inputData$Y, na.rm = TRUE)))
  return(p)
}

d <- add_map_albers(plot_obj = d, map_data = usFortified, dat = inputData_long)

quartz()
print(d)

#remove_sites1 <- plot_biomass_pollen[-c(61,67),]
breaks <-  c(0, 25, 50, 75, 100, 125, 150, 200, 250, 300, 400)
colors <- rev(terrain.colors(length(breaks)-1))
data_binned1 <-  cut(biomass.preds[,3], breaks, include.lowest = TRUE, labels = colors)

quartz()
map('state', xlim=range(plot_biomass_pollen[,2])+c(-2, 2), ylim=range(plot_biomass_pollen[,3])+c(-1, 1))
points(remove_sites1[-sites_rm,2], remove_sites1[-sites_rm,3], pch=19, cex=1, col = as.character(data_binned1))
title("Min List Point Preds")


head(y)
y[y[,1] == input_points[,1],]

breaks <-  c( -400, -300, -200, -100, 100, 200, 300, 400)
colors <- c("dark blue", "blue", "light blue", "white", "pink", "red")
data_binned2 <-  cut(biomass.preds[,3], breaks, include.lowest = TRUE, labels = colors)

quartz()
map('state', xlim=range(plot_biomass_pollen[,2])+c(-2, 2), ylim=range(plot_biomass_pollen[,3])+c(-1, 1))
points(remove_sites1[-sites_rm,2], remove_sites1[-sites_rm,3], pch=19, cex=1, col = as.character(data_binned1))
title("Min List Point Preds")




multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
# 
# 
# 
# 
# 
# 
# 
# 
# blank = exp(Z.new%*%betas)/rowSums(exp(Z.new%*%betas))
# blank[1:30,]
# 
# colnames(blank)<-colnames(counts)
# colMeans(Y[biomass<20,])
# colMeans(Y[biomass>20&biomass<50,])
# 
# Z = bs(biomass,knots=9)
# glm_mod = glm(cbind(counts[,i],total_counts-counts[,i]) ~ Z - 1,family=binomial(link="logit")) 
# new.biomass = seq(1,600,1)
# Z.new = bs(new.biomass)
# colnames(Z.new) <- c("Z1","Z2","Z3","Z4")
# plot(new.biomass,predict(glm_mod,newdata=list(Z=Z.new),type="response"))
# 
# 
# #dbinom(counts[,i], prob = inv.logit(Z.new%*%glm_mod$coefficients), size = total_counts, log = TRUE))
# 
# 
# save_biom = array(0,dim=c(600,ncol(counts),142))
# #par(mfrow=c(3,3))
# for(p in 1:1){
# 
# for(i in 1:ncol(counts)){
#   glm_mod = glm(cbind(counts[,i],total_counts-counts[,i]) ~ Z - 1,family=binomial(link="logit")) 
# 	save_biom[,i,p] <- dbinom(counts[p,i],prob = inv.logit(Z.new%*%glm_mod$coefficients),size = total_counts[p],log = T)	
# }
#     plot(new.biomass,rowSums(save_biom[,,p]))
# 	
# }











