pdf("new.splines.pdf")
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
dev.off()

#post processing
biomass.preds = cbind(final_coors[,3:4],summary(csamp.real.pred)$statistics[,1])
colnames(biomass.preds)<-c("x","y","biomass") #important step #gam has problems with formatting otherwise

biomass_gam_mod = gam(log(biomass) ~ s(x,y,k=80),data = as.data.frame(biomass.preds))
#load("/Users/paleolab/Documents/babySTEPPS/biomass_dat5.Rdata")
biomass_dat_est <- read.csv(paste0(data.dir,"biomass_prediction_v0.2.csv"))
xiao_ests <- rowSums(biomass_dat_est[,4:23])
pred_biomass_gam = exp(predict(biomass_gam_mod,newdata = as.data.frame(biomass_dat_est[,1:2])))

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
#for(i in 2:n.iter){
colnames(biomass.preds)<-c("x","y","biomass")

full.mat <- cbind(biomass_dat_est[,1:2],xiao_ests,as.vector(pred_biomass_gam))
colnames(full.mat) <- c("x","y","Xiao Total Biomass","Smoothed Biomass")
y = as.data.frame(full.mat)

# chris's code for plotting discrete values
# data would be a vector of biomas values

#Regular breaks and colors
breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)-1))

#Difference breaks and colors
breaks <-  c(seq(-100,-25,25),-10,10,seq(25,100,25))
colors <- cm.colors(length(breaks)-1)

legendName <- "Biomass (Mg/ha)"

data_binned <-  cut(y[,3]-y[,4], breaks, include.lowest = TRUE, labels = FALSE)

breaklabels <- apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1,  function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })

inputData <- data.frame(X = y[,1], Y = y[,2], Preds = cbind(data_binned,data_binned))
inputData_long <- melt(inputData, c('X', 'Y'))

input_points <- data.frame(final_coors[,3:4])
colnames(input_points) <- c('lat','lon')

d <- ggplot() + geom_raster(data = inputData_long, aes(x = X, y = Y, fill = factor(value))) + scale_fill_manual(labels = breaklabels, name = legendName, drop = FALSE, values = colors, guide = "legend") + 
  theme(strip.text.x = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16)) +
  geom_point(data = input_points, aes(x=lat,y=lon), pch=16, size=2,colour="black") +
  ggtitle("Xiaoping estimates - smoothed prediciton")

add_map_albers <- function(plot_obj, map_data = usFortified, dat){
  p <- plot_obj + geom_path(data = map_data, aes(x = long, y = lat, group = group), size = 0.1) +
    scale_x_continuous(limits = c(min(inputData$X, na.rm = TRUE), max(inputData$X, na.rm = TRUE))) +
    scale_y_continuous(limits = c(min(inputData$Y, na.rm = TRUE), max(inputData$Y, na.rm = TRUE)))
  return(p)
}

d <- add_map_albers(plot_obj = d, map_data = usFortified, dat = inputData_long)

quartz()
print(d)

data_binned1 <-  cut(biomass - biomass.preds[,3], breaks, include.lowest = TRUE, labels = colors)

quartz()
map('state', xlim=range(final_coors[,2])+c(-2, 2), ylim=range(final_coors[,1])+c(-1, 1))
points(final_coors[,2],final_coors[,1], pch=19, cex=1, col = as.character(data_binned1))
title("Difference in Point Predictions for Subset of Settlement Points")
legend("topright",c("-100 - -75","-75 - -50","-50 - -25","-25 -10", "-10 - 10","10 - 25","25 - 50",
                    "50-75","75-100"),col = colors, pch = 19)
legend("topright",c("0-10","10-20","20-30","30-40","40-50","50-75","75-100",
                    "100-125","125-150","150-175"),col = colors, pch = 19)


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











