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
  
  glm_mod = glm(cbind(counts[,i],total_counts-counts[,i]) ~ Z - 1,family=binomial(link="logit"))   
  points(biomass,counts[,i]/total_counts,pch=19,cex=.4,col='grey')
  new.biomass = seq(1,400,1)
  Z.new = bs(new.biomass,intercept=TRUE)
  lines(new.biomass, predict(glm_mod,newdata=list(Z=Z.new),type="response"),col="blue")  
  points(biomass,exp(Z%*%beta.est.real)[,i]/rowSums(exp(Z%*%beta.est.real)),col="red")
}

#post processing
biomass.preds = cbind(new.site.locs[-61,],summary(csamp.real.pred)$statistics[,1])
colnames(biomass.preds)<-c("x","y","biomass")

biomass_gam_mod = gam(log(biomass) ~ s(x,y,k=100),data = as.data.frame(biomass.preds))
load("/Users/paleolab/Documents/babySTEPPS/biomass_dat5.Rdata")
biomass_dat_est <- read.csv("biomass_estimate_v1.9.csv")
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

thresh = 800
values <- seq(0,800,1)

#for(i in 2:n.iter){
colnames(biomass.preds)<-c("x","y","biomass")
full.mat <- cbind(biomass_dat5[,1:2],xiao_ests,pred_biomass_gam)
colnames(full.mat) <- c("x","y","X.biomass","biomass")
y = as.data.frame(full.mat) #rowSums(biomass_dat_est) to make xiaopings
#id = which(y[,3] > thresh,arr.ind=T)
#y[id,3] = thresh

biomass_dat6 = melt(y, c('x','y'))

d <- ggplot() + geom_raster(data=biomass_dat6, aes(x=x,y=y,fill=value)) +
  scale_fill_gradientn(colours=tim.colors(length(values)),values=values, 
                       rescaler = function(x, ...) x, oob = identity) + 
  coord_fixed()

d <- d + facet_wrap(~ variable, ncol=1)

d <- theme_clean(d) + theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5))) +
  ggtitle("Xiaoping's Biomass")

d <- add_map_albers(plot_obj = d,  map_data = usFortified, dat = biomass_dat6)
print(d)




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











