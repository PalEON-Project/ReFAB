
####
#### Calibration Models
####

Z.knots = bs(biomass,intercept=TRUE,df=4)

beta = matrix(NA,ncol(Z.knots),ncol(Y))
p = matrix(NA,nrow(Y),ncol(Y)) ; phi = p
J = nrow(Y)

data.real.cal = list("R" = ncol(Z.knots), "I" = ncol(Y), "J" = nrow(Y), "Y" = counts , 
                     "n" = rowSums(counts),"Z" =  Z.knots, "beta" = beta, "p" = p)
inits.cal = list(list(beta = matrix(0,ncol(Z.knots),ncol(Y))),list(beta = matrix(3,ncol(Z.knots),ncol(Y))),list(beta = matrix(-3,ncol(Z.knots),ncol(Y))))

mod.real.cal <- jags.model(paste0(model.dir,'biomass_jags1.R'),data = data.real.cal, 
                           n.chains = length(inits.cal), n.adapt = n.adapt, inits = inits.cal)
csamp.real.cal <- coda.samples(mod.real.cal,c("beta"),n.iter = n.iter)

quartz()
plot(c(betas.save),summary(csamp.real.cal)$statistics[,1],main = "Real Data",xlab = "GLM betas",ylab="beta estimates",xlim=c(-9,5),ylim=c(-9,5))
points(c(betas.save),summary(csamp.real.cal)$quantiles[,2],pch=16,col="blue")
lines(seq(-10,8,1),seq(-10,8,1))

print("finished cal1")

save(csamp.real.cal,file ="beta.samps.Rdata")

beta.est.real = matrix(summary(csamp.real.cal)$statistics[,1],ncol(Z.knots),ncol(Y))

plot.betas <- as.matrix(exp(Z.knots%*%beta.est.real)/rowSums(exp(Z.knots%*%beta.est.real)))

pdf(paste0(fig.dir,"scatter.plus.betas.pdf"))
par(mfrow=c(3,3))
for(i in 1:ncol(counts)){
	plot(biomass,counts[,i]/total_counts,pch=19,cex=.7,col='black',ylab="Pollen Proportions",main=colnames(counts)[i],xlab="Biomass")
	points(biomass,plot.betas[,i],col='red',pch=19,cex=1)
  }
dev.off()

pdf(paste0(fig.dir,"calib.trace.pdf"))
plot(csamp.real.cal)
dev.off()
####
#### Prediction Models
####

# x = pol.cal.count[pol.cal.count$Age>=100,]
# x = x[x$Age<=200,]
# 
# x = x[,-which(colnames(x)==c("PINUSX"))]
# trees <- c("ACERX","CUPRESSA","FRAXINUX","FAGUS","CYPERACE","LARIXPSEU","TSUGAX","QUERCUS","TILIA","BETULA","PICEAX","OSTRYCAR","ULMUS","ABIES","POPULUS")
# other.trees <- c("TAXUS","NYSSA","JUGLANSX","CASTANEA","PLATANUS","SALIX","LIQUIDAM","ALNUSX")
# ten.count = matrix(0,nrow(x),length(trees)+3)
# prairie <- c("CORYLUS","ARTEMISIA","ASTERX","POACEAE","AMBROSIA","CHENOAMX")
# ten.count[,1] <- rowSums(x[,prairie])
# ten.count[,2] <- rowSums(x[,other.trees])
# ten.count[,3:(length(trees)+2)] <- as.matrix(x[,trees])
# ten.count[,(length(trees)+3)] <- rowSums(x[,7:ncol(x)]) - rowSums(ten.count)
# colnames(ten.count)<-c("prairie","other trees",trees,"other herbs")
# 
# ten.count <- round(ten.count)
# counts <- ten.count[1,]

J = nrow(counts)#nrow(Z.new)
Zb = matrix(NA,J,ncol(Z.knots))
DFS = ncol(Zb)
p = matrix(NA,J,ncol(counts)); phi.first = p; phi = p
#beta.est.sim = matrix(summary(csamp.sim.cal)$statistics[,1],4,ncol(Y))
new.biomass = seq(1,156,1)
Z.new = bs(new.biomass,intercept=TRUE,knots = attr(Z.knots,"knots"),Boundary.knots = attr(Z.knots,"Boundary.knots"))

#data.sim.pred = list("I" = ncol(Y), "DFS" = DFS, "J" = J, "Zb" = Zb, "Y" = Y , "n" = rowSums(Y), "Z" = Z.new, "beta" = beta.est.sim, "p" = p, "phi.first" = phi.first)
data.real.pred = list("I" = ncol(counts), "DFS" = DFS, "J" = J, "Zb" = Zb, "Y" = counts , "n" = rowSums(counts), "Z" = Z.new, "beta" = beta.est.real, "p" = p, "phi.first" = phi.first)

inits.pred = list(list(b = rep(20,J)),list(b=rep(150,J)))

#mod.sim.pred <- jags.model(paste0(model.dir,'biomass_pred_jags.R'),data = data.sim.pred, n.chains = length(inits.pred), n.adapt = n.adapt, inits = inits.pred)
#csamp.sim.pred <- coda.samples(mod.sim.pred,c("b"),n.iter=n.iter)
#csamp.sim.pred.p <- coda.samples(mod.sim.pred,c("p"),n.iter=n.iter)

mod.real.pred <- jags.model(paste0(model.dir,'biomass_pred_jags.R'),data = data.real.pred, n.chains = length(inits.pred), n.adapt = n.adapt, inits = inits.pred)
csamp.real.pred <- coda.samples(mod.real.pred,c("b"),n.iter=n.iter)
#csamp.real.pred.p <- coda.samples(mod.real.pred,c("p"),n.iter=n.iter)

nopine.max.list1 = summary(csamp.real.pred)$statistics[,1]
pine.max.list1 = summary(pine.max.list)$statistics[,1]
nopine.min.list1 = summary(nopine.min.list)$statistics[,1]
pine.min.list1 = summary(pine.min.list)$statistics[,1]

save(nopine.min.list,file="nopine.min.list.Rdata")

csamp.real.pred = nopine.min.list

quartz()
par(mfrow=c(2,2))

plot(biomass,nopine.min.list1 ,main = "No Pine -- Max List", 
     xlab = "true biomass",ylab="biomass estimates",ylim=c(0,180),xlim=c(0,180),
     pch = 19, cex = 1.2)
mod = lm(biomass~nopine.max.list1)
abline(mod,lty = 2)
text(150,25,paste0("R^2 = ",round(summary(mod)$r.squared,digits = 4)))
#points(biomass[which(biomass<20&sum.real.pred>70)],summary(csamp.real.pred)$statistics[which(biomass<20&sum.real.pred>70),1],pch=16,col="blue")
#points(biomass[which(biomass<20&sum.real.pred<50)],summary(csamp.real.pred)$statistics[which(biomass<20&sum.real.pred<50),1],pch=16,col="green")
lines(seq(0,400,1),seq(0,400,1))

plot(biomass,pine.max.list1,main = "Pine -- Max List", 
     xlab = "true biomass",ylab="biomass estimates",ylim=c(0,180),xlim=c(0,180),
     pch = 19, cex = 1.2)
mod = lm(biomass~pine.max.list1)
abline(mod,lty = 2)
text(150,25,paste0("R^2 = ",round(summary(mod)$r.squared,digits = 4)))
#points(biomass[which(biomass<20&sum.real.pred>70)],summary(csamp.real.pred)$statistics[which(biomass<20&sum.real.pred>70),1],pch=16,col="blue")
#points(biomass[which(biomass<20&sum.real.pred<50)],summary(csamp.real.pred)$statistics[which(biomass<20&sum.real.pred<50),1],pch=16,col="green")
lines(seq(0,400,1),seq(0,400,1))

plot(biomass,nopine.min.list1,main = "No Pine -- Min List", 
     xlab = "true biomass",ylab="biomass estimates",ylim=c(0,180),xlim=c(0,180),
     pch = 19, cex = 1.2)
mod = lm(biomass~nopine.min.list1)
abline(mod,lty = 2)
text(150,25,paste0("R^2 = ",round(summary(mod)$r.squared,digits = 4)))
#points(biomass[which(biomass<20&sum.real.pred>70)],summary(csamp.real.pred)$statistics[which(biomass<20&sum.real.pred>70),1],pch=16,col="blue")
#points(biomass[which(biomass<20&sum.real.pred<50)],summary(csamp.real.pred)$statistics[which(biomass<20&sum.real.pred<50),1],pch=16,col="green")
lines(seq(0,400,1),seq(0,400,1))

plot(biomass,pine.min.list1,main = "Pine -- Min List", 
     xlab = "true biomass",ylab="biomass estimates",ylim=c(0,180),xlim=c(0,180),
     pch = 19, cex = 1.2)
mod = lm(biomass~pine.min.list1)
abline(mod,lty = 2)
text(150,25,paste0("R^2 = ",round(summary(mod)$r.squared,digits = 4)))
#points(biomass[which(biomass<20&sum.real.pred>70)],summary(csamp.real.pred)$statistics[which(biomass<20&sum.real.pred>70),1],pch=16,col="blue")
#points(biomass[which(biomass<20&sum.real.pred<50)],summary(csamp.real.pred)$statistics[which(biomass<20&sum.real.pred<50),1],pch=16,col="green")
lines(seq(0,400,1),seq(0,400,1))

#plot(csamp.sim.pred)

pdf(paste0(dump.dir,"biom_posts.pdf"))
plot(csamp.real.pred)
dev.off()

#if(DRAW == TRUE) pdf(paste0(dump.dir,"pred1.pdf"))

#par(mfrow=c(1,2))
quartz()
par(mfrow=c(2,1))

sum.real.pred = summary(csamp.real.pred)$statistics[,1]

map('state', xlim=range(plot_biomass_pollen[,3])+c(-2, 2), ylim=range(plot_biomass_pollen[,2])+c(-1, 1))
pine.prop = counts[,4]/rowSums(counts)
points(plot_biomass_pollen[which(biomass<20&sum.real.pred>70),3], plot_biomass_pollen[which(biomass<20&sum.real.pred>70),2], pch=19, cex=1,col = "blue")
points(plot_biomass_pollen[which(biomass<20&sum.real.pred<50),3], plot_biomass_pollen[which(biomass<20&sum.real.pred<50),2], pch=19, cex=1,col = "green")

quartz()
plot(biomass,summary(csamp.real.pred)$statistics[,1],main = "Including Alnusx and Juglansx", 
     xlab = "true biomass",ylab="biomass estimates",ylim=c(0,180),xlim=c(0,180),
     pch = 19, cex = 1.2)
mod = lm(biomass~summary(csamp.real.pred)$statistics[,1])
abline(mod,lty = 2)
text(150,25,paste0("R^2 = ",round(summary(mod)$r.squared,digits = 4)))
#points(biomass[which(biomass<20&sum.real.pred>70)],summary(csamp.real.pred)$statistics[which(biomass<20&sum.real.pred>70),1],pch=16,col="blue")
#points(biomass[which(biomass<20&sum.real.pred<50)],summary(csamp.real.pred)$statistics[which(biomass<20&sum.real.pred<50),1],pch=16,col="green")
lines(seq(0,400,1),seq(0,400,1))


props2 = counts/rowSums(counts)

sort(colMeans(props[which(biomass<20&sum.real.pred>70),])-
       colMeans(props[which(biomass<20&sum.real.pred<50),]))

look = rbind(colMeans(props2[which(biomass<20&sum.real.pred>70),]),
         colMeans(props2[which(biomass<20&sum.real.pred<50),]))

rownames(look)<- c("over_est","on_line")

look[,order(look[1,])]

all.preds = summary(csamp.real.pred)

pdf("min.max.lists.pdf")
quartz()
par(mfrow=c(2,2))
plot(biomass,min.list$statistics[,1],main = "Min List", xlab = "true biomass",ylab="biomass estimates",ylim=c(0,400),xlim=c(0,400),pch=16,col="blue")
points(biomass,min.list$quantiles[,2])
lines(seq(0,400,1),seq(0,400,1))
legend("bottomright",c('Mean','Median'),col=c("blue","black"),pch=c(16,1))

plot(biomass,max.list$statistics[,1],main = "Max List", xlab = "true biomass",ylab="biomass estimates",ylim=c(0,400),xlim=c(0,400),pch=16,col="red")
points(biomass,max.list$quantiles[,2])
lines(seq(0,400,1),seq(0,400,1))
legend("bottomright",c('Mean','Median'),col=c("red","black"),pch=c(16,1))

plot(biomass,max.list$statistics[,1],main = "Max and Min List", xlab = "true biomass",ylab="biomass estimates",ylim=c(0,400),xlim=c(0,400),col='red',pch=16)
points(biomass,min.list$statistics[,1],pch=16,col="blue")
lines(seq(0,400,1),seq(0,400,1))
legend("bottomright",c('Max','Min'),col=c("red","blue"),pch=c(16,16))

dev.off()
#abline(h=200)
#abline(v=100)
#abline(h=100)
#if(DRAW == TRUE) dev.off()

#save.image(paste0(dump.dir,"mod_output.Rdata"))

print("finished pred1")

#map.x = plot_biomass_pollen[-sites_rm,2]
#map.y = plot_biomass_pollen[-sites_rm,3]
#biom.ests = summary(csamp.real.pred)$quantiles[,2]
#map('state', xlim=range(plot_biomass_pollen[-sites_rm,2])+c(-2, 2), ylim=range(plot_biomass_pollen[-sites_rm,3])+c(-1, 1))
#points(map.x[biomass>100&biom.ests<100], map.y[biomass>100&biom.ests<100], pch=19, cex=1)





