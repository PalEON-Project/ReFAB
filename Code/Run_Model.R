
####
#### Calibration Models
####

z.seq1 = seq(1,400, length=5)
Z.knots = bs(biomass,intercept=TRUE,df=4)

beta = matrix(NA,ncol(Z.knots),ncol(Y))
p = matrix(NA,nrow(Y),ncol(Y)) ; phi = p
J = nrow(Y)

library("rjags")

#data.sim.cal = list("R" = ncol(Z), "I" = ncol(Y), "J" = nrow(Y), "Y" = Y , "n" = rowSums(Y),"Z" = Z, "beta" = beta, "p" = p)
data.real.cal = list("R" = ncol(Z.knots), "I" = ncol(Y), "J" = nrow(Y), "Y" = counts , "n" = rowSums(counts),"Z" =  Z.knots, "beta" = beta, "p" = p)

inits.cal = list(list(beta = matrix(0,ncol(Z.knots),ncol(Y))),list(beta = matrix(3,ncol(Z.knots),ncol(Y))))

#mod.sim.cal <- jags.model(paste0(model.dir,'biomass_jags1.R'),data = data.sim.cal, n.chains = length(inits.cal), n.adapt = n.adapt,inits = inits.cal)
#csamp.sim.cal <- coda.samples(mod.sim.cal,c("beta"),n.iter = n.iter)
#csamp.p.sim.cal <- coda.samples(mod.sim.cal,c("p"),n.iter = n.iter)

mod.real.cal <- jags.model(paste0(model.dir,'biomass_jags1.R'),data = data.real.cal, n.chains = length(inits.cal), n.adapt = n.adapt, inits = inits.cal)
csamp.real.cal <- coda.samples(mod.real.cal,c("beta"),n.iter = n.iter)
#csamp.p.real.cal <- coda.samples(mod.real.cal,c("p"),n.iter = n.iter)

#load("pollen_lik_explore.Rdata")

#if(DRAW == TRUE) pdf(paste0(dump.dir,"calib.pdf"))
#par(mfrow=c(1,2))
#plot(c(betas),summary(csamp.sim.cal)$statistics[,1],main = "Simulated Data",xlab = "GLM betas",ylab="beta estimates",xlim=c(-9,5),ylim=c(-9,5))
#points(c(betas),summary(csamp.sim.cal)$quantiles[,2],main = "Real Data", pch=16,col="blue")
#lines(seq(-10,8,1),seq(-10,8,1))

plot(c(betas),summary(csamp.real.cal)$statistics[,1],main = "Real Data",xlab = "GLM betas",ylab="beta estimates",xlim=c(-9,5),ylim=c(-9,5))
points(c(betas),summary(csamp.real.cal)$quantiles[,2],pch=16,col="blue")
lines(seq(-10,8,1),seq(-10,8,1))
#if(DRAW == TRUE) dev.off()

print("finished cal1")

####
#### Prediction Models
####

J = nrow(Y)#nrow(Z.new)
Zb = matrix(NA,J,ncol(Z.knots))
DFS = ncol(Zb)
p = matrix(NA,J,ncol(Y)); phi.first = p; phi = p
#beta.est.sim = matrix(summary(csamp.sim.cal)$statistics[,1],4,ncol(Y))
beta.est.real = matrix(summary(csamp.real.cal)$statistics[,1],ncol(Z.knots),ncol(Y))
new.biomass = seq(1,400,1)
Z.new = bs(new.biomass,intercept=TRUE,knots = attr(Z.knots,"knots"),Boundary.knots = attr(Z.knots,"Boundary.knots"))

#data.sim.pred = list("I" = ncol(Y), "DFS" = DFS, "J" = J, "Zb" = Zb, "Y" = Y , "n" = rowSums(Y), "Z" = Z.new, "beta" = beta.est.sim, "p" = p, "phi.first" = phi.first)
data.real.pred = list("I" = ncol(Y), "DFS" = DFS, "J" = J, "Zb" = Zb, "Y" = counts , "n" = rowSums(counts), "Z" = Z.new, "beta" = beta.est.real, "p" = p, "phi.first" = phi.first)

inits.pred = list(list(b = rep(20,J)),list(b=rep(300,J)))

#mod.sim.pred <- jags.model(paste0(model.dir,'biomass_pred_jags.R'),data = data.sim.pred, n.chains = length(inits.pred), n.adapt = n.adapt, inits = inits.pred)
#csamp.sim.pred <- coda.samples(mod.sim.pred,c("b"),n.iter=n.iter)
#csamp.sim.pred.p <- coda.samples(mod.sim.pred,c("p"),n.iter=n.iter)

mod.real.pred <- jags.model(paste0(model.dir,'biomass_pred_jags.R'),data = data.real.pred, n.chains = length(inits.pred), n.adapt = n.adapt, inits = inits.pred)
csamp.real.pred <- coda.samples(mod.real.pred,c("b"),n.iter=n.iter)
#csamp.real.pred.p <- coda.samples(mod.real.pred,c("p"),n.iter=n.iter)

#plot(csamp.sim.pred)
#plot(csamp.real.pred)

#if(DRAW == TRUE) pdf(paste0(dump.dir,"pred1.pdf"))

#par(mfrow=c(1,2))
#plot(biomass,summary(csamp.sim.pred)$statistics[,1],main = "Simulated Data", xlab = "true biomass",ylab="biomass estimates",ylim=c(0,400),xlim=c(0,400))
#points(biomass,summary(csamp.sim.pred)$quantiles[,2],pch=16,col="blue")
#lines(seq(0,400,1),seq(0,400,1))

pine.df4 = summary(csamp.real.pred)

plot(biomass,summary(csamp.real.pred)$statistics[,1],main = "Real Data -- 4 df -- no pine", xlab = "true biomass",ylab="biomass estimates",ylim=c(0,400),xlim=c(0,400))
points(biomass,summary(csamp.real.pred)$quantiles[,2],pch=16,col="blue")
lines(seq(0,400,1),seq(0,400,1))

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





