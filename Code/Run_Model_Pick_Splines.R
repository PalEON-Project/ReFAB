####
#### Large Basis Models
####

####
#### Calibration Models
####

#z.seq = seq(0,400, length=8)
Z.big.gen = bs(biomass, intercept = TRUE,df = 10) # sim with df=10 but model should still look through all 30
delta = 50#bigger
phi.b = matrix(0,nrow(counts),10); p = phi.b
Y = matrix(0,nrow(counts),10)
betas.big = matrix(0,10,10)

for(i in 1:ncol(counts)){
  glm_mod = glm(cbind(counts[,i],total_counts-counts[,i]) ~ Z.big.gen - 1, family=binomial(link="logit"))
  betas.big[,i] = glm_mod$coefficients
}


phi.b = exp(Z.big.gen%*%betas.big)/rowSums(exp(Z.big.gen%*%betas.big))

for(j in 1:nrow(counts)){
  p[j,] = rdirichlet(1,phi.b[j,]*delta)
  Y[j,] = rmultinom(1,prob = p[j,], size = rowSums(counts)[j])
}

beta.big = matrix(NA,30,10)
z.seq = seq(0,400, length=28)
Z.big = bs(biomass, intercept = TRUE, knots = z.seq[2:27]) # sim with df=10 but model should still look through all 30

data.sim.cal1 = list( "Y" = Y , "n" = rowSums(Y),"Z" = Z.big, "beta" = beta.big, "p" = p)
inits.cal1 = list(list(beta = matrix(0,nrow(beta.big),10)),list(beta = matrix(3,nrow(beta.big),10)))

mod.sim.cal1 <- jags.model(paste0(model.dir,'biomass_jags2.R'),data = data.sim.cal1, n.chains = length(inits.cal1), n.adapt = n.adapt,inits = inits.cal1)
csamp.sim.cal1.beta <- coda.samples(mod.sim.cal1,c("beta"),n.iter = n.iter)
csamp.sim.cal1.tau <- coda.samples(mod.sim.cal1,c("tau"),n.iter = n.iter)

beta.est.big = matrix(summary(csamp.sim.cal1.beta)$statistics[,1],30,10)

if(DRAW==TRUE) pdf(paste0(dump.dir,"splines2.pdf"))
par(mfrow=c(3,3))
for(i in 1:ncol(counts)){
  plot(biomass,Y[,i]/rowSums(Y),pch=19,cex=.4,col='grey',ylab="Pollen Prop",main=colnames(counts)[i])
  points(biomass,exp(Z.big.gen%*%betas.big)[,i]/rowSums(exp(Z.big.gen%*%betas.big)),col="blue")
  points(biomass,exp(Z.big%*%beta.est.big)[,i]/rowSums(exp(Z.big%*%beta.est.big)),col="red")
}
if(DRAW==TRUE) dev.off()

data.real.cal1 = list( "Y" = counts , "n" = rowSums(counts),"Z" = Z.big, "beta" = beta.big, "p" = p)
mod.real.cal1 <- jags.model(paste0(model.dir,'biomass_jags2.R'),data = data.real.cal1, n.chains = length(inits.cal1), n.adapt = n.adapt,inits = inits.cal1)
csamp.real.cal1.beta <- coda.samples(mod.real.cal1,c("beta"),n.iter = n.iter)

beta.est.big1 = matrix(summary(csamp.real.cal1.beta)$statistics[,1],30,10)

if(DRAW==TRUE) pdf(paste0(dump.dir,"splines3.pdf"))
par(mfrow=c(3,3))
for(i in 1:ncol(counts)){
  plot(biomass,counts[,i]/rowSums(counts),pch=19,cex=.4,col='grey',ylab="Pollen Prop",main=colnames(counts)[i])
  points(biomass,exp(Z.big%*%beta.est.big1)[,i]/rowSums(exp(Z.big%*%beta.est.big1)),col="red")
}
if(DRAW==TRUE) dev.off()

print("finished cal2")

####
#### Prediction Models
####

J = 141#nrow(Z.new)
Zb.big = matrix(NA,J,30)
p = matrix(NA,J,10); phi.first = p; phi = p
new.biomass = seq(1,400,1)
Z.new = bs(new.biomass,intercept=TRUE, knots = z.seq[2:27])

data.sim.pred1 = list("DFS" = ncol(Zb.big), "J" = J, "Zb" = Zb.big, "Y" = Y , "n" = rowSums(Y), "Z" = Z.new, "beta" = beta.est.big, "p" = p, "phi.first" = phi.first)

inits.pred = list(list(b = rep(20,J)),list(b=rep(300,J)))

mod.sim.pred1 <- jags.model(paste0(model.dir,'biomass_pred_jags.R'),data = data.sim.pred1, n.chains = length(inits.pred), n.adapt = n.adapt, inits = inits.pred)
csamp.sim.pred1 <- coda.samples(mod.sim.pred1,c("b"),n.iter=n.iter)

data.real.pred1 = list("DFS" = ncol(Zb.big), "J" = J, "Zb" = Zb.big, "Y" = counts , "n" = rowSums(counts), "Z" = Z.new, "beta" = beta.est.big1, "p" = p, "phi.first" = phi.first)

mod.real.pred1 <- jags.model(paste0(model.dir,'biomass_pred_jags.R'),data = data.real.pred1, n.chains = length(inits.pred), n.adapt = n.adapt, inits = inits.pred)
csamp.real.pred1 <- coda.samples(mod.real.pred1,c("b"),n.iter=n.iter)

if(DRAW==TRUE) pdf(paste0(dump.dir,"pred2.pdf"))
par(mfrow=c(1,1))
plot(biomass,summary(csamp.real.pred1)$statistics[,1],main = "Real Data", xlab = "true biomass",ylab="biomass estimates",ylim=c(0,400),xlim=c(0,400))
points(biomass,summary(csamp.real.pred1)$quantiles[,2],pch=16,col="blue")
lines(seq(0,400,1),seq(0,400,1))
if(DRAW==TRUE) dev.off()

print("finished pred2")

save.image(paste0(dump.dir,"mod_output.Rdata"))

print("Finished fitting models")