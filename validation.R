ddirchmulti <- nimbleFunction(
  run = function(x = double(1), alpha = double(1), size = double(0), log_value = integer(0)){
    returnType(double(0))
    logProb <- lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum(lgamma(alpha + x)) - lgamma(sum(alpha) + size)
    
    if(log_value) {
      return(logProb)
    } else {
      return(exp(logProb))
    }
    
  }
)

# set up the "r" function
rdirchmulti <- nimbleFunction(
  run = function(n = integer(0), alpha = double(1), size = double(0)) {
    returnType(double(1))
    if(n != 1) nimPrint("rdirchmulti only allows n = 1; using n = 1.")
    p <- rdirch(1, alpha)
    return(rmulti(1, size = size, prob = p))
  })

# tell NIMBLE about the newly available distribution
registerDistributions(list(ddirchmulti = list(BUGSdist = "ddirchmulti(alpha, size)",
                                              types = c('value = double(1)', 'alpha = double(1)'))))

pred_code <- nimbleCode({
  for(j in 1:J){
    b[j] ~ dunif(0,145)
    
    Zb[j,1:5] <- bs_nimble(b[j], u[1:3], N0[1:2], N1[1:3], N2[1:4], N3[1:5])
  }
  
  for(i in 1:I){
    for(j in 1:J){
      phi.first[j,i] <- sum(Zb[j,1:5] %*% beta[1:5,i])
      phi.first1[j,i] <- sum(Zb[j,1:5] %*% beta1[1:5,i])
    }
  }
  
  for(j in 1:J){
    for(i in 1:I){
      exp.phi[j,i] <- exp(phi.first[j,i])
      exp.phi1[j,i] <- exp(phi.first1[j,i])
    }
  }
  
  for(j in 1:J){
    p.true[j,1] ~ dbeta(exp.phi[j,1],exp.phi1[j,1])
    p.rel[j,1] <- p.true[j,1]
    
    for(i in 2:(I-1)){
      p.rel[j,i]  ~ dbeta(exp.phi[j,i],exp.phi1[j,i]) 
      p.true[j,i] <-  p.rel[j,i] * (1 - sum(p.true[j,1:(i-1)]))
    }	
    p.true[j,21] <- 1 - sum(p.true[j,1:20])
  }    
  
  for(j in 1:J){
    Y[j,] ~ dmulti(p.true[j,],n[j])
  }
  
})

J = nrow(counts)#nrow(Z.new)
Zb = matrix(NA,J,ncol(Z.knots))
DFS = ncol(Zb)
phi = matrix(NA,J,ncol(counts)); phi.first = phi;
#beta.est.sim = matrix(summary(csamp.sim.cal)$statistics[,1],4,ncol(Y))
#load("beta.samps.Rdata")
beta.est = matrix(colMeans(samples.mixed[100:nrow(samples.mixed),]),ncol(Z.knots),ncol(Y)) #matrix(summary(csamp.real.cal)$statistics[,1],ncol(Z.knots),ncol(Y))# 
new.biomass = seq(1,u[length(u)],1)
#Z.new = bs(new.biomass,intercept=TRUE,knots = attr(Z.knots,"knots"),Boundary.knots = attr(Z.knots,"Boundary.knots"))
Z.new = matrix(0,nrow=length(new.biomass),ncol=5)
u <- u #should have defined knots in calibration
u<-c(rep(attr(Z.knots,"Boundary.knots")[1],1),attr(Z.knots,"knots"),rep(attr(Z.knots,"Boundary.knots")[2],1))

for(i in 1:length(new.biomass)){
  u_given <- new.biomass[i]
  Z.new[i,] = bs_nimble(u_given, u=u, N0 = rep(0, (length(u)-1)),N1 = rep(0, (length(u))), N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))
}

data.pred = list(Y = as.matrix(counts))

constants.pred = list(beta = beta1.est.real, beta1 = beta2.est.real, I = ncol(counts),
                      DFS = DFS, J = J, n = rowSums(counts), u = u,
                      N0 = rep(0, (length(u)-1)),N1 = rep(0, (length(u))),
                      N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))

inits.pred = list(b=rep(25,J))

dimensions.pred = list(exp.phi = dim(phi), exp.phi1 = dim(phi), phi.first = dim(phi),
                       phi.first1 = dim(phi), Zb = dim(Zb), beta = dim(beta.est),
                       beta1 = dim(beta.est), Y = dim(counts), p.rel = dim(phi), 
                       p.true = dim(phi))

save(pred_code, inits.pred, constants.pred, data.pred, dimensions.pred, file = 'validation.Rdata')

model_pred <- nimbleModel(pred_code, inits = inits.pred, constants = constants.pred,
                          data = data.pred, dimensions = dimensions.pred)

#Cmodel.pred <- compileNimble(model_pred) #opens browser...
spec.pred <- configureMCMC(model_pred, thin = 10, print = TRUE, useConjugacy = FALSE,
                           control = list(log=TRUE))
spec.pred$addMonitors(c('b')) 
Rmcmc.pred <- buildMCMC(spec.pred)
cm <- compileNimble(model_pred)
Cmcmc.pred <- compileNimble(Rmcmc.pred, project = model_pred) #Error in cModel$.nodeValPointers_byGID : $ operator not defined for this S4 class

ptm <- proc.time()
Cmcmc.pred$run(5000)
samples.pred <- as.matrix(Cmcmc.pred$mvSamples)
proc.time() - ptm

pdf('trace.settlement.pdf')
par(mfrow=c(2,2))
for(i in 1:ncol(samples.pred)){
  plot(samples.pred[,i],ylim=c(0,145),typ='l')
  abline(h=biomass[i],col='red')
}
dev.off()


#97,77,76,3

# #samples.pred1<-samples.pred
#load("samples.pred1.Rdata")
save(samples.pred,file="twothirds.pred_horizon_plus.Rdata")

pdf(paste0("pred_validation_two_betas_horizon_plus",Sys.Date(),".pdf"))

par(mfrow=c(1,1))
plot(biomass,
     colMeans(samples.pred[100:nrow(samples.pred),
                           grep('b',colnames(samples.pred))]),
     xlim=c(0,145),ylim=c(0,145),pch=19,xlab="True Biomass",ylab="Predicted Mean Biomass")
abline(a=0,b=1)
abline(lm(biomass~colMeans(samples.pred[100:nrow(samples.pred),grep('b',colnames(samples.pred))])+0),lty=2)
mtext(paste("r-squared",summary(lm(biomass~colMeans(samples.pred[100:nrow(samples.pred),grep('b',colnames(samples.pred))])+0))$r.squared))

arrows(x0 = biomass, y0 = apply(samples.pred,2,FUN = quantile,.05),
       x1 = biomass,y1 = apply(samples.pred,2,FUN = quantile,.975),
       code = 0, lwd=2)

dev.off()
