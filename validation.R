validation_model <- function(Y, Z.knots, samples.mixed, u, Niters,
                             bMax = 143, group_rm){
  
library(nimble)
source(file.path('genPareto','betabin.R')) # code for user-defined beta-binomial distribution
  
pred_code <- nimbleCode({
  for(j in 1:J){
    b[j] ~ dunif(0,bMax)
    
    Zb[j,1:5] <- bs_nimble(b[j], u[1:3], N0[1:2], N1[1:3], N2[1:4], N3[1:5])
  }
  
  for(i in 1:I){
    for(j in 1:J){
      shape1[j,i] <- exp(sum(Zb[j,1:5] %*% beta1[1:5,i]))
      shape2[j,i] <- exp(sum(Zb[j,1:5] %*% beta2[1:5,i]))
    }
  }
  
  for(j in 1:J){
    Y[j, 1] ~ dbetabin(shape1[j, 1], shape2[j, 1], n[j])
    for(i in 2:(I-1)){
      Y[j, i] ~ dbetabin(shape1[j, i], shape2[j, i], n[j] - sum(Y[j,1:(i-1)]))
      
    }
  }
  
  # for(j in 1:J){
  #   for(i in 1:I){
  #     exp.phi[j,i] <- exp(phi.first[j,i])
  #     exp.phi1[j,i] <- exp(phi.first1[j,i])
  #   }
  # }
  # 
  # for(j in 1:J){
  #   p.true[j,1] ~ dbeta(exp.phi[j,1],exp.phi1[j,1])
  #   p.rel[j,1] <- p.true[j,1]
  #   
  #   for(i in 2:(I-1)){
  #     p.rel[j,i]  ~ dbeta(exp.phi[j,i],exp.phi1[j,i]) 
  #     p.true[j,i] <-  p.rel[j,i] * (1 - sum(p.true[j,1:(i-1)]))
  #   }	
  #   p.true[j,I] <- 1 - sum(p.true[j,1:(I-1)])
  # }    
  # 
  # for(j in 1:J){
  #   Y[j,] ~ dmulti(p.true[j,],n[j])
  # }
  
})

i.beta1 <- grep("beta1",colnames(samples.mixed))
i.beta2 <- grep("beta2",colnames(samples.mixed))

burnin <- round(.2 * nrow(samples.mixed))

beta1.est.real = matrix(colMeans(samples.mixed[burnin:nrow(samples.mixed),i.beta1]),ncol(Z.knots),ncol(Y))
beta2.est.real = matrix(colMeans(samples.mixed[burnin:nrow(samples.mixed),i.beta2]),ncol(Z.knots),ncol(Y))

J = nrow(Y)
Zb = matrix(NA,J,ncol(Z.knots))
DFS = ncol(Zb)
phi = matrix(NA,J,ncol(Y)); phi.first = phi;

data.pred = list(Y = as.matrix(Y))

constants.pred = list(bMax = bMax, beta1 = beta1.est.real,
                      beta2 = beta2.est.real, I = ncol(Y),
                      DFS = DFS, J = J, n = rowSums(Y), u = u,
                      N0 = rep(0, (length(u)-1)), N1 = rep(0, (length(u))),
                      N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))

inits.pred = list(b=rep(25,J))

dimensions.pred = list(shape1 = dim(phi), shape2 = dim(phi),
                       Zb = dim(Zb), beta1 = dim(beta1.est.real),
                       beta2 = dim(beta2.est.real), Y = dim(Y))

#save(pred_code, inits.pred, constants.pred, data.pred, dimensions.pred, file = 'validation.Rdata')

model_pred <- nimbleModel(pred_code, inits = inits.pred, constants = constants.pred,
                          data = data.pred, dimensions = dimensions.pred)

#Cmodel.pred <- compileNimble(model_pred) #opens browser...
spec.pred <- configureMCMC(model_pred, thin = 10, print = TRUE, useConjugacy = FALSE,
                           control = list(log=TRUE))
spec.pred$addMonitors(c('b')) 
Rmcmc.pred <- buildMCMC(spec.pred)
cm <- compileNimble(model_pred)
Cmcmc.pred <- compileNimble(Rmcmc.pred, project = model_pred) 

vals <- 1:bMax
outLik = outPost = matrix(NA, bMax, J)
for(j in 1:J){
  calcNodes <-  cm$getDependencies(paste0('b[',j,']'))
  for(val in vals) {
    cm$b[j] <- val
    outPost[val,j] = calculate(cm,calcNodes)# cm$calculate(calcNodes)
    # likelihood portion
    outLik[val,j] =  calculate(cm,calcNodes[grep("Y", calcNodes)]) # cm$calculate(calcNodes[45])  #
  }	
}

save(outLik, file=paste0('outLik.group.',group_rm,'.Rdata'))

bInit <- apply(outLik,2,which.max)

inits_pred = list(b = bInit)
cm$setInits(inits_pred)

ptm <- proc.time()
set.seed(0)
Cmcmc.pred$run(Niters)
samples.pred <- as.matrix(Cmcmc.pred$mvSamples)
proc.time() - ptm

return(samples.pred)

}

