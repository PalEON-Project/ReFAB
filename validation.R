validation_model <- function(Y, Z.knots, samples.mixed, u, Niters,
                             bMax = 150, group_rm = NA){
  
library(nimble)
source(file.path('genPareto','betabin.R')) # code for user-defined beta-binomial distribution
source("Workflow_Code/utils/bs_nimble.R")
  
pred_code <- nimbleCode({
  for(j in 1:J){
    b[j] ~ dunif(0,bMax)
    
    Zb[j,1:knots] <- bs_nimble(b[j], u[1:(knots-2)], N0[1:(knots-3)], N1[1:(knots-2)], N2[1:(knots-1)], N3[1:knots])
  }
  
  shape1.hold[,] <- (Zb[,] %*% beta1[,])
  shape2.hold[,] <- (Zb[,] %*% beta2[,])
  for(j in 1:J){
    shape1[j,] <- linexp(shape1.hold[j,])
    shape2[j,] <- linexp(shape2.hold[j,])
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
                      N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)),
                      knots = ncol(Z.knots))

inits.pred = list(b=rep(25,J))

dimensions.pred = list(shape1 = dim(phi), shape2 = dim(phi),
                       shape1.hold = c(nrow(Y),ncol(Y)),
                       shape2.hold = c(nrow(Y),ncol(Y)),
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

# bInit <- numeric(J)
# for(j in 1:j){
#   bInit[j] <- mean(apply(outLik[,j],2,which.max))
# }

bInit <- apply(outLik,1,which.max)
bInit[is.na(bInit)] <- 25

inits_pred = list(b = bInit)
cm$setInits(inits_pred)

if(FALSE){
vals <- 1:bMax
outLik = outPost = array(NA, dim = c(bMax, J, (ncol(Y)-1)))

for(j in 1:J){
  calcNodes <-  cm$getDependencies(paste0('b[',j,']'))
  for(val in vals) {
    # 1: set b value
    cm$b[j] <- val
    # 2: do calculate on entire model to update all values:
    calculate(cm, calcNodes)
    # 3: for aggregated likelihood, just do calculate on 'Y' 
    #so we get only the likelihood without the prior for b[j] 
    #(though is is flat so probably doesn't matter)
    calculate(cm, paste0("Y[",j,",", 1:ncol(Y),"]"))
    # likelihood portion
    # decompose into taxa Y's
    # look at a matrix of likelihood contributions to each of the taxa
    # one for each qually time points # one for the calibration outlier
    for(s in 1:(ncol(Y)-1)){
      # 4: for likelihood by taxon, need likelihood for each val, j, i triplet:
      outLik[val,j,s] =  calculate(cm,paste0("Y[",j,",", s,"]")) # cm$calculate(calcNodes[45])  #
    }
  }	
}

}

ptm <- proc.time()
set.seed(0)
Cmcmc.pred$run(Niters)
samples.pred <- as.matrix(Cmcmc.pred$mvSamples)
proc.time() - ptm

save(samples.pred, file = paste0('samples.pred.group',group_rm,'.Rdata'))

return(samples.pred)

}

