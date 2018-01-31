
### after get.data

library(nimble)
library(splines)
library(maps)

ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}

code <- nimbleCode({
  
  for(r in 1:R){ 
    for(i in 1:I){
      beta[r,i] ~ dnorm(0,.04)
      beta.pine[r,i] ~ dnorm(0,.04)
    }  
  }
  
  phi.first[,] <- Z[,]%*%beta[,]
  pine.phi[,] <- Z[,]%*%beta.pine[,]
  
  for(j in 1:J){
    for(i in 1:I){
      exp.phi[j,i] <- exp(phi.first[j,i])
      exp.pine.phi[j,i] <- exp(pine.phi[j,i])
    }
  }
  
  for(j in 1:J){
    p.true[j,1] ~ dbeta(exp.phi[j,1],exp.pine.phi[j,1])
    p.rel[j,1] <- p.true[j,1]
    
    for(i in 2:(I-1)){
      p.rel[j,i]  ~ dbeta(exp.phi[j,i],exp.pine.phi[j,i]) 
      p.true[j,i] <-  p.rel[j,i] * (1 - sum(p.true[j,1:(i-1)]))
    }	
    p.true[j,I] <- 1 - sum(p.true[j,1:(I-1)])
  }  
  
  for(j in 1:J){
    Y[j,] ~ dmulti(size = n[j], prob = p.true[j,])
  }
  
  phi.first1[,] <- Z.new[,]%*%beta[,]
  pine.phi1[,] <- Z.new[,]%*%beta.pine[,]
  
  for(j in 1:145){
    for(i in 1:I){
      exp.phi1[j,i] <- exp(phi.first1[j,i])
      exp.pine.phi1[j,i] <- exp(pine.phi1[j,i])
    }
  }
  
  for(j in 1:145){
    p.true1[j,1] ~ dbeta(exp.phi1[j,1],exp.pine.phi1[j,1])
    p.rel1[j,1] <- p.true1[j,1]
    
    for(i in 2:(I-1)){
      p.rel1[j,i]  ~ dbeta(exp.phi1[j,i],exp.pine.phi1[j,i]) 
      p.true1[j,i] <-  p.rel1[j,i] * (1 - sum(p.true1[j,1:(i-1)]))
    }	
    p.true1[j,I] <- 1 - sum(p.true1[j,1:(I-1)])
  }  
  
})

load("2018-01-08calibration.data.Rdata")

source(file.path('Workflow_Code','calibration.model.R'))
calibration_model(Y = Y, biomass = biomass, code = code, Niters = 50000, DRAW = TRUE)

beta.names <- rep(rep(colnames(Y),each=5),2)
beta.i <- grep('beta',colnames(samples.mixed))

pdf('beta.posteriors.pdf')
par(mfrow=c(4,4))
for(i in seq_along(beta.i)){
  plot(samples.mixed[,beta.i[i]],typ='l',main=paste('beta',beta.i[i],'pollen',beta.names[i]))
  hist(samples.mixed[,beta.i[i]],main=paste('beta',beta.i[i]),freq=F)
}
dev.off()

ddirchmulti <- nimbleFunction(
  run = function(x = double(1), alpha = double(1), size = double(0), log = integer(0)){
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

load(file = paste0("nimble.betas_1_2_horiz_plus",Sys.Date(),".Rdata"))

validation_model(counts = counts, Z.knots = Z.knots, 
                 samples.mixed = samples.mixed, u = u,
                 Niters = 5000, biomass = biomass)



