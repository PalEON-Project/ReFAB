library(nimble)
source(file.path('Workflow_Code','utils','linexp.R'))
source(file.path('genPareto','betabin.R')) # code for user-defined beta-binomial distribution
source(file.path('Workflow_Code','utils','bs_nimble.R')) # code for b-spline basis as user-defined function

pred_code <- nimbleCode({
  sigma ~ dunif(0.01,1000) #GELMAN PAPER #5
  # sigma is scale parameter of generalized Pareto, hence uniform prior, but rate parameter of gamma prior on precisions (see dgamma below)
  # .01 lower bound because higher values can get stuck and b's don't move well with very small sigma values; need to revisit to ensure we aren't keeping model way from high posterior region
  
  for(t in 1:order){
    b[t] ~ dunif(0, bMax)
    ##         b[2] ~ dunif(0, bMax) # dnorm(b[1], var = omega[2]) # if try smoothing over baseline...
    ##         b[3] ~ dunif(0, bMax) # dnorm(2*b[2]-b[1], var = omega[3])
  }
  for(t in (order+1):T) {
      b[t] ~ dnorm(3*b[t-1]-3*b[t-2]+b[t-3], var = omega[t])
      constraint[t] ~ dconstraint(b[t] > 0 & b[t] < bMax)
  }
  for(t in (order+1):T) {  # (order+1)  
       omega[t] ~ dexp(lambda[t]*lambda[t]/2)
       lambda[t] ~ dgamma(1, rate = sigma) 
  }
  # these are dummies; should be removed once Chris sorts out remaining issue with model mixing
  for(t in 1:order) {
      omega[t] ~ dunif(0,1000)
      lambda[t] ~ dunif(0,1000)
  }
 
  for(t in 1:T){
    Zb[t,1:5] <- bs_nimble(b[t], u[1:3], N0[1:2], N1[1:3], N2[1:4], N3[1:5])
  }

  # should be best to have each row of shape{1,2} be a node rather than
  # entire shape{1,2} matrices; this syntax accomplishes that
  # this is because b[t] updates have rows of shape{1,2} as dependents
  for(t in 1:T){
    shape1[t,1:I] <- exp(Zb[t,1:5] %*% beta1[1:5,1:I])
    shape2[t,1:I] <- exp(Zb[t,1:5] %*% beta2[1:5,1:I])
  }
  
  for(j in 1:J){
      Y[j, 1] ~ dbetabin(shape1[age_index[j], 1], shape2[age_index[j], 1], n[j])
      for(i in 2:(I-1)){
        Y[j, i] ~ dbetabin(shape1[age_index[j], i], shape2[age_index[j], i], n[j] - sum(Y[j,1:(i-1)]))
        
      }
  }
 
})
pred_code_fix_sigma <- nimbleCode({
  
  #settleMean ~ dnorm(mean = b[1], sd = settleSD)
  
  for(t in 1:order){
    b[t] ~ dunif(0, bMax)
    ##         b[2] ~ dunif(0, bMax) # dnorm(b[1], var = omega[2]) # if try smoothing over baseline...
    ##         b[3] ~ dunif(0, bMax) # dnorm(2*b[2]-b[1], var = omega[3])
  }
  for(t in (order+1):T) {
    b[t] ~ dnorm(3*b[t-1]-3*b[t-2]+b[t-3], var = omega[t])
    constraint[t] ~ dconstraint(b[t] > 0 & b[t] < bMax)
  }
  for(t in (order+1):T) {  # (order+1)  
    omega[t] ~ dexp(lambda[t]*lambda[t]/2)
    lambda[t] ~ dgamma(1, rate = sigma) 
  }
  # these are dummies; should be removed once Chris sorts out remaining issue with model mixing
  for(t in 1:order) {
    omega[t] ~ dunif(0,1000)
    lambda[t] ~ dunif(0,1000)
  }
  
  for(t in 1:T){
    Zb[t,1:5] <- bs_nimble(b[t], u[1:3], N0[1:2], N1[1:3], N2[1:4], N3[1:5])
  }
  
  # should be best to have each row of shape{1,2} be a node rather than
  # entire shape{1,2} matrices; this syntax accomplishes that
  # this is because b[t] updates have rows of shape{1,2} as dependents
    shape1[,] <- linexp(Zb[,] %*% beta1[,],J,I)
    shape2[,] <- linexp(Zb[,] %*% beta2[,],J,I)
  
  for(j in 1:J){
    Y[j, 1] ~ dbetabin(shape1[age_index[j], 1], shape2[age_index[j], 1], n[j])
    for(i in 2:(I-1)){
      Y[j, i] ~ dbetabin(shape1[age_index[j], i], shape2[age_index[j], i], n[j] - sum(Y[j,1:(i-1)]))
      
    }
  }
  
  
  
})

pred_code_fix_b <- nimbleCode({
#  for(t in 1:order){
#   b[t] ~ dunif(0, bMax)
    ##         b[2] ~ dunif(0, bMax) # dnorm(b[1], var = omega[2]) # if try smoothing over baseline...
    ##         b[3] ~ dunif(0, bMax) # dnorm(2*b[2]-b[1], var = omega[3])
#  }
  for(t in (order+1):T) {
#    b[t] ~ dnorm(3*b[t-1]-3*b[t-2]+b[t-3], var = omega[t])
    constraint[t] ~ dconstraint(b[t] > 0 & b[t] < bMax)
  }
  for(t in (order+1):T) {  # (order+1)  
    omega[t] ~ dexp(lambda[t]*lambda[t]/2)
    lambda[t] ~ dgamma(1, rate = sigma) 
  }
  # these are dummies; should be removed once Chris sorts out remaining issue with model mixing
  for(t in 1:order) {
    omega[t] ~ dunif(0,1000)
    lambda[t] ~ dunif(0,1000)
  }
  
  for(t in 1:T){
    Zb[t,1:5] <- bs_nimble(b[t], u[1:3], N0[1:2], N1[1:3], N2[1:4], N3[1:5])
  }
  
  # should be best to have each row of shape{1,2} be a node rather than
  # entire shape{1,2} matrices; this syntax accomplishes that
  # this is because b[t] updates have rows of shape{1,2} as dependents
  for(t in 1:T){
    shape1[t,1:I] <- exp(Zb[t,1:5] %*% beta1[1:5,1:I])
    shape2[t,1:I] <- exp(Zb[t,1:5] %*% beta2[1:5,1:I])
  }
  
  for(j in 1:J){
    Y[j, 1] ~ dbetabin(shape1[age_index[j], 1], shape2[age_index[j], 1], n[j])
    for(i in 2:(I-1)){
      Y[j, i] ~ dbetabin(shape1[age_index[j], i], shape2[age_index[j], i], n[j] - sum(Y[j,1:(i-1)]))
      
    }
  }
  
})

pred_code_no_rand_walk <- nimbleCode({
  # sigma ~ dunif(0.01,1000) #GELMAN PAPER #5
  # sigma is scale parameter of generalized Pareto, hence uniform prior, but rate parameter of gamma prior on precisions (see dgamma below)
  # .01 lower bound because higher values can get stuck and b's don't move well with very small sigma values; need to revisit to ensure we aren't keeping model way from high posterior region
  
  for(j in 1:J){
    b[j] ~ dunif(0, bMax)
    ##         b[2] ~ dunif(0, bMax) # dnorm(b[1], var = omega[2]) # if try smoothing over baseline...
    ##         b[3] ~ dunif(0, bMax) # dnorm(2*b[2]-b[1], var = omega[3])
  }
  # for(t in (order+1):T) {
  #   b[t] ~ dnorm(3*b[t-1]-3*b[t-2]+b[t-3], var = omega[t])
  #   constraint[t] ~ dconstraint(b[t] > 0 & b[t] < bMax)
  # }
  # for(t in (order+1):T) {  # (order+1)  
  #   omega[t] ~ dexp(lambda[t]*lambda[t]/2)
  #   lambda[t] ~ dgamma(1, rate = sigma) 
  # }
  # # these are dummies; should be removed once Chris sorts out remaining issue with model mixing
  # for(t in 1:order) {
  #   omega[t] ~ dunif(0,1000)
  #   lambda[t] ~ dunif(0,1000)
  # }
  
  for(j in 1:J){
    Zb[j,1:5] <- bs_nimble(b[j], u[1:3], N0[1:2], N1[1:3], N2[1:4], N3[1:5])
  }
  
  # should be best to have each row of shape{1,2} be a node rather than
  # entire shape{1,2} matrices; this syntax accomplishes that
  # this is because b[t] updates have rows of shape{1,2} as dependents
  for(j in 1:J){
    shape1[j,1:I] <- exp(Zb[j,1:5] %*% beta1[1:5,1:I])
    shape2[j,1:I] <- exp(Zb[j,1:5] %*% beta2[1:5,1:I])
  }
  
  for(j in 1:J){
    Y[j, 1] ~ dbetabin(shape1[j, 1], shape2[j, 1], n[j])
    for(i in 2:(I-1)){
      Y[j, i] ~ dbetabin(shape1[j, i], shape2[j, i], n[j] - sum(Y[j,1:(i-1)]))
      
    }
  }
  
})
