source(file.path('genPareto','betabin.R')) # code for user-defined beta-binomial distribution
source(file.path('Workflow_Code','utils','bs_nimble.R')) # code for b-spline basis as user-defined function

order1 <- TRUE # for first order model; set to FALSE for 2nd order

pred_code <- nimbleCode({
  
  sigma ~ dunif(0,50) #GELMAN PAPER #5
  
  if(order1) {
    b[1] ~ dunif(0,145)
    for(t in 2:T){
      b[t] ~ T(dnorm(b[t-1],1/sigma^2),0,145)
    }
  } else  {
    b[1] ~ dunif(0,145)
    b[2] ~ dunif(0,145)
    for(t in 3:T){
      b[t] ~ dnorm(2*b[t-1] - b[t-2],1/sigma^2)
    }
  }
  
  for(t in 1:T){
    Zb[t,1:5] <- bs_nimble(b[t], u[1:3], N0[1:2], N1[1:3], N2[1:4], N3[1:5])
  }
  
  for(i in 1:I){
    for(t in 1:T){
      phi.first[t,i] <- sum(Zb[t,1:5] %*% beta1[1:5,i])
      phi.first1[t,i] <- sum(Zb[t,1:5] %*% beta2[1:5,i])
    }
  }
  
  for(t in 1:T){
    for(i in 1:I){
      shape1[t,i] <- exp(phi.first[t,i])
      shape2[t,i] <- exp(phi.first1[t,i])
    }
  }
  
  for(j in 1:J){
    Y[j, 1] ~ dbetabin(shape1[age_index[j], 1], shape2[age_index[j], 1], n[j])
    for(i in 2:(I-1)){
      Y[j, i] ~ dbetabin(shape1[age_index[j], i], shape2[age_index[j], i], n[j] - sum(Y[j,1:(i-1)]))
    }
  }
  
})
