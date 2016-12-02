pred_code <- nimbleCode({
  
 # sigma ~ dunif(0,1)#GELMAN PAPER
  
 # for(t in 1:T){
  #	log(b[t]) ~ dnorm(log(b[t-1]),1/sigma^2)
  #}
  
  for(j in 1:J){
    b[j] ~ dunif(1, 400)
    #b_trunc <- b[j]
 
    Zb[j,1:5] <- bs_nimble(b[j], u[1:3], N0[1:2], N1[1:3], N2[1:4], N3[1:5])
    }

  for(i in 1:I){
  	for(j in 1:J){
  		phi.first[j,i] <- sum(Zb[j,1:5] %*% beta[1:5,i])
  	}
  }
 # phi.first[,] <- Zb[,]%*%beta[,]
  
  for(j in 1:J){
    for(i in 1:I){
      exp.phi[j,i] <- exp(phi.first[j,i])
    }
  }
  
  for(j in 1:J){
   # for(i in 1:10){
   #  phi[j,i] <- exp.phi[j,i]/sum(exp.phi[j,])
   #} 
   #p[j,] ~ ddirch(exp.phi[j,]) #http://sourceforge.net/p/mcmc-jags/discussion/610037/thread/c21ef62a
   #Y[j,] ~ dmulti(prob = p[j,], size = n[j])
   Y[j,] ~ ddirchmulti(exp.phi[j,],n[j])
  }
  
})

pred_code <- nimbleCode({
  
  sigma ~ dunif(0,1) #GELMAN PAPER
  
  for(p in 1:P){
  	for(t in 1:T){
  	   log(b[p,t]) ~ dnorm(log(b[p,t-1]),1/sigma^2)
    }
  }
  
  
  for(p in 1:P){
  	for(t in 1:T){
    Zb[p,t,1:5] <- bs_nimble(b[p,t], u[1:8], N0[1:8], N1[1:7], N2[1:6], N3[1:5])
    }
    }

  for(i in 1:I){
  	for(p in 1:P){
  	for(t in 1:T){
  		phi.first[p,t,i] <- sum(Zb[p,t,1:5] %*% beta[1:5,i])
  	}
  }
  }
 # phi.first[,] <- Zb[,]%*%beta[,]
  
  for(p in 1:P){
  	for(t in 1:T){
    for(i in 1:I){
      exp.phi[p,t,i] <- exp(phi.first[p,t,i])
    }
  }
  }
  
  for(p in 1:P){
  	for(t in 1:T){
   # for(i in 1:10){
   #  phi[j,i] <- exp.phi[j,i]/sum(exp.phi[j,])
   #} 
   #p[j,] ~ ddirch(exp.phi[j,]) #http://sourceforge.net/p/mcmc-jags/discussion/610037/thread/c21ef62a
   #Y[j,] ~ dmulti(prob = p[j,], size = n[j])
   Y[j,] ~ ddirchmulti(exp.phi[p,t,],n[p,t])
  }
  }
  
})