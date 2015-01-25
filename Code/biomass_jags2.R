model{
	
  #delta ~ dunif(0,1000)
  tau ~ dgamma(.001,.001)
  
  for(i in 1:10){
    beta[1,i] ~ dnorm(0,.0001)
  }

	for(j in 2:30){ 
		for(i in 1:10){
			beta[j,i] ~ dnorm(beta[j-1,i],tau)
		}
	}

  phi.first <- Z%*%beta
	
	for(j in 1:141){
	  for(i in 1:10){
	    exp.phi[j,i] <- exp(phi.first[j,i])
	  }
 	  row.sums[j] <- sum(exp.phi[j,])
# 	  for(i in 1:10){
# 	    phi[j,i] <- exp.phi[j,i]/(row.sums[j])
# 	  }
	}
	
	for(j in 1:141){
    p[j,] ~ ddirch(exp.phi[j,]+.025) #http://sourceforge.net/p/mcmc-jags/discussion/610037/thread/c21ef62a
		Y[j,] ~ dmulti(p[j,],n[j])
	}
	
}
