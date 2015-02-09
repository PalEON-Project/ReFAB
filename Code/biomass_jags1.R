model{
	
  #delta ~ dunif(0,1000)

	for(r in 1:R){ 
		for(i in 1:I){
			beta[r,i] ~ dnorm(0,.04)
		}
	}

  phi.first <- Z%*%beta
	
	for(j in 1:J){
	  for(i in 1:I){
	    exp.phi[j,i] <- exp(phi.first[j,i])
	  }
 	  row.sums[j] <- sum(exp.phi[j,])
# 	  for(i in 1:10){
# 	    phi[j,i] <- exp.phi[j,i]/(row.sums[j])
# 	  }
	}
	
	for(j in 1:J){
    p[j,] ~ ddirch(exp.phi[j,]+.025) #http://sourceforge.net/p/mcmc-jags/discussion/610037/thread/c21ef62a
		Y[j,] ~ dmulti(p[j,],n[j])
	}
	
}
