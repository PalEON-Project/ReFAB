model{
  
  for(j in 1:J){
    b[j] ~ dunif(1, 400)
    b_trunc[j] <- trunc(b[j])
    
    for(c in 1:DFS) {
      Zb[j,c] <- Z[b_trunc[j], c] + (b[j] - b_trunc[j]) * (Z[b_trunc[j]+1, c] - Z[b_trunc[j], c])
    }
  }
  
  phi.first <- Zb%*%beta
  
  for(j in 1:J){
    for(i in 1:I){
      exp.phi[j,i] <- exp(phi.first[j,i])
    }
  }
  
  for(j in 1:J){
   # for(i in 1:10){
   #  phi[j,i] <- exp.phi[j,i]/sum(exp.phi[j,])
   #} 
    p[j,] ~ ddirch(exp.phi[j,] + .025) #http://sourceforge.net/p/mcmc-jags/discussion/610037/thread/c21ef62a
    Y[j,] ~ dmulti(p[j,],n[j])
  }
  
}