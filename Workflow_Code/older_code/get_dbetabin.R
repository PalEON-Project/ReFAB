

get_dbetabin <- function(biomass, u, N0, N1, N2, N3, beta1, beta2, n, Y) {
  I <- length(Y)
  
  Zb <- numeric(5)
  
  shape.hold1 <- shape.hold2 <- shape1 <- shape2 <- numeric(I)
  
  save_beta_dens <- numeric(I)
  
  Zb[1:5] <- bs_nimble(biomass, u[1:3], N0[1:2], N1[1:3], N2[1:4], N3[1:5])
  
  shape.hold1[1:I] <- Zb[1:5] %*% beta1[1:5, 1:I]
  shape.hold2[1:I] <- Zb[1:5] %*% beta2[1:5, 1:I]
  
  for (i in 1:I) {
    shape1[i] <- linexp_vec(shape.hold1[i])
    shape2[i] <- linexp_vec(shape.hold2[i])
  }
  
  
  #for (j in 1:J) {
    save_beta_dens[1] <- dbetabin(Y[ 1],
                                     shape1[1], shape2[1],
                                     n, log=T)
    for (i in 2:(I - 1)) {
      save_beta_dens[ i] <-dbetabin(Y[ i],
                                      shape1[i], shape2[i],
                                      n - sum(Y[ 1:(i - 1)]),log=T)
      
    }
  #}

    
  return(save_beta_dens)
}
