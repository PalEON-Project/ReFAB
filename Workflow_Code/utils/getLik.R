getLik <-  function(Z,u,beta,bMax,Y) {
  source(file.path('Workflow_Code','utils','linexp.R'))

    knots = 5
    #u <- c(rep(attr(Z,"Boundary.knots")[1],1),attr(Z,"knots"),rep(attr(Z,"Boundary.knots")[2],1))
    N0 = rep(0, (length(u)-1))
    N1 = rep(0, (length(u)))
    N2 = rep(0, (length(u)+1))
    N3 = rep(0, (length(u)+2))
    
    if(FALSE) {  ## this is if you want to use your nimble spline code
        Zn=Z
        for(i in 1:nrow(Z))
            Zn[i,] =bs_nimble(i-1, u[1:(knots-2)], N0[1:(knots-3)], N1[1:(knots-2)], N2[1:(knots-1)], N3[1:knots])
        Z = Zn
    }
    
    shape1 <- linexp(Z%*%matrix(beta[1:110], nrow = 5), I = ncol(Y),
                     J = nrow(Z.new))
    shape2 <- linexp(Z%*%matrix(beta[111:220], nrow = 5), I = ncol(Y),
                     J = nrow(Z.new))
    
    # for(j in 1:nrow(shape1.hold)){
    #   shape1[j,] <- linexp(shape1.hold[j,])
    #   shape2[j,] <- linexp(shape2.hold[j,])
    # }

    J <- nrow(Y)
    liks <- matrix(0, J, bMax)
    n <- rowSums(Y)
    I <- ncol(Y)
    ## this calculates full likelihoods;
    ## if done by species you'd need 'liks' to be 
    ## a 3-d array and not add in the line with the 
    ## second call to dbetabin
    for(bb in 1:bMax) { 
      print(paste('biomass =',bb))
        for(j in 1:J) {
            liks[j,bb] <- dbetabin(Y[j,1], shape1[bb, 1], shape2[bb, 1], n[j], log = TRUE)
            for(i in 2:(I-1)){
                liks[j,bb] <- liks[j,bb] + dbetabin(Y[j,i], shape1[bb, i], shape2[bb, i], n[j] - sum(Y[j,1:(i-1)]), log = TRUE)
            }
        }
    }
    return(liks)
}
