fit_fix_sigma <- function(locn, pred_code_fix_sigma, pred_code_fix_b, 
                          order = 3, Z, u, x.meta, ten_count_use, beta1, beta2,
                          minAge = 0, maxAge = 10000, sigmaInit = 1, nIts = 10000, 
                          nItsSave = 1000, ageInterval = 100, seed = 1, bMax = 150, 
                          nbhd = 5, lik.only = NULL, control.pts, 
                          sigma, group = NULL, group.mat, override = TRUE, 
                          Nbeta=NA, ID = NA, liks.by.taxa = TRUE) {

  site_number = unique(x.meta[x.meta$site.name == locn,1])
  
  x.meta.use <- x.meta[x.meta$site.name == locn,]
  
  source(file.path('Workflow_Code','utils','test_site.R'))
  test_site(x.meta.use)
  
  ten_count_use = ten.count[which(x.meta$site.name == locn), ]
  ten_count_use[which(is.na(ten_count_use))] <- 0
  
  Y = as.matrix(ten_count_use)
  
  sample_ages <- x.meta.use$age_bacon
  age_bins <- seq(minAge, maxAge, ageInterval)
  age_index <- as.matrix(as.numeric(
    cut(sample_ages, breaks = age_bins, labels=seq(1:(length(age_bins)-1)))
  ))
  
  tmp <- data.frame(cbind(age_index, Y))
  names(tmp)[1] <- 'age_index'
  
  Y2 <- aggregate(tmp, by = list(tmp$age_index), FUN = sum)
  
  if(!is.null(group)){
    # Y2 <- Y2[-group.mat[group,],]
  }
  
  Y <- as.matrix(Y2[ , -c(1,2)])
  age_index <- Y2[,1] #group.1 is age_index because age_index gets summed in aggregate()
  
  Z_knots <- Z
  TT <- length(age_bins)-1
  I <- ncol(Y)
  K <- ncol(Z_knots)
  J <- length(age_index)
  n <- rowSums(Y)
  Zb <- matrix(NA,TT,K)
  # new_biomass <- seq(1, bMax, 1)  # needed?
  # Z_new <- matrix(0, nrow=length(new_biomass), ncol=K) # needed?
  
  #settleMean <- x.meta[x.meta[,1] == site_number, ]$SettleBiomassMean[1]
  #settleSD <- x.meta[x.meta[,1] == site_number, ]$SettleBiomassSD[1]
  
  data_pred = list(Y = Y, sigma = sigma)#, settleMean = settleMean, settleSD = settleSD)
  
  constants_pred = list(order = order, beta1 = beta1, beta2 = beta2, I = I, J = J,
                        T = TT, n = n, u = u, N0 = rep(0, (length(u)-1)), 
                        N1 = rep(0, (length(u))), N2 = rep(0, (length(u)+1)),
                        N3 = rep(0, (length(u)+2)), age_index = age_index, bMax = bMax)
  
  dimensions_pred = list(shape1 = c(TT,I), shape2 = c(TT,I), Zb = dim(Zb), Y = dim(Y),
                         beta1 = dim(beta1), beta2=dim(beta2),shape.hold1 = c(TT,I), shape.hold2 = c(TT,I))
  
  inits.pred = list(b=rep(25,TT))

  locnClean <- gsub(' ', '-', locn)
  workFile <- paste0('workInfo_', ID, '_', locnClean, '_Beta_', Nbeta, '.Rdata')

  model_pred <- nimbleModel(pred_code_fix_sigma, constants = constants_pred,
                              data = c(data_pred, list(constraint = rep(1,TT))),
                              dimensions = dimensions_pred, inits = inits.pred)
  
  # get normal approx to likelihood for all samples for the location

  source(file.path('genPareto','calc_lik_approx.R'))
  if(TRUE){
    calc_lik_approx(model = model_pred, bName = 'b', dataName = 'Y',
                    age_index, J, I, bMin = 5, bMax =  bMax-5,
                    workFile = workFile)
  }
  
  load(workFile)
  
  Cmodel_pred <- compileNimble(model_pred)
      
  set.seed(seed)
      
      # initial values
      # set precision and variance values based on priors
      lambdaInit <- rgamma(TT, 1, rate = sigmaInit)
      omegaInit <- rexp(lambdaInit^2/2)
      
      # set up initial b values based on normal approximation to likelihood,
      # and current value of omega (variances)
      M <- matrix(0, nrow = J, ncol = TT)
      M[cbind(1:J, age_index)] <- 1
      work_prec <- diag(t(M) %*% diag(1 / workVars) %*% M)
      work_vinvm <- t(M) %*% (workData / workVars)
      if(order == 1) contrasts = c(-1, 1)
      if(order == 2) contrasts = c(-1, 2, -1)
      if(order == 3) contrasts = c(-1,3,-3,1)
      if(!order %in% 1:3) stop("'order' must be 1, 2, or 3")
      D = matrix(c(rep(c(contrasts,rep(0,TT-order)), TT-order-1), contrasts), nrow =
                   TT-order, ncol = TT, byrow=TRUE)
      tmp = D / sqrt(omegaInit[(order+1):TT])
      Q = crossprod(tmp)
      diag(Q) = diag(Q) + work_prec
      mStar = solve(Q, work_vinvm)
      bInit = mStar+(t(chol(solve(Q))))%*%rnorm(TT)
      bInit[bInit > bMax] <- bMax - 5
      
      bInit[bInit<0] <- runif(length(bInit[bInit<0]),5,10)
      
      inits_pred = list(b = bInit, omega = omegaInit, lambda = lambdaInit)
      
      Cmodel_pred$setInits(inits_pred)
      
      source(file.path('genPareto','sampler_local.R'))
      mcmcConf_pred <- configureMCMC(model_pred, thin = nIts/nItsSave, nodes = 'b')
      
      for(i in 1:TT) {  # eventually should be order+1 once sort out issue with omitted lambda/omega values
        mcmcConf_pred$addSampler(target = paste0('lambda[',i,']'), type = 'RW', control = list(log=TRUE))
        mcmcConf_pred$addSampler(target = paste0('omega[',i,']'), type = 'RW', control = list(log=TRUE))
        # specialized joint sampler for omega[t], b[t], and b[t] neighbors
        if(i > 3){ 
          mcmcConf_pred$addSampler(target= c(paste0('lambda[',i,']'), paste0('omega[',i,']')), type = 'local',
                                   control = list(log = TRUE, Dmat = D, order = order,
                                                  procName = 'b', varName = 'omega',
                                                  work_prec = work_prec, work_vinvm = c(work_vinvm), focal = i, nbhd = nbhd,
                                                  propCov = matrix(c(.1,-.07,-.07,.1),2)))
        }
      }
      
      #mcmcConf_pred$addSampler('sigma', 'RW', control=list(log=TRUE))
      
      mcmcConf_pred$addMonitors(c("b", "omega", "lambda")) 
      Rmcmc_pred <- buildMCMC(mcmcConf_pred)
      
      model_pred$setInits(inits_pred)
      
      Cmcmc_pred <- compileNimble(Rmcmc_pred, project = model_pred,resetFunctions = T)
      
      Cmcmc_pred$run(nIts)
      
      samplesList = as.matrix(Cmcmc_pred$mvSamples)
      
      save(samplesList,file = paste0('samplesList_',workFile))
      # or if we want multiple runs: but need to change seed and generate different initial values
      #  samplesList <- runMCMC(mcmc = cm$Rmcmc.pred, niter = 50000, nchains = ...,
      #                      inits = ...
  

}
