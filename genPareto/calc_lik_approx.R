optimFun <- nimbleFunction(
        setup = function(model, targetNode) {
        calcNodes <- model$getDependencies(targetNode, self = FALSE)
        dataNodes <- model$getDependencies(targetNode, dataOnly = TRUE)
    },
    run = function(value = double(0)) {
        model[[targetNode]] <<- value
        model$calculate(calcNodes)
        returnType(double(0))
        return(-model$getLogProb(dataNodes))
    } 
)

calclik <- nimbleFunction(
  setup = function(model, calcnode, parentnodes, timebin) {
    deps = model$getDependencies(parentnodes)
  },
  run = function(biomasses=double(1), nSims=double()) {
    n <- length(biomasses)
    returnType(double(1))
    out = numeric(n)
    for(i in 1:n) {
      model$b[timebin] <<- biomasses[i]
      model$calculate()# so deterministic dependencies of b are updated
      
      out[i] = 0
      for(j in 1:nSims) {
        model$simulate(parentnodes) # simulate p.rel and fill in resulting p.true
        model$calculate(deps) # fill in deterministic dependencies of parentnodes and calculate densities for dependencies (which should include calcnode)
        out[i] <- out[i] + exp(model$getLogProb(calcnode)) # add density values across iterations
      }
      out[i] = out[i]/nSims
    }
    return(out)
  })


calc_lik_approx <- function(model, bName, dataName, age_index, J, I, bMin, bMax, varMax = 25^2, workFile) {
  cmodel <- compileNimble(model)
  calcNodes <- cmodel$getDependencies(bName)

  # set up grid to figure out initial values
  bNodes <- cmodel$expandNodeNames(bName)
  nNodes <- length(bNodes)
  bGrid <- seq(bMin, bMax, by = 2)
  out <- matrix(0, length(bGrid), J)
  for(b in seq_along(bGrid)) {
    values(cmodel, bNodes) <- rep(bGrid[b], nNodes)
    cmodel$calculate(calcNodes)
    for(j in seq_len(J))
        out[b, j] <- cmodel$calculate(paste0(dataName, '[', j, ', 1:', I-1, ']'))
  }

  # optimize separately at all times for which there are data
  workVars <- workData <- rep(0, J)
  for(j in seq_len(J)) {
    roptimFun <- optimFun(model, paste0(bName, '[',age_index[j],']'))
    coptimFun <- compileNimble(roptimFun, project = model)
    output = optim(bGrid[which.max(out[ , j])], coptimFun$run, lower = bMin, upper = bMax, hessian = TRUE, method = 'L-BFGS-B', )
    workData[j] <- output$par
    workVars[j] <- 1/output$hessian[1,1]
    # can't trust hessian when on boundary so fix at large-ish value
    if(workData[j] > bMax-.01 || workData[j] < bMin + .01){
      workVars[j] <- varMax
    }
    if(j > 50) {
      nimble:::clearCompiled(model)  # trying to avoid too many loaded DLLs
      cmodel <- compileNimble(model)
    }
  }
  
  # lik_mat <- matrix(NA,length(age.index),length(bGrid))
  # log_lik_mat <- matrix(NA,length(age.index),length(bGrid))
  
  # for(i in 1:J){
  #   rcalclik <- calclik(model = model_pred,
  #                       calcnode = paste0('Y[',i,',1:21]'),
  #                       parentnodes = c(paste0('shape1[',age.index[i],',1:20]'),
  #                                       paste0('shape2[',age.index[i],',1:20]')),
  #                       timebin = age.index[i])
  #   ccalclik <- compileNimble(rcalclik, project = model_pred)
  #   
  #   lik_mat[i,] <- ccalclik$run(biomasses = bGrid, nSims = 10000)
  #   
  #   a <- log(lik_mat[i,])
  #   a[is.na(a)] <- 0
  #   log_lik_mat[i,] <- exp(a - max(a))/-sum(a)
  #   
  # }
  save(workData, workVars, out, file = workFile)
  nimble:::clearCompiled(model)
}  
  
calc_lik_approx_no_rand_walk <- function(model, bName, dataName, age_index, J, I, bMin, bMax, varMax = 25^2, workFile) {
  cmodel <- compileNimble(model)
  calcNodes <- cmodel$getDependencies(bName)
  
  # set up grid to figure out initial values
  bNodes <- cmodel$expandNodeNames(bName)
  nNodes <- length(bNodes)
  bGrid <- seq(bMin, bMax, by = 2)
  out <- matrix(0, length(bGrid), J)
  for(b in seq_along(bGrid)) {
    values(cmodel, bNodes) <- rep(bGrid[b], nNodes)
    cmodel$calculate(calcNodes)
    for(j in seq_len(J))
      out[b, j] <- cmodel$calculate(paste0(dataName, '[', j, ', 1:', I-1, ']'))
  }
  
  # lik_mat <- matrix(NA,length(age.index),length(bGrid))
  # log_lik_mat <- matrix(NA,length(age.index),length(bGrid))
  
  # for(i in 1:J){
  #   rcalclik <- calclik(model = model_pred,
  #                       calcnode = paste0('Y[',i,',1:21]'),
  #                       parentnodes = c(paste0('shape1[',age.index[i],',1:20]'),
  #                                       paste0('shape2[',age.index[i],',1:20]')),
  #                       timebin = age.index[i])
  #   ccalclik <- compileNimble(rcalclik, project = model_pred)
  #   
  #   lik_mat[i,] <- ccalclik$run(biomasses = bGrid, nSims = 10000)
  #   
  #   a <- log(lik_mat[i,])
  #   a[is.na(a)] <- 0
  #   log_lik_mat[i,] <- exp(a - max(a))/-sum(a)
  #   
  # }
  save(out, file = workFile)
  nimble:::clearCompiled(model)
}  
