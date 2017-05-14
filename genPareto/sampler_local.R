sampler_local <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        logScale    <- control$log
        adaptive      <- control$adaptive
        adaptScaleOnly <- control$adaptScaleOnly
        adaptInterval  <- control$adaptInterval
        propCov        <- control$propCov
        scale         <- control$scale
        order <- control$order
	focal <- control$focal
	nbhd <- control$nbhd
	nbrs <- seq(focal-nbhd, focal+nbhd, by = 1)
        procName <- control$procName
        varName <- control$varName
        varFocal <- target[2]
        varNodes <- model$expandNodeNames(varName)
        precFocal <- target[1]
        Dmat <- control$Dmat        
        work_prec <- control$work_prec
        work_vinvm <- control$work_vinvm

	d <- 2
	procNodes <- model$expandNodeNames(control$procName)
        TT <- as.double(length(procNodes))
	allowed <- 1:TT
        nbrs <- intersect(nbrs, allowed)	
	lower <- if(focal-nbhd-1 >= 1) 1:(focal-nbhd-1) else NULL
	upper <- if(focal+nbhd+1 <= TT) (focal+nbhd+1):TT else NULL
	nonnbrs <- as.double(c(lower,upper))
	procTarget <- paste0(procName, '[', nbrs, ']')
	procNonTarget <- paste0(procName, '[', nonnbrs, ']')
        calcNodes  <- model$getDependencies(c(procTarget, varFocal, precFocal))
        ## numeric value generation
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        ##scaleHistory  <- c(0, 0)   ## scaleHistory
        optimalAR     <- 0.44
        gamma1        <- 0
	zeroes <- rep(0, d)
        if(is.character(propCov) && propCov == 'identity')     propCov <- diag(d)
        propCovOriginal <- propCov
        chol_propCov <- chol(propCov)
        chol_propCov_scale <- scale * chol_propCov
        empirSamp <- matrix(0, nrow=adaptInterval, ncol=d)
       	my_calcAdaptationFactor <- calcAdaptationFactor(d)
        ## checks
        if(class(propCov) != 'matrix')        stop('propCov must be a matrix\n')
        if(class(propCov[1,1]) != 'numeric')  stop('propCov matrix must be numeric\n')
        if(!all(dim(propCov) == d))           stop('propCov matrix must have dimension ', d, 'x', d, '\n')
        if(!isSymmetric(propCov))             stop('propCov matrix must be symmetric')
 
    },
    run = function() {
      currentValue <- c(model[[precFocal]], model[[varFocal]])
      jumpBack <- -sum(log(currentValue))
      
      tmp = Dmat * matrix(sqrt(1/values(model,varNodes)[(order+1):TT]), nrow = TT-order, ncol = TT)
      Q = t(tmp) %*% tmp
      diag(Q)=diag(Q)+work_prec
      V = inverse(Q)
      m = (V %*% work_vinvm)[,1] # solve(Q, work_vinvm) # might use V instead
      
      proc2 = values(model,procNonTarget)
      mu2 = m[nonnbrs]
      Sigma22 = V[nonnbrs,nonnbrs]
      Sigma11 = V[nbrs,nbrs]
      Sigma12 = V[nbrs,nonnbrs]

      # could get cholesky and use fwd/backsolve
      mu1_2 = m[nbrs]+(Sigma12 %*% solve(Sigma22, proc2-mu2))[,1]
      ch_V1_2 = chol(Sigma11 - Sigma12 %*%solve(Sigma22, t(Sigma12)))
      jumpBack <- jumpBack + dmnorm_chol(values(model,procTarget), mu1_2, ch_V1_2, log = TRUE, prec_param = FALSE)
      
      if(logScale) { propLogScale <- rmnorm_chol(1, zeroes, chol_propCov_scale, prec_param = FALSE)
                     propValue <- currentValue * exp(propLogScale)
                   } else         propValue <- rmnorm_chol(1, mean = currentValue, chol_propCov_scale, prec_param = FALSE) 
      
      model[[precFocal]] <<- propValue[1]
      model[[varFocal]] <<- propValue[2]
      jumpForward <- -sum(log(propValue))
      
      tmp = Dmat * matrix(sqrt(1/values(model,varNodes)[(order+1):TT]), nrow = TT-order, ncol = TT)
      Q = t(tmp)%*%tmp
      diag(Q)=diag(Q)+work_prec
      V = inverse(Q)
      m = (V %*% work_vinvm)[,1] # solve(Q, work_vinvm)

      mu2 = m[nonnbrs]
      Sigma22 = V[nonnbrs,nonnbrs]
      Sigma11 = V[nbrs,nbrs]
      Sigma12 = V[nbrs,nonnbrs]
      
      mu1_2 = m[nbrs]+(Sigma12 %*% solve(Sigma22, proc2-mu2))[,1]
      ch_V1_2 = chol(Sigma11 - Sigma12 %*%solve(Sigma22, t(Sigma12)))
      values(model, procTarget) <<- rmnorm_chol(1, mu1_2, ch_V1_2, prec_param = FALSE)
      
      jumpForward <- jumpForward + dmnorm_chol(values(model, procTarget), mu1_2, ch_V1_2, log = TRUE, prec_param = FALSE)

      logMHR <- calculateDiff(model, calcNodes) - jumpForward + jumpBack
      jump <- decide(logMHR)
      if(jump) nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
      else     nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
      if(adaptive)     adaptiveProcedure(jump)
    },
    methods = list(
        adaptiveProcedure = function(jump = logical()) {
            timesRan <<- timesRan + 1
            if(jump)     timesAccepted <<- timesAccepted + 1
            if(timesRan %% adaptInterval == 0) {
                acceptanceRate <- timesAccepted / timesRan
                timesAdapted <<- timesAdapted + 1
                gamma1 <<- 1/((timesAdapted + 3)^0.8)
                gamma2 <- 10 * gamma1
                adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
                scale <<- scale * adaptFactor
                timesRan <<- 0
                timesAccepted <<- 0
            }
        },
       reset = function() {
            scale <<- scaleOriginal
            timesRan      <<- 0
            timesAccepted <<- 0
            timesAdapted  <<- 0
            gamma1 <<- 0
        }
    ), where = getLoadingNamespace()
)
