sampler_RW_trunc <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        logScale      <- control$log
        reflective    <- control$reflective
        adaptive      <- control$adaptive
        adaptInterval <- control$adaptInterval
        scale         <- control$scale
        range <- control$range
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes  <- model$getDependencies(target)
        ## numeric value generation
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        ##scaleHistory  <- c(0, 0)   ## scaleHistory
        optimalAR     <- 0.44
        gamma1        <- 0
        ## checks
        if(length(targetAsScalar) > 1)   stop('cannot use RW sampler on more than one target; try RW_block sampler')
        if(model$isDiscrete(target))     stop('cannot use RW sampler on discrete-valued target; try slice sampler')
        if(logScale & reflective)        stop('cannot use reflective RW sampler on a log scale (i.e. with options log=TRUE and reflective=TRUE')
    },
    run = function() {
        currentValue <- model[[target]]
        propLogScale <- 0
        if(logScale) { propLogScale <- rnorm(1, mean = 0, sd = scale)
                       propValue <- currentValue * exp(propLogScale)
                   } else         propValue <- rnorm(1, mean = currentValue,  sd = scale)
        if(reflective) {
            lower <- model$getBound(target, 'lower')
            upper <- model$getBound(target, 'upper')
            while(propValue < lower | propValue > upper) {
                if(propValue < lower) propValue <- 2*lower - propValue
                if(propValue > upper) propValue <- 2*upper - propValue
            }
        }
        if(propValue < range[1] | propValue > range[2]) {
            jump <- FALSE
        } else {
            model[[target]] <<- propValue
            logMHR <- calculateDiff(model, calcNodes) + propLogScale
            jump <- decide(logMHR)
        }
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
                ##setSize(scaleHistory, timesAdapted)         ## scaleHistory
                ##scaleHistory[timesAdapted] <<- scale        ## scaleHistory
                gamma1 <<- 1/((timesAdapted + 3)^0.8)
                gamma2 <- 10 * gamma1
                adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
                scale <<- scale * adaptFactor
                timesRan <<- 0
                timesAccepted <<- 0
            }
        },
        ##getScaleHistory = function() { returnType(double(1)); return(scaleHistory) },          ## scaleHistory
        ##getScaleHistoryExpanded = function() {                                                 ## scaleHistory
        ##    scaleHistoryExpanded <- numeric(timesAdapted*adaptInterval, init=FALSE)            ## scaleHistory
        ##    for(iTA in 1:timesAdapted)                                                         ## scaleHistory
        ##        for(j in 1:adaptInterval)                                                      ## scaleHistory
        ##            scaleHistoryExpanded[(iTA-1)*adaptInterval+j] <- scaleHistory[iTA]         ## scaleHistory
        ##    returnType(double(1)); return(scaleHistoryExpanded) },                             ## scaleHistory
        reset = function() {
            scale <<- scaleOriginal
            timesRan      <<- 0
            timesAccepted <<- 0
            timesAdapted  <<- 0
            ##scaleHistory  <<- scaleHistory * 0    ## scaleHistory
            gamma1 <<- 0
        }
    ), where = getLoadingNamespace()
)


sampler_RWt_trunc <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        logScale      <- control$log
        reflective    <- control$reflective
        adaptive      <- control$adaptive
        adaptInterval <- control$adaptInterval
        scale         <- control$scale
                range <- control$range
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes  <- model$getDependencies(target)
        ## numeric value generation
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        ##scaleHistory  <- c(0, 0)   ## scaleHistory
        optimalAR     <- 0.44
        gamma1        <- 0
        ## checks
        if(length(targetAsScalar) > 1)   stop('cannot use RW sampler on more than one target; try RW_block sampler')
        if(model$isDiscrete(target))     stop('cannot use RW sampler on discrete-valued target; try slice sampler')
        if(logScale & reflective)        stop('cannot use reflective RW sampler on a log scale (i.e. with options log=TRUE and reflective=TRUE')
    },
    run = function() {
        currentValue <- model[[target]]
        propLogScale <- 0
        if(logScale) { propLogScale <- rt_nonstandard(1, df = 1, mu = 0, sigma = scale)
                       propValue <- currentValue * exp(propLogScale)
                   } else         propValue <- rt_nonstandard(1, df = 1, mu = currentValue, sigma = scale)
        if(reflective) {
            lower <- model$getBound(target, 'lower')
            upper <- model$getBound(target, 'upper')
            while(propValue < lower | propValue > upper) {
                if(propValue < lower) propValue <- 2*lower - propValue
                if(propValue > upper) propValue <- 2*upper - propValue
            }
        }
        if(propValue < range[1] | propValue > range[2]) {
            jump <- FALSE
        } else {
            model[[target]] <<- propValue
            logMHR <- calculateDiff(model, calcNodes) + propLogScale
            jump <- decide(logMHR)
        }
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
                ##setSize(scaleHistory, timesAdapted)         ## scaleHistory
                ##scaleHistory[timesAdapted] <<- scale        ## scaleHistory
                gamma1 <<- 1/((timesAdapted + 3)^0.8)
                gamma2 <- 10 * gamma1
                adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
                scale <<- scale * adaptFactor
                timesRan <<- 0
                timesAccepted <<- 0
            }
        },
        ##getScaleHistory = function() { returnType(double(1)); return(scaleHistory) },          ## scaleHistory
        ##getScaleHistoryExpanded = function() {                                                 ## scaleHistory
        ##    scaleHistoryExpanded <- numeric(timesAdapted*adaptInterval, init=FALSE)            ## scaleHistory
        ##    for(iTA in 1:timesAdapted)                                                         ## scaleHistory
        ##        for(j in 1:adaptInterval)                                                      ## scaleHistory
        ##            scaleHistoryExpanded[(iTA-1)*adaptInterval+j] <- scaleHistory[iTA]         ## scaleHistory
        ##    returnType(double(1)); return(scaleHistoryExpanded) },                             ## scaleHistory
        reset = function() {
            scale <<- scaleOriginal
            timesRan      <<- 0
            timesAccepted <<- 0
            timesAdapted  <<- 0
            ##scaleHistory  <<- scaleHistory * 0    ## scaleHistory
            gamma1 <<- 0
        }
    ), where = getLoadingNamespace()
)

sampler_jointb <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        reflective    <- control$reflective
        logScale    <- control$log
        adaptive      <- control$adaptive
        adaptInterval <- control$adaptInterval
        scale         <- control$scale
        weights <- control$weights
                range <- control$range

        tmp <- strsplit(target, split = '[\\[\\]]', perl = TRUE)[[1]]
        focal <- as.numeric(tmp[2])
        var <- tmp[1]
        allNodes <- model$expandNodeNames(var)
        times <- seq(focal-length(weights), focal+length(weights), by = 1)
        weights <- c(rev(weights), 1, weights)
        include <- times > 0 & times <= length(allNodes)
        times <- times[include]
        weights <- weights[include]
        focal <- paste0(var, '[', times[weights == 1], ']')
        target <- paste0(var, '[', times, ']')
        sumwgts <- sum(weights)
        calcNodes  <- model$getDependencies(target)
        ## numeric value generation
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        ##scaleHistory  <- c(0, 0)   ## scaleHistory
        optimalAR     <- 0.44
        gamma1        <- 0
        if(logScale & reflective)        stop('cannot use reflective RW sampler on a log scale (i.e. with options log=TRUE and reflective=TRUE')

    },
    run = function() {
        currentValue <- model[[focal]]
        propLogScale <- 0
        # log scale is simply sum of weights times increment
        if(logScale) { propLogScale <- rt_nonstandard(1, df = 1, mu = 0, sigma = scale)
                       propValue <- currentValue * exp(propLogScale)
                   } else         propValue <- rt_nonstandard(1, df = 1, mu = currentValue, sigma = scale)
        if(reflective) {
            lower <- model$getBound(focal, 'lower')
            upper <- model$getBound(focal, 'upper')
            while(propValue < lower | propValue > upper) {
                if(propValue < lower) propValue <- 2*lower - propValue
                if(propValue > upper) propValue <- 2*upper - propValue
            }
        }
        incr <- propValue - currentValue
        values(model, target) <<- values(model, target) + incr*weights

        if(max(values(model,target)) > range[2] | min(values(model,target)) < range[1]) {
            jump <- FALSE
        } else {
            logMHR <- calculateDiff(model, calcNodes) + sumwgts * propLogScale
            jump <- decide(logMHR)
        }
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
                ##setSize(scaleHistory, timesAdapted)         ## scaleHistory
                ##scaleHistory[timesAdapted] <<- scale        ## scaleHistory
                gamma1 <<- 1/((timesAdapted + 3)^0.8)
                gamma2 <- 10 * gamma1
                adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
                scale <<- scale * adaptFactor
                timesRan <<- 0
                timesAccepted <<- 0
            }
        },
        ##getScaleHistory = function() { returnType(double(1)); return(scaleHistory) },          ## scaleHistory
        ##getScaleHistoryExpanded = function() {                                                 ## scaleHistory
        ##    scaleHistoryExpanded <- numeric(timesAdapted*adaptInterval, init=FALSE)            ## scaleHistory
        ##    for(iTA in 1:timesAdapted)                                                         ## scaleHistory
        ##        for(j in 1:adaptInterval)                                                      ## scaleHistory
        ##            scaleHistoryExpanded[(iTA-1)*adaptInterval+j] <- scaleHistory[iTA]         ## scaleHistory
        ##    returnType(double(1)); return(scaleHistoryExpanded) },                             ## scaleHistory
        reset = function() {
            scale <<- scaleOriginal
            timesRan      <<- 0
            timesAccepted <<- 0
            timesAdapted  <<- 0
            ##scaleHistory  <<- scaleHistory * 0    ## scaleHistory
            gamma1 <<- 0
        }
    ), where = getLoadingNamespace()
)

