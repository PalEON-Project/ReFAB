dbetabin <- nimbleFunction(
    run = function(x = double(0), alpha = double(0), beta = double(0), size = double(0), 
        log = integer(0, default = 0)) {
        returnType(double(0))
        logProb <- lgamma(size+1) - lgamma(x+1) - lgamma(size - x + 1) +
            lgamma(alpha + beta) - lgamma(alpha) - lgamma(beta) +
            lgamma(x + alpha) + lgamma(size - x + beta) - lgamma(size + alpha + beta)
        if(log) return(logProb)
        else return(exp(logProb))
    })

rbetabin <- nimbleFunction(
    run = function(n = integer(0), alpha = double(0), beta = double(0), size = double(0)) {
        returnType(double(0))
        if(n != 1) print("rbetabin only allows n = 1; using n = 1.")
        p <- rbeta(1, alpha, beta)
        return(rbinom(1, size = size, prob = p))
    })
