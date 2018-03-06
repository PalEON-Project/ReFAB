linexp <- nimbleFunction(
  run = function(x = double(1)) {
    returnType(double(1))
    negvals <- which(x < 0)
    x[negvals] <- exp(x[negvals])
    x[!negvals] <- x[!negvals] + 1
    return(x)
  })