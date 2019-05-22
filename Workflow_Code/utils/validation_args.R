
#### These are important constants used in all calibrations and predictions

Niters <- 10000
bMax <- 294#228
median_use <- 54
u <- c(0,median_use,bMax)

### For 10 fold cross validation only
seed <- 5
set.seed(seed)
sets10 <- matrix(sample(x = 1:150,size = 150),15,10)


### function to draw credible intervals on a plot
ciEnvelope <- function (x, ylo, yhi, ...) {
  m <- rbind(x, ylo, yhi)
  nas <- which(apply(is.na(m), 2, sum) > 0)
  if (length(nas) > 0) {
    sub.m <- list()
    for (i in seq_along(nas)) {
      if (i == 1) {
        if (nas[i] > 1) {
          sub.m[[i]] <- m[, 1:(nas[i] - 1)]
        }
      }
      else {
        if (nas[i] > (nas[i - 1] + 1)) {
          sub.m[[i]] <- m[, (nas[i - 1] + 1):(nas[i] - 
                                                1)]
        }
      }
    }
  }
  else {
    sub.m <- list(m = m)
  }
  for (i in seq_along(sub.m)) {
    x <- sub.m[[i]]["x", ]
    ylo <- sub.m[[i]]["ylo", ]
    yhi <- sub.m[[i]]["yhi", ]
    polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi), ylo[1])), 
            border = NA, ...)
  }
}
