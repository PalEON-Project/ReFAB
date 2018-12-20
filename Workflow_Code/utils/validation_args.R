
#### These are important constants used in all calibrations and predictions

Niters <- 50000
bMax <- 228
median_use <- 45
u <- c(0,median_use,bMax)

### For 10 fold cross validation only
seed <- 5
set.seed(seed)
sets10 <- matrix(sample(x = 1:150,size = 150, replace = F),15,10)
