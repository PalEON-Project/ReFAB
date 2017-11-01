RWfit <- function(locn, pred_code, Z,
                          u, x.meta, ten.count, beta1, beta2,
                          minAge = 0, maxAge = 11000, nIts = 10000, nItsSave = 1000,
                          ageInterval = 100, seed = 1, bMax = 150, control.pts, 
                          override = TRUE) {
  
#Here's where you pick which lake you want to run
site_number = unique(x.meta[x.meta$site.name==locn,1])
ten_count_use = ten.count[which(x.meta$site.id==site_number),]
Y = as.matrix(ten_count_use)

sample_ages <- x.meta[x.meta[,1] == site_number, ]$age_bacon
age_bins <- seq(minAge, maxAge, ageInterval)
age_index <- as.matrix(as.numeric(
  cut(sample_ages, breaks = age_bins, labels=seq(1:(length(age_bins)-1)))
))

tmp <- data.frame(cbind(age_index, Y))
names(tmp)[1] <- 'age_index'

Y2 <- aggregate(tmp, by = list(tmp$age_index), FUN = sum)

Y <- as.matrix(Y2[ , -c(1,2)])
age_index <- Y2[,1]

Z_knots <- Z
TT <- length(age_bins)-1
I <- ncol(Y)
K <- ncol(Z_knots)
J <- length(age_index)
n <- rowSums(Y)
Zb <- matrix(NA,TT,K)

data.pred = list(Y = Y)

constants.pred = list(beta1 = beta1.est.real, beta2 = beta2.est.real, I = I, J = J,
                      T = TT, n = n, u = u, N0 = rep(0, (length(u)-1)), 
                      N1 = rep(0, (length(u))), N2 = rep(0, (length(u)+1)),
                      N3 = rep(0, (length(u)+2)), age_index = age_index)

inits.pred = list(b = rep(10, TT),sigma = 4.5)#logb = matrix(log(10),1,T) #b = matrix(10,1,T), 

dimensions.pred = list(exp.phi = c(TT,I), exp.phi1 = c(TT,I), phi.first = c(TT,I),
                       phi.first1 = c(TT,I), Zb = dim(Zb), Y = dim(Y))

set.seed(0)

source('~/babySTEPPS/Workflow Code/samplers/samplers.R')
model_pred <- nimbleModel(pred_code, inits = inits.pred, constants = constants.pred,
                          data = data.pred, dimensions = dimensions.pred)
spec.pred <- configureMCMC(model_pred, thin = 10,control = list(log=TRUE))#, print = FALSE,control = list(log=TRUE)


smp <- spec.pred$getSamplers()
for(i in 1:length(smp)) {
  if(smp[[i]]$name == 'RW sampler' && smp[[i]]$target != 'sigma') {
    spec.pred$removeSamplers(smp[[i]]$target)
    spec.pred$addSampler(smp[[i]]$target, type = 'RWt_trunc', control = list(log=TRUE, range = c(0,145)))
    spec.pred$addSampler(smp[[i]]$target, type = 'jointb', control = list(log = TRUE, range = c(0,145), weights = c(.7,.2)))  # this seems to help avoid getting stuck at low-lik values early in chain and leads to higher ESS, but sampling does take longer... 
  }
}

spec.pred$addMonitors(c("b")) 
Rmcmc.pred <- buildMCMC(spec.pred)
cm <- compileNimble(model_pred, Rmcmc.pred)

# don't initialize all b's at same value as that can lead to samples for sigma being driven to be very small, at least for a while
# b1 <- rnorm(T, 25, 10)
# b2 <- rnorm(T, 75, 10)
# b3 <- rnorm(T, 125, 10)
# b1[b1 < 0] <- 2
# b3[b3 > 145] <- 144

samplesList <- runMCMC(mcmc = cm$Rmcmc.pred, niter = 10000)

locnClean <- gsub(' ', '-', locn)
workFile <- paste0('workInfo_', locnClean, 'RW', '.Rda')

save(samplesList,file = paste0('samplesList_',workFile,'.Rda'))

return(apply(samplesList[200:1000,1:100],2,quantile,c(.025,.5,.975)))

}
