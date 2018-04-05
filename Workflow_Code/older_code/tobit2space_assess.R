library(nimble)
ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}
sampler_toggle <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    type <- control$type
    nested_sampler_name <- paste0('sampler_', type)
    control_new <- nimbleOptions('MCMCcontrolDefaultList')
    control_new[[names(control)]] <- control
    nested_sampler_list <- nimbleFunctionList(sampler_BASE)
    nested_sampler_list[[1]] <- do.call(nested_sampler_name, list(model, mvSaved, target, control_new))
    toggle <- 1
  },
  run = function() {
    if(toggle == 1)
      nested_sampler_list[[1]]$run()
  },
  methods = list(
    reset = function()
      nested_sampler_list[[1]]$reset()
  )
)

assessParams <- function(dat, Xt, mu_f_TRUE = NULL, P_f_TRUE = NULL){  
  #mu_f_TRUE and P_f_TRUE used for simulation
  
  
  #* page 6 looks more like I expected, but I’m not sure how we’re getting negative variances
  
  #* In general, the first 3 pages of pairs plots doesn’t seem to be producing anything too absurd — there’s no estimates going off to large negative values and nothing TOO far from the sample mean. That said, I do find some of the estimates to be surprising (e.g. why is the posterior for mu12 at t=2 lower than the sample mean when the ensemble shouldn’t include any zeros)
  
  imuf   <- grep("muf", colnames(dat))
  muf <- colMeans(dat[, imuf])
  mufT <- apply(Xt,2,mean)
  PfT <- cov(Xt)
  
  mufCI <- apply(dat[,imuf],2,quantile,c(0.025,0.975))
  mufTCI <- apply(Xt,2,quantile,c(0.025,0.975))
  
  hist(Xt[,1],freq=FALSE,main='Hist Xt[,1] with density line of dat[,1]',xlim=c(-1,1))
  lines(density(dat[,1]))
  
  par(mfrow=c(1,1))
  plot(mufT,muf,pch=19,ylim=range(mufCI),xlim=range(mufTCI))
  abline(a=0,b=1,lty=2)
  for(i in 1:length(muf)){
    lines(mufTCI[,i],rep(as.vector(muf)[i],2),col=i,lwd=2)
    lines(rep(as.vector(mufT)[i],2),mufCI[,i],col=i,lwd=2)
  }
  
  #muf mufT scatter plot
  par(mfrow=c(2,2))
  for(i in 1:(length(imuf)-1)){
    plot(dat[,i],dat[,i+1],xlab=paste('mu', i),ylab=paste('mu', i+1))
    #points(mu_f_TRUE[i],mu_f_TRUE[i+1],cex=3,col=2,pch=18)
    points(muf[i],muf[i+1],cex=3,col=3,pch=19)
    points(mufT[i],mufT[i+1],cex=3,col=4,pch=20)
  }
  plot.new()
  legend("topleft",legend=c("post","sampT"),col=3:4,pch = 19:20)
  #legend("topleft",legend=c("TRUE","post","sampT"),col=2:4,pch = 18:20)
  
  boxplot(Xt,xlab='State Variables',ylab='X')
  points(muf,col='red',pch=19)
  legend("topleft",legend=c("muf"),col='red',pch = 19)
  
  #cor(dat[,1:6])
  
  iPf   <- grep("pf", colnames(dat))
  Pf <- matrix(colMeans(dat[, iPf]),ncol(X),ncol(X))
  
  PfCI <- apply(dat[,iPf],2,quantile,c(0.025,0.975))
  
  diag.stopper <- diag(length(muf))
  
  par(mfrow=c(1,1))
  plot(PfT,Pf,ylim=range(PfCI),pch=19,xlab='Pf Ensemble (True)',ylab='Pf Estimated (tobit2space)')
  abline(0,1,lty=2)
  for(i in 1:length(Pf)){
    lines(rep(as.vector(PfT)[i],2),PfCI[,i],col=i,lwd=2)
    if(diag.stopper[i]==1){
      points(PfT[i],Pf[i],cex=2,pch = 7)
    }
  }
  legend('topleft','variance',pch = 7,cex=2)
  
  diag.stopper2 <- diag.stopper+1
  diag(diag.stopper2) <- 0
  
  plot(cov2cor(PfT)[which(diag.stopper2==1)],
       cov2cor(Pf)[which(diag.stopper2==1)],pch=19,
       ylab = 'Pf', xlab = 'Pft', main = 'Correlations')
  abline(a=0,b=1,lty=2)
  
  corrCI <- apply(dat[,iPf[which(diag.stopper2!=0)]],2,quantile,c(0.025,0.975))
  
  par(mfrow=c(1,1))
  plot(PfT[which(diag.stopper2!=0)],Pf[which(diag.stopper2!=0)],
       ylim=range(corrCI),pch=19,xlab='Pf Ensemble (True)',
       ylab='Pf Estimated (tobit2space)',
       main='Non-Diagonal Covariance')
  abline(a=0,b=1,lty=2)
  for(i in 1:length(Pf)){
    if(diag.stopper2[i]==1){
      lines(rep(as.vector(PfT)[i],2),PfCI[,i],col=i,lwd=2)
    }
  }
  
  par(mfrow=c(1,1))
  plot(diag(PfT)-diag(Pf),xlab='State Variable',pch=19,
       cex=2,main='Which variance changed the most?')
  
}

## Simulate Data with No Zeros

Posdef <- function (n, ev = runif(n, 0, 10)) 
{
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

nens <- 100
mu_f_TRUE <- c(-1,5,3)
D  <- sqrt(diag(0.5,length(mu_f_TRUE)))
R <- Posdef(length(mu_f_TRUE))
P_f_TRUE  <- D%*%R%*%D

X <- Xt <- mvtnorm::rmvnorm(nens,mu_f_TRUE,P_f_TRUE)
X[X<0] = 0
plot(Xt,col=3,pch=".",cex=3)
points(X)

X_all <- X
Xt_all <- Xt

intervalX <- matrix(c(0,10000),14,2,byrow = TRUE)

x.ind <- x.censored <- matrix(NA, ncol=ncol(X), nrow=nrow(X))

mu.x <- colMeans(X)
pf.x <- cov(X)

for(i in 1:ncol(pf.x)){
  if(diag(pf.x)[i]==0) diag(pf.x)[i] = min(diag(pf.x))/2
}

J <- ncol(X)
N <- nrow(X)

for(j in 1:J){
  for(n in 1:N){
    x.ind[n,j] <- as.numeric(X[n,j] > 0)
    x.censored[n,j] <- as.numeric(ifelse(X[n,j] > intervalX[j,2], 0, X[n,j]))
  }
}

## Define Model

tobit2space.model <- nimbleCode({
  for(i in 1:N){
    y.censored[i,1:J] ~ dmnorm(muf[1:J], cov = pf[1:J,1:J])
    for(j in 1:J){
      y.ind[i,j] ~ dinterval(y.censored[i,j], c[j])
    }
  }
  
  cov[1:J,1:J] <- k_0 * pf[1:J,1:J]
  muf[1:J] ~ dmnorm(mean = mu_0[1:J], cov = pf[1:J,1:J])
  
  Sigma[1:J,1:J] <- lambda_0[1:J,1:J]/nu_0
  pf[1:J,1:J] ~ dinvwish(S = Sigma[1:J,1:J], df = J)
  
  
})

## Compile Model

constants.tobit2space = list(N = N,
                             J = J)

data.tobit2space = list(y.ind = x.ind,
                        y.censored = x.censored,
                        mu_0 = rep(0,ncol(X)),
                        lambda_0 = diag(10,ncol(X)),
                        k_0 = 3, #or 3 or some measure of prior obs
                        nu_0 = 3,
                        c = rep(0,J))#or 3 or some measure of prior obs

inits.tobit2space = list(pf = diag(1,ncol(X)), muf = colMeans(X)) #
dimensions.tobit2space <- list(S = c(ncol(X),ncol(X)), cov = c(ncol(X),ncol(X)))
#set.seed(0)
#ptm <- proc.time()
tobit2space_pred <- nimbleModel(tobit2space.model, data = data.tobit2space,
                                constants = constants.tobit2space,
                                inits = inits.tobit2space,
                                dimensions = dimensions.tobit2space)
## Adding X.mod,q,r as data for building model.
conf_tobit2space <- configureMCMC(tobit2space_pred, print=TRUE, thin = 10, useConjugacy = TRUE)
conf_tobit2space$addMonitors(c("pf", "muf","y.censored")) 
## [1] conjugate_dmnorm_dmnorm sampler: X[1:5]
## important!
## this is needed for correct indexing later
samplerNumberOffset_tobit2space <- length(conf_tobit2space$getSamplers())

#for(n in 1:nens){
for(j in 1:J){
  for(n in 1:N){
    node <- paste0('y.censored[',n,',',j,']')
    conf_tobit2space$addSampler(node, 'toggle', control=list(type='RW'))
    ## could instead use slice samplers, or any combination thereof, e.g.:
    ##conf$addSampler(node, 'toggle', control=list(type='slice'))
  }
}

#}

#conf_tobit2space$addSampler('pf',type = 'conjugate')

#conf_tobit2space$printSamplers() #use to check for congugacy

## can monitor y.censored, if you wish, to verify correct behaviour
#conf_tobit2space$addMonitors('y.censored')

Rmcmc_tobit2space <- buildMCMC(conf_tobit2space)

Cmodel_tobit2space <- compileNimble(tobit2space_pred)
Cmcmc_tobit2space <- compileNimble(Rmcmc_tobit2space, project = tobit2space_pred)

for(i in 1:length(X)) {
  ## ironically, here we have to "toggle" the value of y.ind[i]
  ## this specifies that when y.ind[i] = 1,
  ## indicator variable is set to 0, which specifies *not* to sample
  valueInCompiledNimbleFunction(Cmcmc_tobit2space$samplerFunctions[[samplerNumberOffset_tobit2space+i]], 'toggle', 1-x.ind[i])
}

set.seed(0)
dat <- runMCMC(Cmcmc_tobit2space, niter = 10000, progressBar=TRUE)
pdf('tobit2space_test.pdf')
assessParams(dat = dat, Xt = X)
dev.off()
