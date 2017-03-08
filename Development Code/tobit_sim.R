library(mvtnorm)
library(rjags)
library(coda)

ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}

wish.df <- function(Om,X,i,j,col){
  n = (Om[i,j]^2 + Om[i,i]*Om[j,j])/var(X[,col])
  return(n)
}


################################   FULL MODEL

nt = 50
m = c(0.01,.9,.2)
model = matrix(0,nt,3)
Y.dat = matrix(0,nt,2)
model[1,] = c(0,0,1)
Q = diag(3) #process variance
r = diag(2)*2 #observation error

for(t in 2:nt){
  model[t,] = rmvnorm(1,m*model[t-1,],Q)
}

model[model<0]<-0 #not sure if this is right
mu.f = t(rmvnorm(1,m*model[1,],Q))
mu.f[mu.f<0]<-0
Pf   = solve(cov(model))

for(t in 1:nt){
  Y.dat[t,] = rmvnorm(1,model[t,-3],r)
}

Y.dat[Y.dat<0]<-0

#### Plot data
plot(Y.dat[,1],ylim=range(Y.dat),pch=19)
lines(model[,1],lwd=2)
points(Y.dat[,2],col="blue",pch=18)
lines(model[,2],col="blue",lwd=2)

#### Storage arrays  #### need to change dimensions
aqq = array(0,dim=c(3,3,nt+1)); Sbar.save = aqq; Pf.save = aqq; q.bar.save = aqq; Pa.save = aqq
bqq = numeric(nt+1);
Sbar.CI = array(0,dim=c(3,4,nt)); q.bar.CI = Sbar.CI
dat.save = array(0,dim=c(7001,10,nt))
CI.X1 <- matrix(0,3,nt) ; CI.X2 = CI.X1


#### Tobit Model
tobit <- "
model{ 

## Analysis
y.censored  ~ dmnorm(X[1:2],r) ##cannot be partially observed -- JAGS Manual
for(i in 1:N){
     y.ind[i] ~ dinterval(y.censored[i], interval[i,])
}

X.mod ~ dmnorm(muf,pf) ## Model Forecast

## add process error
q  ~ dwish(aq,bq)
X  ~ dmnorm(X.mod,q)
Q <- inverse(q)

}"

#### initial conditions
interval=matrix(0,2,2)
interval[1,]=c(0,10^100)
interval[2,]=c(0,10^100)

bqq[1] <- 3
aqq[,,1] <- Q*bqq[1]

tobit.inits <- function() {
  sigma <- runif(1)
  list(y.censored=as.numeric(ifelse(y.ind, NA, -abs(rnorm(1, 0, sigma)))))
}

for(t in 1:50){
Y = Y.dat[t,]
y.ind <- as.numeric(Y > interval[,1])
y.censored <- as.numeric(ifelse(Y > interval[,1], Y, 0))

### analysis of model and data
#update = list(Y=Y.dat[t,], muf=mu.f, pf=Pf, aq=aqq[,,t], bq=bqq[t], r=solve(r))
update = list(interval=interval, N=length(y.ind), y.ind=y.ind,
              y.censored=y.censored, muf=mu.f, pf=Pf,
              aq=aqq[,,t], bq=bqq[t], r=solve(r))

mod <- jags.model(file=textConnection(tobit),
                  data=update,
                  n.adapt=1000,
                  n.chains=3)

jdat <- coda.samples(mod,variable.names=c("X","X.mod","q","Q"),n.iter=10000)
summary(jdat)

dat = as.matrix(jdat)
dat = dat[3000:10000,]
#dat.save[,,t] = dat
iq = grep("q",colnames(dat))
iQ = grep("Q",colnames(dat))
iX = grep("X[",colnames(dat),fixed=TRUE)
mu.a  = colMeans(dat[,iX])
mu.a[mu.a<0]<-0 #posthoc truncation #Here or when we draw ensembles??
Pa  = cov(dat[,iX])
Pa.save[,,t] = Pa
Pa[is.na(Pa)]<- 0 

CI.X1[,t] = quantile(dat[,iX[1]],c(0.025,0.5,0.975))
CI.X2[,t] = quantile(dat[,iX[2]],c(0.025,0.5,0.975))

mQ = dat[,iQ] #Sigma, Variance
mq = dat[,iq] #Omega, Precision

Sbar = matrix(apply(mQ,2,mean),2,2) #Mean Sigma, Variance
q.bar = matrix(apply(mq,2,mean),2,2) #Mean Omega, Precision

col = matrix(1:4,2,2)
WV = matrix(0,2,2)
for(i in 1:2){
  for(j in 1:2){
    WV[i,j] <- wish.df(q.bar, X = mq, i=i, j=j, col=col[i,j])
  }
}

n = mean(WV) #n + 1
if(n < 2) n = 2
V = solve(q.bar)*n

aqq[,,t+1] = V
bqq[t+1] = n

q.bar.save[,,t] = q.bar
q.bar.CI[,,t] = apply(mq,2,quantile,c(0.025,0.5,0.975))
Sbar.save[,,t] = Sbar
Sbar.CI[,,t] = apply(mQ,2,quantile,c(0.025,0.5,0.975))

## Ensemble forward simulation
Xf = rmvnorm(1000,m*mu.a,Pa)
Xf[Xf<0]<-0
mu.f = t(colMeans(Xf))
Pf = solve(cov(Xf))
Pf.save[,,t] = Pf
}

plot(bqq,xlab="Time",ylab="Degrees of Freedom of Wishart",pch=16)

### Data assimilation time series
par(mfrow=c(1,1))
plot(Y.dat[,1],ylim=range(Y.dat),pch=19,xlab="Time",ylab="Xs")
lines(model[,1],lwd=2)
col=col2rgb("darkgrey")
col1=rgb(col[1],col[2],col[3],0.4*256,maxColorValue=256)
ciEnvelope(1:nt,CI.X1[1,],CI.X1[3,],col=col1)

points(Y.dat[,2],col="blue",pch=18)
lines(model[,2],col="blue",lwd=2)
col=col2rgb("lightblue")
col1=rgb(col[1],col[2],col[3],0.4*256,maxColorValue=256)
ciEnvelope(1:nt,CI.X2[1,],CI.X2[3,],col=col1)

## how well are we estimating process error (Q) ---> should be diag(2)
par(mfrow=c(2,2))
plot(Sbar.save[1,1,],ylim=range(Sbar.CI[,1,]))
abline(h=q[1,1])
ciEnvelope(1:nt,Sbar.CI[1,1,],Sbar.CI[3,1,],col=col1)

plot(Sbar.save[1,2,],ylim=range(Sbar.CI[,2,]))
abline(h=q[1,2])
ciEnvelope(1:nt,Sbar.CI[1,2,],Sbar.CI[3,2,],col=col1)

plot(Sbar.save[2,1,],ylim=range(Sbar.CI[,3,]))
abline(h=q[2,1])
ciEnvelope(1:nt,Sbar.CI[1,3,],Sbar.CI[3,3,],col=col1)

plot(Sbar.save[2,2,],ylim=range(Sbar.CI[,4,]))
abline(h=q[2,2])
ciEnvelope(1:nt,Sbar.CI[1,4,],Sbar.CI[3,4,],col=col1)

#### Looking for autocorrelation between process covariance and forecast covariance
par(mfrow=c(2,2))
plot(Pa.save[1,1,seq(2,50,2)],Sbar.save[1,1,seq(2,50,2)],pch=16,xlab="Pa",
     ylab="Sbar",main="Element [1,1]",ylim=c(-1,9),xlim=c(.5,1.5))
points(Pa.save[1,1,seq(1,50,2)],Sbar.save[1,1,seq(1,50,2)],col="blue",pch=16)
abline(h=1)
abline(0,1)
plot(Pa.save[1,2,],Sbar.save[1,2,],pch=16,xlab="Pa",ylab="Sbar",main="Element [1,2]")
abline(h=0)
plot(Pa.save[2,1,],Sbar.save[2,1,],pch=16,xlab="Pa",ylab="Sbar",main="Element [2,1]")
abline(h=0)
plot(Pa.save[2,2,],Sbar.save[2,2,],pch=16,xlab="Pa",ylab="Sbar",main="Element [2,2]")
abline(h=1)

######## To estimate mu.f and Pf

#### Tobit Model
tobit.first <- "
model{ 
  mu ~ dmnorm(a,b)
  P ~ dwish(aq,bq) #this is what Chris is worried about... this might not be uninformative
  for(i in 1:nens){
 y.censored[i,] ~ dmnorm(mu,P) ##cannot be partially observed -- JAGS Manual
 
      for(j in 1:n.var){
        y.ind[i,j] ~ dinterval(y.censored[i,j], interval[j,])
      }
   }
}"

update = list(interval=interval, nens=50,n.var = 2, y.ind=matrix(as.numeric(model[,1:2]>0),50,2),
              y.censored=model[,1:2],a=rep(0,2),b=diag(2),aq=diag(2),bq=2)

mod <- jags.model(file=textConnection(tobit.first),
                  data=update,
                  n.adapt=1000,
                  n.chains=3)

jdat <- coda.samples(mod,variable.names=c('mu','P'),n.iter=10000)
summary(jdat)

