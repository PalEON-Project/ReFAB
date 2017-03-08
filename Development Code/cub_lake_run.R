library(nimble)
library(splines)
library(maps)
library(plyr)
library(oce)
library(RCurl)

data.dir = "/Users/paleolab/babySTEPPS/Data/"
fig.dir = "/Users/paleolab/babySTEPPS/Figures/"
model.dir = "/Users/paleolab/babySTEPPS/Code/"

setwd("/Users/paleolab/babySTEPPS/")
load("add.bacon2.Rdata")
#load("nimble.betas2016-11-22.Rdata") #load("2016-05-31nimble.betas.Rdata") 
source(paste0(model.dir,"bs_nimble.R"))

load(file = 'nimble.betas_1_22016-12-02.Rdata')

beta1.est.real = matrix(colMeans(samples.mixed[100:nrow(samples.mixed),1:105]),ncol(Z.knots),ncol(Y))
beta2.est.real = matrix(colMeans(samples.mixed[100:nrow(samples.mixed),106:210]),ncol(Z.knots),ncol(Y))

#plots a confidence interval around an x-y plot (e.g. a timeseries)
ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}
pred_code <- nimbleCode({
  
  sigma ~ dunif(0,5) #GELMAN PAPER #5
  
  # Notes -- log biomass evolution
  #logb[1,1] ~ dunif(log(1),log(145)) #this makes it hard to jump modes because it's on a different scale
  # for(t in 2:T){
  #   logb[1,t] ~ T(dnorm(logb[1,t-1]),1/sigma^2),-Inf,log(145))
  # }
  
  #First Order -- biomass evolution
  b[1,1] ~ dunif(0,145)
  for(t in 2:T){
    b[1,t] ~ T(dnorm(b[1,t-1],1/sigma^2),0,145)
  }
  
  #Second Order - Ann Version
  # b[1,1] ~ dunif(0,145)
  # b[1,2] ~ dunif(0,145)
  # for(t in 3:T){
  #   b[1,t] ~ T(dlnorm(log(2*b[1,t-1] - b[1,t-2]),1/sigma^2),0,145)
  # }
  
  #Second Order - Chris Version
  # logb[1,1] ~ dunif(log(2),log(145))
  # logb[1,2] ~ dunif(log(2),log(145))
  # for(t in 3:T){
  # logb[1,t] ~ T(dnorm(2*logb[1,t-1] - logb[1,t-2], 1/sigma^2), log(2), log(145))
  # }
  # 
  # for(t in 1:T){
  #   b[1,t] <- exp(logb[1,t])
  # }
  
  for(t in 1:T){
    Zb[t,1:5] <- bs_nimble(b[1,t], u[1:3], N0[1:2], N1[1:3], N2[1:4], N3[1:5])
  }
  
  for(i in 1:I){
    for(t in 1:T){
      phi.first[t,i] <- sum(Zb[t,1:5] %*% beta[1:5,i])
      phi.first1[t,i] <- sum(Zb[t,1:5] %*% beta1[1:5,i])
    }
  }
  
  for(t in 1:T){
    for(i in 1:I){
      exp.phi[t,i] <- exp(phi.first[t,i])
      exp.phi1[t,i] <- exp(phi.first1[t,i])
    }
  }
  
  for(t in 1:T){
    p.true[t,1] ~ dbeta(exp.phi[t,1],exp.phi1[t,1])
    p.rel[t,1] <- p.true[t,1]
    
    for(i in 2:(I-1)){
      p.rel[t,i]  ~ dbeta(exp.phi[t,i],exp.phi1[t,i]) 
      p.true[t,i] <-  p.rel[t,i] * (1 - sum(p.true[t,1:(i-1)]))
    }	
    p.true[t,21] <- 1 - sum(p.true[t,1:20])
  }    
  
  for(j in 1:J){
    Y[j,] ~ dmulti(p.true[age.index[j,1],],n[j])
  }
  
})

u<-c(rep(attr(Z,"Boundary.knots")[1],1),attr(Z,"knots"),rep(attr(Z,"Boundary.knots")[2],1))

x = new.pol1[new.pol1$age_bacon>=200,]
x = x[x$age_bacon<=10000,]

x.meta = x[,c('site.id','lat',"long","dataset.id","site.name","age_bacon")]

trees <- c("PINUSX","ALNUSX","JUGLANSX","ACERX","CUPRESSA","FRAXINUX","FAGUS","CYPERACE","LARIXPSEU","TSUGAX","QUERCUS","TILIA","BETULA","PICEAX","OSTRYCAR","ULMUS","ABIES","POPULUS")
other.trees <- c("TAXUS","NYSSA","CASTANEA","PLATANUS","SALIX","LIQUIDAM")
ten.count = matrix(0,nrow(x),length(trees)+3)
prairie <- c("ARTEMISIA","ASTERX","POACEAE","AMBROSIA","CHENOAMX","CORYLUS")
ten.count[,1] <- unlist(rowSums(x[,prairie]))
ten.count[,2] <- unlist(rowSums(x[,other.trees]))
ten.count[,3:(length(trees)+2)] <- as.matrix(x[,trees])
ten.count[,(length(trees)+3)] <- rowSums(x[,20:99]) - rowSums(ten.count)
colnames(ten.count)<-c("prairie","other trees",trees,"other herbs")

ten.count.save = ten.count
ten.count = round(ten.count.save)

counts <- Y[,rev(order(colMeans(Y)))]

ten.count <- ten.count[,colnames(counts)]

site_number = unique(x.meta[x.meta$site.name=='Cub Lake',1])
ten.count.use = ten.count[which(x.meta$site.id==site_number),]

Y = as.matrix(ten.count.use)

sample.ages <- x.meta[x.meta[,1]==site_number,]$age_bacon
age.bins <- seq(0,10000,100)
age.index <- as.matrix(as.numeric(cut(sample.ages,breaks = age.bins,labels=seq(1:(length(age.bins)-1)))))

Z.knots = Z
T = length(age.bins)-1
I = ncol(Y)
K = ncol(Z.knots)
J = length(age.index)
n = rowSums(Y)
Zb = matrix(NA,T,K)
phi.first = matrix(NA,T,I); exp.phi = phi.first
#beta.est = matrix(colMeans(samples1[100:nrow(samples1),]),K,I) 
new.biomass = seq(1,200,1)
Z.new = matrix(0,nrow=length(new.biomass),ncol=K)

u<-c(rep(attr(Z.knots,"Boundary.knots")[1],1),attr(Z.knots,"knots"),rep(attr(Z.knots,"Boundary.knots")[2],1))

data.pred = list(Y = Y)

constants.pred = list(beta = beta1.est.real, beta1 = beta2.est.real, I = I, J = J,
                      T = T, n = n, u = u, N0 = rep(0, (length(u)-1)), 
                      N1 = rep(0, (length(u))), N2 = rep(0, (length(u)+1)),
                      N3 = rep(0, (length(u)+2)), age.index = age.index)

inits.pred = list(b = matrix(10,1,T),sigma = 4.5)#logb = matrix(log(10),1,T) #b = matrix(10,1,T), 

dimensions.pred = list(exp.phi = c(T,I), exp.phi1 = c(T,I), phi.first = c(T,I),
                       phi.first1 = c(T,I), Zb = dim(Zb), Y = dim(Y),
                       p.rel = c(T,I), p.true = c(T,I)) #  b = dim(inits.pred$b),

set.seed(0)

model_pred <- nimbleModel(pred_code, inits = inits.pred, constants = constants.pred,
                          data = data.pred, dimensions = dimensions.pred)
spec.pred <- configureMCMC(model_pred, thin = 10, print = FALSE,
                           useConjugacy = FALSE,control = list(log=TRUE))#,control = list(log=TRUE)
spec.pred$removeSamplers('b')
for(i in 1:T){
  spec.pred$addSampler(paste0('b[1,',i,']'),'slice')
}
spec.pred$addMonitors(c("b")) 
Rmcmc.pred <- buildMCMC(spec.pred)
cm <- compileNimble(model_pred, Rmcmc.pred)

samplesList <- runMCMC(mcmc = cm$Rmcmc.pred, niter = 10000, nchains = 3,
                       inits = list(list(b = matrix((25),1,T), sigma = 4.5),
                                    list(b = matrix((100),1,T), sigma = 4.5),
                                    list(b = matrix((140),1,T), sigma = 4.5)))
