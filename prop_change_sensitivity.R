### From validation 51, 95 are bimodal sites I am going to test on
### 58, 15 are going to be the unimodal sites
### Two things we're doing:

### 2. make total pollen count = 10000 then change the proportions of 
### each taxa up and down

library(nimble)
library(splines)
library(maps)
library(methods)

source(file.path('Workflow_Code','utils','getLik.R'))
source("Workflow_Code/utils/bs_nimble.R")

ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}

load(file = paste0("beta.est.group.inTWOTHIRDS_150.Rdata"))

burnin <- round(.2 * nrow(samples.mixed))
bMax = 150
u <- c(0,30,bMax)
new.biomass <- 1:bMax
Z.new = matrix(0,nrow=length(new.biomass),ncol=5)
for(i in 1:length(new.biomass)){
  u_given <- new.biomass[i]
  Z.new[i,] = bs_nimble(u_given, u=u, N0 = rep(0, (length(u)-1)),
                        N1 = rep(0, (length(u))), 
                        N2 = rep(0, (length(u)+1)), 
                        N3 = rep(0, (length(u)+2)))
}

load("twothirds_v1.0.Rdata")

for(k in c(93,94,95,15,58)){
which_site <- k
dat51 <- Y[which_site,]
source(file.path('Workflow_Code','utils','getLik.R'))
outLik.real<- getLik(Z = Z.new, u = u, beta = colMeans(samples.mixed),
                      bMax = bMax, Y = dat51)

prop.seq <- seq(0,1,.01)
prop.mat <- array(NA,dim=c(length(prop.seq),22,22))

for(j in 1:22){
  for(i in 1:length(prop.seq)){
    props <- dat51/rowSums(dat51)
    
    val <- props[j] - prop.seq[i] 
    props[j] <- prop.seq[i]
    
    props[-j] <- props[-j] * as.numeric((1-props[j]) / (sum(props[-j])))
                                             
    prop.mat[i,,j] <- as.numeric(props)
  }
  print(j)
}

prop.mat[prop.mat<0] <- 0 

outLik <- dat51.mat <- list()
for(i in 1:22){
  dat51.mat[[i]] <- round(prop.mat[,,i] * 10000)
  
  source(file.path('Workflow_Code','utils','getLik.R'))
  outLik[[i]] <- getLik(Z = Z.new, u = u, beta = colMeans(samples.mixed),
                   bMax = bMax, Y = dat51.mat[[i]])
}

pdf(paste0('site_',which_site,'.liks.prop.change.pdf'))
par(mfrow=c(3,2))
for(i in 1:22){
  plot(new.biomass,outLik.real,typ='l',main=colnames(Y)[i],ylim=range(outLik),lwd=2)
  for(j in 1:nrow(outLik[[i]])){
    points(new.biomass,outLik[[i]][j,],typ='l',col=rainbow(length(prop.seq),start = .5,end = .99)[j])
  }
  points(new.biomass,outLik.real,typ='l',main=colnames(Y)[i],ylim=range(outLik),lwd=2)
  
  plot(dat51.mat[[i]][,i]/rowSums(dat51.mat[[i]]), apply(outLik[[i]],1,which.max),
       pch=19,main=colnames(Y)[i],ylab = 'Max Lik Biomass',xlab='Proportions',
       ylim=c(0,150),col=rainbow(length(prop.seq),start = .5,end = .99))
  abline(v = dat51[i]/sum(dat51),lwd=2)
}

dev.off()
}


