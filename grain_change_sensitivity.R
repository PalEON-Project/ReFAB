### From validation 51, 95 are bimodal sites I am going to test on
### 58, 15 are going to be the unimodal sites
### Two things we're doing:
### 1. add and subtract pollen grains from each taxa and plot the max 
### lik estimate vs. the pollen counts


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

load("twothirds_v1.0.Rdata")

for(k in c(93,94,95,15,58)){
which_site <- k
dat51 <- Y[which_site,]
source(file.path('Workflow_Code','utils','getLik.R'))
outLik.real<- getLik(Z = Z.new, u = u, beta = colMeans(samples.mixed),
                     bMax = bMax, Y = dat51)
dat51.mat <- matrix(as.numeric(dat51), ncol(Y)*20, ncol(Y),byrow=T)
plus.minus.vec <- c(seq(-10,-1,1),seq(1,10,1))
index <- 1
for(i in 1:ncol(Y)){
    index <- seq(index,index+(length(plus.minus.vec)-1),1)
    dat51.mat[index,i] <- plus.minus.vec + as.numeric(dat51[i])
    index <- index[length(index)] + 1
}

dat51.mat[dat51.mat<0]<-0

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
source(file.path('Workflow_Code','utils','getLik.R'))
outLik <- getLik(Z = Z.new, u = u, beta = colMeans(samples.mixed),
                 bMax = bMax, Y = dat51.mat)

vec.name <- seq(1,440,1)
col.do1 <- colorRampPalette(colors = c('blue','grey'))
col.do2 <- colorRampPalette(colors = c('grey','red'))
colors1 <- col.do1(length(plus.minus.vec)/2)
colors2 <- col.do2(length(plus.minus.vec)/2)
all.colors <- rep(c(colors1,colors2),22)
vec.name[seq(1,440,20)]<- colnames(Y)
save.which.max <- rep(NA,440)
pdf(paste0('site_',which_site,'.liks.grain.change.pdf'))
par(mfrow=c(3,2))
for(i in seq(1,440,20)){
  plot(new.biomass,outLik[i,],typ='l',main=vec.name[i],col=all.colors[i],ylim=range(outLik))
  for(j in i:(i+19)){
    points(new.biomass,outLik[j,],typ='l',col=all.colors[j])
    save.which.max[j] <- which.max(outLik[j,])
  }
  points(new.biomass,outLik.real,typ='l',main=colnames(Y)[i],ylim=range(outLik),lwd=2)
  
  plot(dat51.mat[i:(i+19),which(seq(1,440,20)==i)], save.which.max[i:(i+19)],
       pch=19,main=vec.name[i],ylab = 'Max Lik Biomass',xlab='Grains + -',
       ylim=c(0,150),col=all.colors[i:(i+19)])
  abline(v = dat51[which(seq(1,440,20)==i)],lwd=2)
}
dev.off()
}



