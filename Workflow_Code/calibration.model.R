calibration_model <- function(Y, biomass, code, Niters, DRAW){
  
#load("Data/calibration.data.Rdata")
model.dir <- c('Workflow Code/')
fig.dir <- c("Figures/")
source("Workflow_Code/utils/bs_nimble.R")

counts <- Y

Z.knots = bs(biomass,intercept=TRUE,df=5)
u <- c(rep(attr(Z.knots,"Boundary.knots")[1],1),attr(Z.knots,"knots"),rep(attr(Z.knots,"Boundary.knots")[2],1))

new.biomass = seq(1,145,1)
Z.new = bs(new.biomass,intercept=TRUE,df = ncol(Z.knots))
# 0 to 150 grid points rows are biomass and each column is basis function
# 5 basis functions 
#plot emp props of key taxa based on sampling dates by site
#get a sense of the raw data how long it took for the transition

beta = matrix(NA,ncol(Z.knots),ncol(Y))
beta.pine = matrix(NA,ncol(Z.knots),1)
p = matrix(NA,nrow(Y),ncol(Y)) ; phi = p
J = nrow(Y)

data = list(Y = as.matrix(counts) ,  Z =  Z.knots, Z.new = Z.new)

constants = list(n = rowSums(counts), R = ncol(Z.knots), I = ncol(Y),
                 J = nrow(Y))

inits = list(beta = matrix(1, ncol(Z.knots), ncol(Y)),
             beta.pine = matrix(1, ncol(Z.knots), ncol(Y)), 
             p.rel = matrix(1/ncol(Y), nrow(Y), ncol(Y)),
             p.rel1 = matrix(1/ncol(Y), 145, ncol(Y)))

dimensions = list(exp.phi = dim(phi), exp.pine.phi = dim(phi),
                  phi.first = dim(phi), pine.phi = dim(phi), 
                  Z = dim(Z.knots), beta = dim(beta), beta.pine = dim(beta),
                  p.true = dim(p), Y = dim(counts), n = nrow(Y),
                  pine.dirch = dim(phi),p.rel = dim(p),
                  phi.first1 = c(145,ncol(Y)), pine.phi1 = c(145,ncol(Y)),
                  Z.new = dim(Z.new),p.true1 = c(145,ncol(Y)),p.rel1 = c(145,ncol(Y)))

# in BUGS code, to calculate the vector of basis matrix values for a given biomass, pass that biomass in as 'u_given', pass in the vector of u values for the knots and pass in N0,N1,N2,N3 of correct length - you can do this simply by providing N0,N1,N2,N3 as part of the 'constants' argument given to the 'nimbleModel' function

model <- nimbleModel(code, inits = inits, constants = constants, data = data, dimensions = dimensions)

# compiled version of the model
Cmodel <- compileNimble(model)

# set up MCMC
#4:27pm
spec <- configureMCMC(model, thin = 10, print = TRUE, useConjugacy = FALSE)
spec$addMonitors(c('beta','beta.pine','p.true1','p.rel1')) 

# set up monitoring of whatever
# model variables you want posterior samples for - by default, top level
# parameters are already included, so 'mu' in the above example would by
# default be monitored. 'psi' and 'theta' are just for illustration -
# obviously they are not part of my toy model above

# create MCMC algorithm for the model
Rmcmc <- buildMCMC(spec)

# compiled version of the MCMC
Cmcmc <- compileNimble(Rmcmc, project = model)

# run MCMC for 2000 iterations
set.seed(0)
Cmcmc$run(Niters)#50000
samples.mixed <- as.matrix(Cmcmc$mvSamples)

save(samples.mixed,file = paste0("nimble.betas_1_2_horiz_plus",Sys.Date(),".Rdata"))

i.beta <- grep("beta",colnames(samples.mixed))
i.beta.pine <- grep("beta.pine",colnames(samples.mixed))
i.beta1 <- i.beta[-i.beta.pine]

burnin <- round(.2 * nrow(samples.mixed))

beta1.est.real = matrix(colMeans(samples.mixed[burnin:nrow(samples.mixed),i.beta1]),ncol(Z),ncol(Y))
beta2.est.real = matrix(colMeans(samples.mixed[burnin:nrow(samples.mixed),i.beta.pine]),ncol(Z),ncol(Y))

save(beta1.est.real,beta2.est.real,file = paste0("simple.betas_1_2_horiz_plus",Sys.Date(),".Rdata"))

prop.quants <- matrix(NA,ncol(samples.mixed),3)
for(i in 1:ncol(samples.mixed)){
  prop.quants[i,]<-quantile(samples.mixed[,i],c(.025,.5,.975),na.rm=TRUE)
}
rownames(prop.quants)<-colnames(samples.mixed)

plot.help<- seq(0,3045,145)

if(DRAW == TRUE) pdf(file.path(fig.dir,paste0(Sys.Date(),'tele.betas.calib_horiz_plus.pdf')))
par(mfrow=c(2,2))
for(i in 1:ncol(Y)){
  plot(1:145,prop.quants[grep('p.true1',
                              colnames(samples.mixed))[(plot.help[i]+1):plot.help[i+1]],2],
       pch=21,bg='gray',main=colnames(counts)[i],ylim=range(prop.quants[grep('p.true1',colnames(samples.mixed))[(plot.help[i]+1):plot.help[i+1]],]),ylab='pollen prop',xlab='biomass')
  
  ciEnvelope(x = 1:145,
             ylo = prop.quants[grep('p.true1',colnames(samples.mixed))[(plot.help[i]+1):plot.help[i+1]],1],
             yhi = prop.quants[grep('p.true1',colnames(samples.mixed))[(plot.help[i]+1):plot.help[i+1]],3],
             col = 'lightblue') 
  points(1:145,prop.quants[grep('p.true1',
                                colnames(samples.mixed))[(plot.help[i]+1):plot.help[i+1]],2],
         pch=21,bg='gray')
  points(biomass,counts[,i]/total_counts,cex=.8,pch=19,col='blue')
  points(biomass[95],counts[95,i]/total_counts[95],cex=.8,pch=19,col='red')
}
if(DRAW == TRUE) dev.off()

plot.betas1 <- as.matrix(exp(Z.knots%*%beta1.est.real)/rowSums(exp(Z.knots%*%beta1.est.real)))
plot.betas2 <- as.matrix(exp(Z.knots%*%beta2.est.real)/rowSums(exp(Z.knots%*%beta2.est.real)))

if(DRAW == TRUE) pdf(paste0(fig.dir,paste0("tele.betas",Sys.Date(),".pdf")))
par(mfrow=c(3,3))
for(i in 1:ncol(counts)){
  plot(biomass,counts[,i]/total_counts,pch=19,cex=.7,col='black',ylab="Pollen Proportions",main=colnames(counts)[i],xlab="Biomass",ylim=c(0,max(c(plot.betas1[,i],plot.betas2[,i],counts[,i]/total_counts))))
  points(biomass,plot.betas1[,i],col='red',pch=19,cex=1)
  points(biomass,plot.betas2[,i],col='blue',pch=19,cex=1)
  abline(v=u)
}
if(DRAW == TRUE) dev.off()

}

