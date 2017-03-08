## look specifically at wood and kellners
## add the calibration biomass to time series plots
## calculate sum of the alphas

##pick a pond with prediction switching between t and t+1 (Chip Bog 2000 to 3000 years) and calculate likelihood for weighed pollen proportion by 1*y_t/n ... .5 yt/n + .5 yt/n ... 0 yt/n 1yt+1/n

pred_code <- nimbleCode({
   for(j in 1:J){
    b[j] ~ dunif(0,156)
 
    Zb[j,1:5] <- bs_nimble(b[j], u[1:3], N0[1:2], N1[1:3], N2[1:4], N3[1:5])
    }

  for(i in 1:I){
  	for(j in 1:J){
  		phi.first[j,i] <- sum(Zb[j,1:5] %*% beta[1:5,i])
  	}
  }
  
  for(j in 1:J){
    for(i in 1:I){
      exp.phi[j,i] <- exp(phi.first[j,i])
    }
  }
  
  for(j in 1:J){
   Y[j,] ~ ddirchmulti(exp.phi[j,],n[j])
  }
  
})

site_number = unique(x.meta[x.meta$site.name=='Kellners Lake',1])

ten.count.use = ten.count[which(x.meta$site.id==site_number),]

Y = as.matrix(ten.count.use)

sample.ages <- x.meta[x.meta[,1]==site_number,]$age_bacon
age.bins <- seq(0,10000,100)
age.index <- as.matrix(as.numeric(cut(sample.ages,breaks = age.bins,labels=seq(1:(length(age.bins)-1)))))



Y.save <- Y[14:15,]
w <- c(rev(seq(.1,1,length.out=10)),0)
Y.keep <- matrix(0,length(w),20)
for(i in 1:length(w)){
	Y.keep[i,] <- round(Y.save[1,] * w[i] + Y.save[2,] * (1 - w[i]))
}

Y = Y.keep


J = nrow(Y)#nrow(Z.new)
Zb = matrix(NA,J,ncol(Z.knots))
DFS = ncol(Zb)
phi = matrix(NA,J,ncol(counts)); phi.first = phi;
#beta.est.sim = matrix(summary(csamp.sim.cal)$statistics[,1],4,ncol(Y))
#load("beta.samps.Rdata")
beta.est = matrix(colMeans(samples1[100:nrow(samples1),]),ncol(Z.knots),ncol(Y)) #matrix(summary(csamp.real.cal)$statistics[,1],ncol(Z.knots),ncol(Y))# 
new.biomass = seq(1,u[length(u)],1)
#Z.new = bs(new.biomass,intercept=TRUE,knots = attr(Z.knots,"knots"),Boundary.knots = attr(Z.knots,"Boundary.knots"))
Z.new = matrix(0,nrow=length(new.biomass),ncol=5)
u <- u #should have defined knots in calibration
u<-c(rep(attr(Z.knots,"Boundary.knots")[1],1),attr(Z.knots,"knots"),rep(attr(Z.knots,"Boundary.knots")[2],1))

for(i in 1:length(new.biomass)){
    u_given <- new.biomass[i]
	Z.new[i,] = bs_nimble(u_given, u=u, N0 = rep(0, (length(u)-1)),N1 = rep(0, (length(u))), N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))
}

data.pred = list(Y = as.matrix(Y))

constants.pred = list(beta = beta.est, I = ncol(Y), DFS = DFS, J = J, n = rowSums(Y),  Z =  Z.new, u = u, N0 = rep(0, (length(u)-1)),N1 = rep(0, (length(u))), N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))

inits.pred = list(b=rep(100,J))

dimensions.pred = list(exp.phi = dim(phi), phi.first = dim(phi), Zb = dim(Zb), beta = dim(beta.est), Y = dim(Y))
set.seed(0)

model_pred <- nimbleModel(pred_code, inits = inits.pred, constants = constants.pred, data = data.pred, dimensions = dimensions.pred)

cm <- compileNimble(model_pred)

vals <- 1:157
outLik = outPost = matrix(NA, 157, J)
for(j in 1:J){
	calcNodes <-  cm$getDependencies(paste0('b[',j,']'))
for(val in vals) {
    cm$b[j] <- val
    outPost[val,j] = calculate(cm,calcNodes)# cm$calculate(calcNodes)
    # likelihood portion
    outLik[val,j] =  calculate(cm,calcNodes[grep("Y", calcNodes)]) # cm$calculate(calcNodes[45])  #
}	
}

pdf('kellners_lik_validation_biomass_30-32.pdf')
par(mfrow=c(3,3))
   for(j in 1:J){
        plot(vals,exp(outLik[,j] - max(outLik[,j]))/-sum(outLik[,j])
        ,typ='l',ylab=NA,main=w[j])
    }
dev.off()

b_1 <- which.max(outLik[,1])
b_0 <- which.max(outLik[,11])

pdf('kellners_splines_30-32.pdf')
par(mfrow=c(3,3))
for(i in 1:ncol(counts)){
	plot(biomass,counts[,i]/total_counts,pch=19,cex=.7,col='black',ylab="Pollen Proportions",main=colnames(counts)[i],xlab="Biomass")
	points(biomass,plot.betas[,i],col='red',pch=19,cex=1)
	points(b_1,Y[1,i]/rowSums(Y)[1],pch=19,cex=1.5,col='green')
	points(b_0,Y[11,i]/rowSums(Y)[11],pch=19,cex=1.5,col='green')
	abline(v=u)
  }
dev.off()
