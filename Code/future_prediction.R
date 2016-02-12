library(nimble)
library(splines)

setwd("/Users/paleolab/babySTEPPS/")

load("min.list.june24.Rdata")
load("nimble.betas.Rdata")
load(file="/Users/paleolab/babySTEPPS/Data/pol.cal.count.mnwi1.csv") 

bs_nimble <-  nimbleFunction(
    run = function(u_given = double(0), u = double(1), N0 = double(1),
                    N1 = double(1), N2 = double(1), N3 = double(1)) {
        returnType(double(1))
        
        for(i in 1:7)
            N0[i] = 0
        for(i in 1:7)
            N1[i] = 0
        for(i in 1:6)
            N2[i] = 0
        for(i in 1:5)
            N3[i] = 0

        if(u_given < u[5]){
            N0[4] = 1
        }  else {
            N0[5] = 1
        }

        for(i in 1:7){
            p = 1
            if(N0[i]==0 & N0[i+1]==0){
                N1[i] = 0
            }
            if(N0[i]!=0 & N0[i+1]==0){
                N1[i] = ((u_given - u[i]) / (u[i+p] - u[i])) * N0[i]
            }
            if(N0[i]==0 & N0[i+1]!=0){
                N1[i] = ((u[i+p+1] - u_given) / (u[i+p+1] - u[i+1]) )* N0[i+1]
            }

        }

        for(i in 1:6){
            p = 2
            if(N1[i]==0 & N1[i+1]==0){
                N2[i] = 0
            }
            if(N1[i]!=0 & N1[i+1]!=0){
                N2[i] = ((u_given - u[i]) / (u[i+p] - u[i])) * N1[i] +
                    ((u[i+p+1] - u_given) / (u[i+p+1] - u[i+1]) )* N1[i+1]
            }
            if(N1[i]!=0 & N1[i+1]==0){
                N2[i] = ((u_given - u[i]) / (u[i+p] - u[i])) * N1[i]
            }
            if(N1[i]==0 & N1[i+1]!=0){
                N2[i] = ((u[i+p+1] - u_given) / (u[i+p+1] - u[i+1]) )* N1[i+1]
            }
            
        }

        for(i in 1:5){
            p = 3
            if(N2[i]==0 & N2[i+1]==0){
                N3[i] = 0
            }
            if(N2[i]!=0 & N2[i+1]!=0){
                N3[i] = ((u_given - u[i]) / (u[i+p] - u[i])) * N2[i] +
                    ((u[i+p+1] - u_given) / (u[i+p+1] - u[i+1]) )* N2[i+1]
            }
            if(N2[i]!=0 & N2[i+1]==0){
                N3[i] = ((u_given - u[i]) / (u[i+p] - u[i])) * N2[i]
            }
            if(N2[i]==0 & N2[i+1]!=0){
                N3[i] = ((u[i+p+1] - u_given) / (u[i+p+1] - u[i+1]) )* N2[i+1]
            }
            
        }

        return(N3)
    })
    
pred_code <- nimbleCode({
  
  for(j in 1:J){
    b[j] ~ dunif(1, 400)
    #b_trunc <- b[j]
 
    Zb[j,1:5] <- bs_nimble(b[j], u[1:8], N0[1:8], N1[1:7], N2[1:6], N3[1:5])
    }

  for(i in 1:I){
  	for(j in 1:J){
  		phi.first[j,i] <- sum(Zb[j,1:5] %*% beta[1:5,i])
  	}
  }
 # phi.first[,] <- Zb[,]%*%beta[,]
  
  for(j in 1:J){
    for(i in 1:I){
      exp.phi[j,i] <- exp(phi.first[j,i])
    }
  }
  
  for(j in 1:J){
   # for(i in 1:10){
   #  phi[j,i] <- exp.phi[j,i]/sum(exp.phi[j,])
   #} 
    p[j,] ~ ddirch(exp.phi[j,]) #http://sourceforge.net/p/mcmc-jags/discussion/610037/thread/c21ef62a
    Y[j,] ~ dmulti(prob = p[j,], size = n[j])
  }
  
})

    
Z.knots = matrix(0,nrow=length(biomass),ncol=5)
u <- c(rep(.0001606565,4), quantile(biomass,c(.5)), rep(190.0118,4))

for(i in 1:length(biomass)){
    u_given <- biomass[i]
	Z.knots[i,] = bs_nimble(u_given, u, rep(0, 8), rep(0, 7), rep(0, 6), rep(0, 5))
}

x = pol.cal.count[pol.cal.count$Age>=200,]
x = x[x$Age<=3000,]

x.meta = x[,1:6]
x = x[,7:ncol(x)]

x = x[,-which(colnames(x)==c("PINUSX"))]
trees <- c("ALNUSX","JUGLANSX","ACERX","CUPRESSA","FRAXINUX","FAGUS","CYPERACE","LARIXPSEU","TSUGAX","QUERCUS","TILIA","BETULA","PICEAX","OSTRYCAR","ULMUS","ABIES","POPULUS")
other.trees <- c("TAXUS","NYSSA","CASTANEA","PLATANUS","SALIX","LIQUIDAM")
ten.count = matrix(0,nrow(x),length(trees)+3)
prairie <- c("ARTEMISIA","ASTERX","POACEAE","AMBROSIA","CHENOAMX","CORYLUS")
ten.count[,1] <- as.numeric(rowSums(x[,prairie]))
ten.count[,2] <- as.numeric(rowSums(x[,other.trees]))
ten.count[,3:(length(trees)+2)] <- as.matrix(x[,trees])
ten.count[,(length(trees)+3)] <- as.numeric(rowSums(x)) - as.numeric(rowSums(ten.count))
colnames(ten.count)<-c("prairie","other trees",trees,"other herbs")

ten.count.save = ten.count
ten.count = round(ten.count[1:500,])


J = nrow(ten.count)#nrow(Z.new)
Zb = matrix(NA,J,ncol(Z.knots))
DFS = ncol(Zb)
p = matrix(NA,J,ncol(ten.count)); phi.first = p; phi = p
#beta.est.sim = matrix(summary(csamp.sim.cal)$statistics[,1],4,ncol(Y))


beta.est = matrix(colMeans(samples[100:nrow(samples),]),ncol(Z.knots),ncol(Y)) #matrix(summary(csamp.real.cal)$statistics[,1],ncol(Z.knots),ncol(Y))# 
new.biomass = seq(1,190,1)
#Z.new = bs(new.biomass,intercept=TRUE,knots = attr(Z.knots,"knots"),Boundary.knots = attr(Z.knots,"Boundary.knots"))
Z.new = matrix(0,nrow=length(new.biomass),ncol=5)
u <- c(rep(0,4), 26, rep(190,4))

for(i in 1:length(new.biomass)){
    u_given <- new.biomass[i]
	Z.new[i,] = bs_nimble(u_given, u, rep(0, 8), rep(0, 7), rep(0, 6), rep(0, 5))
}

data.pred = list(Y = ten.count)

constants.pred = list(beta = beta.est, I = ncol(ten.count), DFS = DFS, J = J, n = rowSums(ten.count),  Z =  Z.new, u = u, N0 = rep(0, 8), N1 = rep(0, 7), N2 = rep(0, 6), N3 = rep(0, 5))

inits.pred = list(b=rep(100,J), p = matrix(1/ncol(ten.count),nrow(ten.count),ncol(ten.count)))

dimensions.pred = list(exp.phi = dim(phi), phi.first = dim(phi), p = dim(p), Zb = dim(Zb), beta = dim(beta.est), Y = dim(ten.count))

model_pred <- nimbleModel(pred_code, inits = inits.pred, constants = constants.pred, data = data.pred, dimensions = dimensions.pred)

#Cmodel.pred <- compileNimble(model_pred) #opens browser...
spec.pred <- configureMCMC(model_pred, thin = 1, print = TRUE)
spec.pred$addMonitors(c('b')) 
Rmcmc.pred <- buildMCMC(spec.pred)
cm <- compileNimble(model_pred)
Cmcmc.pred <- compileNimble(Rmcmc.pred, project = model_pred) #Error in cModel$.nodeValPointers_byGID : $ operator not defined for this S4 class

ptm <- proc.time()
Cmcmc.pred$run(5000)
samples.pred <- as.matrix(Cmcmc.pred$mvSamples)
proc.time() - ptm

save(samples.pred,file="future.samples.pred.Rdata")







