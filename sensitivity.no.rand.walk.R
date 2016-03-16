library(nimble)
library(splines)
library(maps)

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
    
# set up the "d" function for the distribution
ddirchmulti <- nimbleFunction(
  run = function(x = double(1), alpha = double(1), size = double(0), log_value = integer(0)){
  returnType(double(0))
  logProb <- lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum(lgamma(alpha + x)) - lgamma(sum(alpha) + size)
    
    if(log_value) {
      return(logProb)
    } else {
      return(exp(logProb))
    }
    
   }
)

# set up the "r" function
rdirchmulti <- nimbleFunction(
run = function(n = integer(0), alpha = double(1), size = double(0)) {
returnType(double(1))
if(n != 1) nimPrint("rdirchmulti only allows n = 1; using n = 1.")
p <- rdirch(1, alpha)
return(rmulti(1, size = size, prob = p))
})

# tell NIMBLE about the newly available distribution
registerDistributions(list(ddirchmulti = list(BUGSdist = "ddirchmulti(alpha, size)", 
  types = c('value = double(1)', 'alpha = double(1)'))))
    

Z.knots = matrix(0,nrow=length(biomass),ncol=5)
u <- c(rep(.0001606565,4), quantile(biomass,c(.5)), rep(190.0118,4))

for(i in 1:length(biomass)){
    u_given <- biomass[i]
	Z.knots[i,] = bs_nimble(u_given, u, rep(0, 8), rep(0, 7), rep(0, 6), rep(0, 5))
}


x = pol.cal.count[pol.cal.count$Age>=200,]
x = x[x$Age<=10000,]

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
ten.count = round(ten.count.save)

site_number = 1128#plot_sites[i]

ten.count.use = ten.count[which(x.meta[,1]==site_number),]
Y = as.matrix(ten.count.use)

if(prairie.sensitivity ==TRUE){
	prairie.vec = ten.count.use[3,]
	forest.vec = ten.count.use[34,]
	
	prairie.mat = matrix(prairie.vec,172,20,byrow=TRUE)
	
    colnames(prairie.mat) <- names(prairie.vec)
	
	prairie.mat[1:21,1] <- c(prairie.vec[1],seq(40,49,1),seq(51,60,1))
	prairie.mat[22:28,2] <- c(1,2,seq(4,8,1))
	prairie.mat[29:35,3] <- c(1,2,seq(4,8,1))
	prairie.mat[36:40,4] <- c(seq(1,5,1))
	prairie.mat[41:45,5] <- c(seq(1,5,1))
	prairie.mat[46:50,6] <- c(seq(1,5,1))
	prairie.mat[51:56,7] <- c(1,seq(3,7,1))
	prairie.mat[57:61,8] <- c(seq(1,5,1))
	prairie.mat[62:81,9] <- c(seq(4,13,1),seq(15,24,1))
	prairie.mat[82:86,10] <- c(seq(1,5,1))
	prairie.mat[87:91,11] <- c(seq(1,5,1))
	prairie.mat[92:111,12] <- c(seq(4,13,1),seq(15,24,1))
	prairie.mat[112:116,13] <- c(seq(1,5,1))
	prairie.mat[117:122,14] <- c(seq(1,3,1),seq(5,7,1))
	prairie.mat[123:128,15] <- c(0,seq(2,6,1))
	prairie.mat[129:134,16] <- c(0,seq(2,6,1))
	prairie.mat[135:140,17] <- c(0,seq(2,6,1))
	prairie.mat[141:145,18] <- c(seq(1,5,1))
	prairie.mat[146:152,19] <- c(0,1,seq(3,7,1))
	prairie.mat[153:172,20] <- c(seq(12,21,1),seq(23,32,1))
	
	Y = prairie.mat
	counts = Y
	
	
}

J = nrow(counts)#nrow(Z.new)
Zb = matrix(NA,J,ncol(Z.knots))
DFS = ncol(Zb)
phi = matrix(NA,J,ncol(counts)); phi.first = phi;
#beta.est.sim = matrix(summary(csamp.sim.cal)$statistics[,1],4,ncol(Y))
#load("beta.samps.Rdata")
beta.est = matrix(colMeans(samples[100:nrow(samples),]),ncol(Z.knots),ncol(Y)) #matrix(summary(csamp.real.cal)$statistics[,1],ncol(Z.knots),ncol(Y))# 
new.biomass = seq(1,190,1)
#Z.new = bs(new.biomass,intercept=TRUE,knots = attr(Z.knots,"knots"),Boundary.knots = attr(Z.knots,"Boundary.knots"))
Z.new = matrix(0,nrow=length(new.biomass),ncol=5)
u <- c(rep(0,4), 26, rep(190,4))

for(i in 1:length(new.biomass)){
    u_given <- new.biomass[i]
	Z.new[i,] = bs_nimble(u_given, u, rep(0, 8), rep(0, 7), rep(0, 6), rep(0, 5))
}

data.pred = list(Y = counts)

constants.pred = list(beta = beta.est, I = ncol(counts), DFS = DFS, J = J, n = rowSums(counts),  Z =  Z.new, u = u, N0 = rep(0, 8), N1 = rep(0, 7), N2 = rep(0, 6), N3 = rep(0, 5))

inits.pred = list(b=rep(100,J))

dimensions.pred = list(exp.phi = dim(phi), phi.first = dim(phi), Zb = dim(Zb), beta = dim(beta.est), Y = dim(counts))

model_pred <- nimbleModel(pred_code, inits = inits.pred, constants = constants.pred, data = data.pred, dimensions = dimensions.pred)

#Cmodel.pred <- compileNimble(model_pred) #opens browser...
spec.pred <- configureMCMC(model_pred, thin = 1, print = TRUE)
spec.pred$addMonitors(c('b')) 
Rmcmc.pred <- buildMCMC(spec.pred)
cm <- compileNimble(model_pred)
Cmcmc.pred <- compileNimble(Rmcmc.pred, project = model_pred) #Error in cModel$.nodeValPointers_byGID : $ operator not defined for this S4 class

ptm <- proc.time()
Cmcmc.pred$run(50000)
samples.pred <- as.matrix(Cmcmc.pred$mvSamples)
proc.time() - ptm

ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}
samples.pred = samples.pred[10000:50000,]
biomassCI = apply(samples.pred,2,quantile,c(0.025,0.5,0.975))


par(mfrow=c(4,5))
for(i in 1:20){
	plot(Y[Y[,i]!=prairie.vec[i],i],colMeans(samples.pred)[Y[,i]!=prairie.vec[i]],ylim=c(0,max(biomassCI)+1),pch=16,cex=1,xlab=paste("counts",names(prairie.vec)[i]),ylab="biomass")
	ciEnvelope(Y[Y[,i]!=prairie.vec[i],i],biomassCI[1,Y[,i]!=prairie.vec[i]],biomassCI[3,Y[,i]!=prairie.vec[i]],col="lightblue")
	points(Y[Y[,i]!=prairie.vec[i],i],colMeans(samples.pred)[Y[,i]!=prairie.vec[i]],pch=16,cex=1)
	abline(h=mean(samples.pred[,1]),col="red")
	abline(h=biomassCI[1,1],col="red",lty=2)
	abline(h=biomassCI[3,1],col="red",lty=2)
	abline(v=Y[1,i],col="red")
}


