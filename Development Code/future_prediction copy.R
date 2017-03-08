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
ten.count = round(ten.count.save)


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

save(samples.pred,file="future.samples.pred.1.15.16.Rdata")

hist(samples.pred[1000:5000,500])
hist(colMeans(samples.pred))
all.preds1 = cbind(x.meta,t(samples.pred))

plot_sites = c(1988,269,1593,1820,2933,8574,318,1101)
map('state', xlim=range(all.preds1[,3])+c(-2, 2), ylim=range(all.preds1[,2])+c(-1, 1))
for(i in 1:length(plot_sites)){
points(all.preds1[all.preds1[,1]==plot_sites[i],3],
       all.preds1[all.preds1[,1]==plot_sites[i],2], pch=19, cex=1)
text(unique(all.preds1[all.preds1[,1]==plot_sites[i],3]),
     unique(all.preds1[all.preds1[,1]==plot_sites[i],2]+.2),
     labels = unique(all.preds1[all.preds1[,1]==plot_sites[i],1]),cex=.6)
}

####
#### Plots time series of estimates and site on map ####
####

plot_sites = unique(all.preds1[order(all.preds1[,3]),1])

pdf('all_sites_ests.pdf')

par(mfrow=c(4,4),mar=c(1,3,3,0))

for(i in 1:10){
  #Calculate Mean and Quantiles for Estiamtes
  mean.vec = rowMeans(all.preds1[all.preds1[,1]==plot_sites[i],1000:5006])
  
  if(length(mean.vec)>10){
    sd.vec = apply(all.preds1[all.preds1[,1]==plot_sites[i],1000:5006], 1, quantile, probs = c(0.05, 0.9))

    #Biomass Estimate Time Series
    plot(all.preds1[all.preds1[,1]==plot_sites[i],6],mean.vec,ylim=c(0,max(mean.vec)+10),pch=21,cex=.5,xlab="AGE",ylab="Biomass",col="blue")
    segments(all.preds1[all.preds1[,1]==plot_sites[i],6],sd.vec[1,],all.preds1[all.preds1[,1]==plot_sites[i],6],sd.vec[2,])
    title(paste0("Mean Range",signif(range(mean.vec),digits=4)))

    #Map of Site
    map('state', xlim=range(all.preds1[,3])+c(-2, 2), ylim=range(all.preds1[,2])+c(-1, 1))
    points(all.preds1[all.preds1[,1]==plot_sites[i],3],
       all.preds1[all.preds1[,1]==plot_sites[i],2], pch=19, cex=1)  
    title(plot_sites[i])

    #Pie Charts for Min/Max Biomass Ests
    pie.vec = ten.count[x.meta[,1]==plot_sites[i],]
    pie(pie.vec[which(mean.vec==min(mean.vec)),],cex=.9,labels=c("PRA","OT","ALN","JUG","ACE","CUP","FRA","FAGU","CYP","LAR","TSU","QUE","TIL","BET","PIC","OST","ULM","ABI","POP","OH"),col=rainbow(20),main="MIN")

    pie(pie.vec[which(mean.vec==max(mean.vec)),],cex=.9,labels=c("PRA","OT","ALN","JUG","ACE","CUP","FRA","FAGU","CYP","LAR","TSU","QUE","TIL","BET","PIC","OST","ULM","ABI","POP","OH"),col=rainbow(20),main="MAX")
}
}

dev.off()



