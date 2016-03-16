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
    
    
pred_code <- nimbleCode({
	
  sigma ~ dunif(0,5) #GELMAN PAPER
  
  b[1,1] ~ dunif(0,400)
  
  for(t in 2:T){
  	   b[1,t] ~ dlnorm(log(b[1,t-1]),1/sigma^2)
    }
 
  	for(t in 1:T){
      Zb[t,1:5] <- bs_nimble(b[1,t], u[1:8], N0[1:8], N1[1:7], N2[1:6], N3[1:5])
    }

  for(i in 1:I){
  	  for(t in 1:T){
  		phi.first[t,i] <- sum(Zb[t,1:5] %*% beta[1:5,i])
  		exp.phi[t,i] <- exp(phi.first[t,i])
  	  }
  }

  	for(j in 1:J){
        Y[j,1:20] ~ ddirchmulti(exp.phi[age.index[j,1], 1:20], n[j])
     }
     
  
})

    
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

choosen <- c(2550,1652,10149,1128,1598,2849,8552,1483,1410,1117,666,197,2956,2814,824,1544)

map('state', xlim=c(-98,-81), ylim=c(41,50))
    points(x.meta[,3],
       x.meta[,2], pch=19, cex=.5)
    points(x.meta[x.meta[,1]==choosen,3],
       x.meta[x.meta[,1]==choosen,2], pch=19, cex=.75,col="red")


count(x.meta[x.meta[,3]>-94&x.meta[,3]<=-93&x.meta[,2]<46&x.meta[,2]>45&x.meta$Age>9000,1])


num = count(x.meta[,1])
length(x.meta[x.meta[,1]==1128,1])

num[order(num[,2]),]

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

plot_sites = c(2550,2548,1991,7538,9953,10211,7542,1979,10132,7531,292,1701,1652,10149,269,29,824,1772,1901,9870,9876,9874,2235,8552,9826,238,2519,1771,1119,1814,1590,2933,10169,1891,10156,9873,8551,1729,9825,9867,9871,245,1548,1483,2524,1128,34,39,9768,10240,26,59,1815,318,2293,8572,23,1593,531,2231,1645,2534,1,8567,1598,1101,8573,24,2937,666,196,8574,25,10160,1679,680,10036,10034,1944,325,2956,1411,9738,9993,2864,10239,10242,10251,2239,661,518,10163,1541,8559,360)

#sigma1 = numeric(4001)
all.preds1 = numeric(5006)

#for(i in 1:16){
site_number = 1544#plot_sites[i]

ten.count.use = ten.count[which(x.meta[,1]==site_number),]
Y = as.matrix(ten.count.use)

sample.ages <- x.meta[x.meta[,1]==site_number,]$Age
age.bins <- seq(0,max(x.meta$Age),100)
age.index <- as.matrix(as.numeric(cut(sample.ages,breaks = age.bins,labels=seq(1:(length(age.bins)-1)))))

T = length(age.bins)-1
I = ncol(Y)
K = ncol(Z.knots)
J = length(age.index)

n = rowSums(Y)

Zb = matrix(NA,T,K)
phi.first = matrix(NA,T,I); exp.phi = phi.first

beta.est = matrix(colMeans(samples[100:nrow(samples),]),K,I) 
new.biomass = seq(1,190,1)
Z.new = matrix(0,nrow=length(new.biomass),ncol=K)
u <- c(rep(0,4), 26, rep(190,4))

for(i in 1:length(new.biomass)){
    u_given <- new.biomass[i]
	Z.new[i,] = bs_nimble(u_given, u, rep(0, 8), rep(0, 7), rep(0, 6), rep(0, 5))
}

data.pred = list(Y = Y)

constants.pred = list(beta = beta.est, I = I, J=J, T=T, n = n, u = u, N0 = rep(0, 8), N1 = rep(0, 7), N2 = rep(0, 6), N3 = rep(0, 5),age.index = age.index)

inits.pred = list(b=matrix(100,1,T), sigma = 4.5)

dimensions.pred = list(exp.phi = dim(exp.phi), phi.first = dim(exp.phi), Zb = dim(Zb), beta = dim(beta.est), Y = dim(Y),  b = dim(inits.pred$b))

model_pred <- nimbleModel(pred_code, inits = inits.pred, constants = constants.pred, data = data.pred, dimensions = dimensions.pred)

#Cmodel.pred <- compileNimble(model_pred) #opens browser...
spec.pred <- configureMCMC(model_pred, thin = 1, print = TRUE)
spec.pred$addMonitors(c('b')) 
Rmcmc.pred <- buildMCMC(spec.pred)
cm <- compileNimble(model_pred)
Cmcmc.pred <- compileNimble(Rmcmc.pred, project = model_pred) #Error in cModel$.nodeValPointers_byGID : $ operator not defined for this S4 class
#dyn.unload(myproject$cppProjects[[1]]$getSOName()) #cant figure it out. 

ptm <- proc.time()
Cmcmc.pred$run(30000)
samples.pred <- as.matrix(Cmcmc.pred$mvSamples)
proc.time() - ptm

dim(samples.pred)
par(mfrow=c(3,2))
for(i in 1:5){
	plot(samples.pred[,i],typ='l')
	hist(samples.pred[,i],col='gray')
}
plot(samples.pred[10000:30000,ncol(samples.pred)],typ="l")

ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}
biomassCI = apply(samples.pred[10000:30000,],2,quantile,c(0.025,0.5,0.975))

#jpeg(paste0(site_number,".jpeg"),quality=50)

par(mfrow=c(1,2))
plot(seq(1,10000,100),colMeans(samples.pred[5000:30000,1:(ncol(samples.pred)-1)]),cex=.01,ylim=range(biomassCI),ylab="Biomass",xlab="Years BP",main=site_number)
ciEnvelope(seq(1,10000,100),biomassCI[1,1:100],biomassCI[3,1:100],col="lightblue")
points(seq(1,10000,100),colMeans(samples.pred[5000:30000,1:(ncol(samples.pred)-1)]),cex=.5,pch=16)

map('state', xlim=c(-98,-81), ylim=c(41,50))
    points(x.meta[x.meta[,1]==site_number,3],
       x.meta[x.meta[,1]==site_number,2], pch=19, cex=1)  
    title(site_number)

#dev.off()
all.preds = as.matrix(cbind(matrix(0,30,6),t(samples.pred[,1:ncol(samples.pred)-1])))
all.preds[,6]<-seq(50,2950,100)
all.preds[,1]<-x.meta[x.meta[,1]==site_number,]$SiteID[1]
all.preds[,2]<-x.meta[x.meta[,1]==site_number,]$LatitudeNorth[1]
all.preds[,3]<-x.meta[x.meta[,1]==site_number,]$LongitudeWest[1]

all.preds1 <- rbind(all.preds1,all.preds)

#sigma = samples.pred[1000:5000,ncol(samples.pred)]
#sigma1 = cbind(sigma1,sigma)
}

####
#### Plots time series of estimates and site on map ####
####

#all.preds2 = all.preds1
all.preds1 = all.preds1[-1,]

par(mfrow=c(4,2))
for(i in 1:nrow(all.preds1)){
	plot(all.preds1[i,2000:5000],type="l")
	hist(all.preds1[i,2000:5000])
}


plot_sites = unique(all.preds1[order(all.preds1[,3]),1])

pdf("sites_rand_walk1.pdf")
par(mfrow=c(4,2),mar=c(1,3,3,0))

for(i in 1:length(plot_sites)){
  #Calculate Mean and Quantiles for Estiamtes
  mean.vec = rev(rowMeans(all.preds1[all.preds1[,1]==plot_sites[i],1000:5000]))
  
  #if(length(mean.vec)<10){
    sd.vec = apply(all.preds1[all.preds1[,1]==plot_sites[i],1000:5000], 1, quantile, probs = c(0.025, 0.975))

    #Biomass Estimate Time Series
    plot(all.preds1[all.preds1[,1]==plot_sites[i],6],mean.vec,ylim=c(0,max(mean.vec+sd.vec)+10),pch=21,cex=.5,xlab="AGE",ylab="Biomass",col="blue")
    segments(all.preds1[all.preds1[,1]==plot_sites[i],6],rev(sd.vec[1,]),all.preds1[all.preds1[,1]==plot_sites[i],6],rev(sd.vec[2,]))
    title(paste0("Mean Range",signif(range(mean.vec)[1],digits=4),-signif(range(mean.vec)[2],digits=4)))
    #mtext(paste("mean sigma =",signif(mean(sigma1[,i+1]),digits=3)))  

    #Map of Site
    map('state', xlim=c(-98,-81), ylim=c(41,50))
    points(all.preds1[all.preds1[,1]==plot_sites[i],3],
       all.preds1[all.preds1[,1]==plot_sites[i],2], pch=19, cex=1)  
    title(plot_sites[i])

    #Pie Charts for Min/Max Biomass Ests
    # pie.vec = ten.count[x.meta[,1]==plot_sites[i],]
    # pie(pie.vec[which(mean.vec==min(mean.vec)),],cex=.9,labels=c("PRA","OT","ALN","JUG","ACE","CUP","FRA","FAGU","CYP","LAR","TSU","QUE","TIL","BET","PIC","OST","ULM","ABI","POP","OH"),col=rainbow(20),main="MIN")

    # pie(pie.vec[which(mean.vec==max(mean.vec)),],cex=.9,labels=c("PRA","OT","ALN","JUG","ACE","CUP","FRA","FAGU","CYP","LAR","TSU","QUE","TIL","BET","PIC","OST","ULM","ABI","POP","OH"),col=rainbow(20),main="MAX")
# #}
}
#dev.off()