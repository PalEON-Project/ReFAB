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
if(length(which(colnames(Y)=='PINUSX'))>0){
	Y = Y[,-3]
	counts = counts[,-3]
}
load("nimble.betas2016-11-22.Rdata") #load("2016-05-31nimble.betas.Rdata") 
source(paste0(model.dir,"bs_nimble.R"))

#load("min.list.june24.Rdata")
#load("nimble.betas.Rdata")
#load(file="/Users/paleolab/babySTEPPS/Data/pol.cal.count.mnwi1.csv") 

#plots a confidence interval around an x-y plot (e.g. a timeseries)
ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}

# bs_nimble <-  nimbleFunction(
    # run = function(u_given = double(0), u = double(1), N0 = double(1),
                    # N1 = double(1), N2 = double(1), N3 = double(1)) {
        # returnType(double(1))
        
        # for(i in 1:7)
            # N0[i] = 0
        # for(i in 1:7)
            # N1[i] = 0
        # for(i in 1:6)
            # N2[i] = 0
        # for(i in 1:5)
            # N3[i] = 0

        # if(u_given < u[5]){
            # N0[4] = 1
        # }  else {
            # N0[5] = 1
        # }

        # for(i in 1:7){
            # p = 1
            # if(N0[i]==0 & N0[i+1]==0){
                # N1[i] = 0
            # }
            # if(N0[i]!=0 & N0[i+1]==0){
                # N1[i] = ((u_given - u[i]) / (u[i+p] - u[i])) * N0[i]
            # }
            # if(N0[i]==0 & N0[i+1]!=0){
                # N1[i] = ((u[i+p+1] - u_given) / (u[i+p+1] - u[i+1]) )* N0[i+1]
            # }

        # }

        # for(i in 1:6){
            # p = 2
            # if(N1[i]==0 & N1[i+1]==0){
                # N2[i] = 0
            # }
            # if(N1[i]!=0 & N1[i+1]!=0){
                # N2[i] = ((u_given - u[i]) / (u[i+p] - u[i])) * N1[i] +
                    # ((u[i+p+1] - u_given) / (u[i+p+1] - u[i+1]) )* N1[i+1]
            # }
            # if(N1[i]!=0 & N1[i+1]==0){
                # N2[i] = ((u_given - u[i]) / (u[i+p] - u[i])) * N1[i]
            # }
            # if(N1[i]==0 & N1[i+1]!=0){
                # N2[i] = ((u[i+p+1] - u_given) / (u[i+p+1] - u[i+1]) )* N1[i+1]
            # }
            
        # }

        # for(i in 1:5){
            # p = 3
            # if(N2[i]==0 & N2[i+1]==0){
                # N3[i] = 0
            # }
            # if(N2[i]!=0 & N2[i+1]!=0){
                # N3[i] = ((u_given - u[i]) / (u[i+p] - u[i])) * N2[i] +
                    # ((u[i+p+1] - u_given) / (u[i+p+1] - u[i+1]) )* N2[i+1]
            # }
            # if(N2[i]!=0 & N2[i+1]==0){
                # N3[i] = ((u_given - u[i]) / (u[i+p] - u[i])) * N2[i]
            # }
            # if(N2[i]==0 & N2[i+1]!=0){
                # N3[i] = ((u[i+p+1] - u_given) / (u[i+p+1] - u[i+1]) )* N2[i+1]
            # }
            
        # }

        # return(N3)
    # })
    
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
  
  b[1,1] ~ dunif(0,157)
  
  for(t in 2:T){
  	   b[1,t] ~ T(dlnorm(log(b[1,t-1]),1/sigma^2),0,157)
    }
 
  	for(t in 1:T){
      Zb[t,1:5] <- bs_nimble(b[1,t], u[1:3], N0[1:2], N1[1:3], N2[1:4], N3[1:5])
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

ten.count <- ten.count[,-c(3)]

#ten.count = x[,7:26]

hist(x.meta$age_bacon,breaks=100)

# choosen <- c(2550,1652,10149,1128,1598,2849,8552,1483,1410,1117,666,197,2956,2814,824,1544)

# map('state', xlim=c(-98,-81), ylim=c(41,50))
    # points(x.meta[,3],
       # x.meta[,2], pch=19, cex=.5)
    # points(x.meta[x.meta[,1]==choosen,3],
       # x.meta[x.meta[,1]==choosen,2], pch=19, cex=.75,col="red")


# count(x.meta[x.meta[,3]>-94&x.meta[,3]<=-93&x.meta[,2]<46&x.meta[,2]>45&x.meta$Age>9000,1])


# num = count(x.meta[,1])
# length(x.meta[x.meta[,1]==1128,1])

# num[order(num[,2]),]

# x = x[,-which(colnames(x)==c("PINUSX"))]
# trees <- c("ALNUSX","JUGLANSX","ACERX","CUPRESSA","FRAXINUX","FAGUS","CYPERACE","LARIXPSEU","TSUGAX","QUERCUS","TILIA","BETULA","PICEAX","OSTRYCAR","ULMUS","ABIES","POPULUS")
# other.trees <- c("TAXUS","NYSSA","CASTANEA","PLATANUS","SALIX","LIQUIDAM")
# ten.count = matrix(0,nrow(x),length(trees)+3)
# prairie <- c("ARTEMISIA","ASTERX","POACEAE","AMBROSIA","CHENOAMX","CORYLUS")
# ten.count[,1] <- as.numeric(rowSums(x[,prairie]))
# ten.count[,2] <- as.numeric(rowSums(x[,other.trees]))
# ten.count[,3:(length(trees)+2)] <- as.matrix(x[,trees])
# ten.count[,(length(trees)+3)] <- as.numeric(rowSums(x)) - as.numeric(rowSums(ten.count))
# colnames(ten.count)<-c("prairie","other trees",trees,"other herbs")

ten.count.save = ten.count
ten.count = round(ten.count.save)

# plot_sites = c(2550,2548,1991,7538,9953,10211,7542,1979,10132,7531,292,1701,1652,10149,269,29,824,1772,1901,9870,9876,9874,2235,8552,9826,238,2519,1771,1119,1814,1590,2933,10169,1891,10156,9873,8551,1729,9825,9867,9871,245,1548,1483,2524,1128,34,39,9768,10240,26,59,1815,318,2293,8572,23,1593,531,2231,1645,2534,1,8567,1598,1101,8573,24,2937,666,196,8574,25,10160,1679,680,10036,10034,1944,325,2956,1411,9738,9993,2864,10239,10242,10251,2239,661,518,10163,1541,8559,360)

#sigma1 = numeric(4001)
all.preds1 = numeric(5006)

samples.pred.save = list()
biomassCI=list()
inits.keep <- list()
inits.keep.approx <- list()
outLik.save<-list()
outPost.save<-list()
#load("biomassCI.Rdata")

#do.sites<-c('WintergreenLake','LakeMendota','RadtkeLake','GassLake','KellnersLake','CubLake','GreenLake','NelsonLake')

#)2:length(do.sites) #### length(unique(x.meta[,1]))
for(s in 99:length(unique(x.meta[,1]))){
	print(paste('working on number',s,'of 183',(s/183)*100,'% complete'))
#site_number = unique(x.meta[x.meta$sitename==do.sites[s],1])
#site_number = unique(x.meta[x.meta$sitename=='GassLake',1])
#site_number = unique(x.meta[x.meta$sitename=='KellysHollow',1])
site_number = unique(x.meta[,1])[s]

ten.count.use = ten.count[which(x.meta$site.id==site_number),]
if(length(ten.count.use)>25*20){

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
beta.est = matrix(colMeans(samples1[100:nrow(samples1),]),K,I) 
new.biomass = seq(1,200,1)
Z.new = matrix(0,nrow=length(new.biomass),ncol=K)

u<-c(rep(attr(Z.knots,"Boundary.knots")[1],1),attr(Z.knots,"knots"),rep(attr(Z.knots,"Boundary.knots")[2],1))

data.pred = list(Y = Y)

constants.pred = list(beta = beta.est, I = I, J=J, T=T, n = n, u = u, N0 = rep(0, (length(u)-1)),N1 = rep(0, (length(u))), N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)), age.index = age.index)

#add random normal noise
inits.pred = list(b=matrix(100,1,T), sigma = 4.5)

dimensions.pred = list(exp.phi = dim(exp.phi), phi.first = dim(exp.phi), Zb = dim(Zb), beta = dim(beta.est), Y = dim(Y),  b = dim(inits.pred$b))

set.seed(0)

model_pred <- nimbleModel(pred_code, inits = inits.pred, constants = constants.pred, data = data.pred, dimensions = dimensions.pred)

cm <- compileNimble(model_pred)

vals <- 1:157
outLik = outPost = matrix(NA, 157, max(age.index))
for(age in age.index){
	calcNodes <-  cm$getDependencies(paste0('b[1, ',age,']'))
for(val in vals) {
    cm$b[1, age] <- val
    outPost[val,age] = calculate(cm,calcNodes)# cm$calculate(calcNodes)
    # likelihood portion
    outLik[val,age] =  calculate(cm,calcNodes[grep("Y", calcNodes)]) # cm$calculate(calcNodes[45])  #
}	
}

outLik.save[[s]] <- outLik
outPost.save[[s]] <- outPost

inits.keep[[s]] <- matrix(NA,1,T)

for(age in age.index){
	normliks <- exp(outLik[,age] - max(outLik[,age]))/-sum(outLik[,age])
    inits.keep[[s]][,age]<-which(normliks==max(normliks))
	}
inits.keep.approx[[s]] <- matrix(approx(x=1:100,y=inits.keep[[s]],xout=1:100,rule=2)$y,1,T)

inits.pred = list(b=inits.keep.approx[[s]], sigma = 4.5)

set.seed(0)

model_pred <- nimbleModel(pred_code, inits = inits.pred, constants = constants.pred, data = data.pred, dimensions = dimensions.pred)

cm <- compileNimble(model_pred)

spec.pred <- configureMCMC(model_pred, thin = 10, print = FALSE)
spec.pred$getSamplers()
spec.pred$addMonitors(c('b')) 
Rmcmc.pred <- buildMCMC(spec.pred)

Cmcmc.pred <- compileNimble(Rmcmc.pred, project = model_pred) 
Cmcmc.pred$run(10000)

samples.pred <- as.matrix(Cmcmc.pred$mvSamples)
samples.pred.save[[s]]<-samples.pred


# rug(sample.ages/100)
dyn.unload(model_pred$nimbleProject$cppProjects[[1]]$getSOName())
biomassCI[[s]] <- apply(samples.pred[250:1000,1:99],2,quantile,c(0.025,0.5,0.975))
#samples.pred.good.inits.save[[s]] <- samples.pred.good.inits

}
save(biomassCI, file="biomass.CI9.Rdata")
save(samples.pred.save, file="samples.pred.save1.Rdata")
save(inits.keep,file = 'inits.keep1.Rdata')
save(outLik.save, file = 'outLik.save1.Rdata')
save(outPost.save,file = 'outPost.save1.Rdata')
save(inits.keep.approx, file= 'inits.keep.approx1.Rdata')
}

#####
##### Plot Normalized Likelihood #####
#####

pdf(paste0('all.likes',Sys.Date(),'.pdf'))
for(s in 1:length(biomassCI)){

	if(length(biomassCI[[s]])>1){
	
	site_number <- unique(x.meta[,1])[s]
	site_name <- x.meta[x.meta[,1]==site_number,'site.name'][1]
    sample.ages <- x.meta[x.meta[,1]==site_number,]$age_bacon
    age.bins <- seq(0,10000,100)
    age.index <- as.matrix(as.numeric(cut(sample.ages,breaks = age.bins,labels=seq(1:(length(age.bins)-1)))))
    
par(mfrow=c(4,4))
   for(age in age.index){
   	if(age < 100){
   		plot(vals,exp(outLik.save[[s]][,age] - max(outLik.save[[s]][,age]))/-sum(outLik.save[[s]][,age])
        ,main=paste( site_name ,', Age =',age),typ='l',ylab=NA)
        abline(v=biomassCI[[s]][2,age],col='red')
        abline(v=biomassCI[[s]][c(1,3),age])
   	}
        
    }
}

}
dev.off()

#load(file="biomass.CI.Rdata")

list.buddy <- seq(1,length(unique(x.meta[,1])),1)
lat.save <- list()
for(i in 1:length(unique(x.meta[,1]))){
	lat.save[[i]]<-unique(x.meta[x.meta[,1]==unique(x.meta[,1])[i],2])
}
#lat.save[[1]]<-45.78650
lat.save <- unlist(lat.save)

site.id.list <- unique(x.meta[,1])
dataset.id.list <- unique(x.meta[,4])
name.list <-unique(x.meta[,5])

pdf(paste0('all.sites.',Sys.Date(),'.pdf'))
for(i in 1:length(biomassCI)){
	if(length(biomassCI[[i]])>1){

      map('state', xlim=c(-97.3,-83), ylim=c(41.5,50))
      site_number = unique(x.meta[,1])[i]
      points(x.meta[x.meta$site.id==site_number,3],      
      x.meta[x.meta$site.id==site_number,2], pch=19, cex=1)  
      #title(unique(x.meta[,1])[i])
      title(x.meta[x.meta$site.id==site_number,'site.name'][1])
      
     fig.mat <- matrix(1,26,1)
     fig.mat[1:6,]<-1
     fig.mat[7:26,]<-seq(2,21,1)
     
     layout(fig.mat)
     
     par(mar = c(0, 4, 0, 3), oma = c(5, 0, 4, 0), tcl=.5)

      plot(seq(100,9900,100),biomassCI[[i]][2,],
      cex=.1,ylim=c(0,200),xlim=c(-10,10000),ylab="Biomass",xlab="Years BP",main=NA,, xaxt='n')
      
      title(x.meta[x.meta$site.id==site_number,'site.name'][1],outer=TRUE)
      axis(3,at=seq(0,10000,1000),labels=seq(0,10000,1000),padj=1)

      ciEnvelope(seq(100,9900,100),biomassCI[[i]][1,],biomassCI[[i]][3,],col="lightblue")

      points(seq(100,9900,100),biomassCI[[i]][2,],cex=.5,pch=16)
      points(seq(100,10000,100),inits.keep.approx[[i]],cex=.5,pch=16,col='blue')
      points(seq(100,10000,100),inits.keep[[i]],cex=.5,pch=16,col='red')
      
      if(any(cast.x$site.id==site_number)){
      	points(0,cast.x[cast.x$site.id==site_number,ncol(cast.x)],col='purple',cex=2,pch=19)
      }
      
      keep.dataset.id <- unique(x.meta[x.meta$site.id==site_number,4])
      
      rug(x.meta[x.meta[,1]==site_number,]$age_bacon,lwd=2)
      rug(control.pts[which(control.pts[,1]%in%keep.dataset.id),]$geo_age,lwd=3,col="red")
      
       ten.count.use = ten.count[which(x.meta[,1]==site_number),]
      prop.use <- prop.table(as.matrix(ten.count.use),margin=1)    
     
      for(p in rev(match(names(sort(colMeans(prop.use))),colnames(prop.use)))){
        prop.plot<- cbind(as.vector(x.meta[which(x.meta$site.id==site_number),]$age_bacon),as.matrix(prop.use[,p]))      	
        prop.plot<-prop.plot[order(prop.plot[,1]),]
        plot(x=prop.plot[,1],y=prop.plot[,2],type="l",xlim=c(0,10000),ylim=c(0,max(prop.use[,p])),ylab=NA,yaxt='n', xaxt='n')
      	 #axis(2,at=signif(max(prop.use)/2),labels=signif(max(prop.use[,p])/2,digits=1))
      	 ciEnvelope(prop.plot[,1],rep(0,nrow(prop.use)),prop.plot[,2],col="darkblue")
      	 legend('topright',colnames(prop.use)[p],bty="n")
      	 legend('topleft',paste(signif(mean(prop.use[,p]),digits=2)),bty="n")
      } 

		
	}
      
      
}
dev.off()

####
#### PLOT ALL SITES ####
####
control.pts<-read.csv(text=getURL('https://raw.githubusercontent.com/PalEON-Project/stepps-baconizing/master/data/bacon_age_markers_v1.csv'))

plot.which=numeric(500)

for(i in list.buddy[order(lat.save)]){#
	if(length(biomassCI[[i]])>1){
    plot.which[i] <- i
    }
}

plot.which.keep<-plot.which[order(lat.save)]
plot.which.keep<-plot.which.keep[plot.which.keep!=0]

      
pdf("allpreds.max.lik.inits1.pdf")
for(i in plot.which.keep){
	  # par(mfrow=c(1,1))
	  # map('state', xlim=c(-97.3,-83), ylim=c(41.5,50))
      # points(x.meta[x.meta[,1]==unique(x.meta[,1])[i],3],      
      # x.meta[x.meta[,1]==unique(x.meta[,1])[i],2], pch=19, cex=1)  
      # #title(unique(x.meta[,1])[i])	
      	# title(x.meta[x.meta[,1]== site.id.list[i],5][1])
      	
     fig.mat <- matrix(1,26,1)
     fig.mat[1:6,]<-1
     fig.mat[7:26,]<-seq(2,21,1)
     
     layout(fig.mat)
     
      par(mar = c(0, 4, 0, 3), oma = c(5, 0, 4, 0),tcl=.5)
      
	  ########## plot biomass time series
      plot(seq(100,9900,100),biomassCI[[i]][2,],cex=.1,ylim=c(0,200),xlim=c(-10,10000),ylab="Biomass",xlab="Years BP",main=NA,, xaxt='n')
      title(x.meta[x.meta[,1]== site.id.list[i],5][1],outer=TRUE)
      axis(3,at=seq(0,10000,1000),labels=seq(0,10000,1000),padj=1)

      ciEnvelope(seq(100,9900,100),biomassCI[[i]][1,],biomassCI[[i]][3,],col="lightblue")

      points(seq(100,9900,100),biomassCI[[i]][2,],cex=.5,pch=16)
      
      ########### plot biomass change time series
       # plot(seq(150,9850,100),biomassCI[[i]][2,1:98]-biomassCI[[i]][2,2:99],cex=.1,xlim=c(-10,10000),ylab="Biomass",xlab="Years BP",main=NA,, xaxt='n')
      # title(x.meta[x.meta[,1]== site.id.list[i],5][1],outer=TRUE)
     
      # ciEnvelope(seq(150,9850,100),biomassCI[[i]][1,1:98]-biomassCI[[i]][1,2:99],biomassCI[[i]][3,1:98]-biomassCI[[i]][3,2:99],col="lightblue")
      # abline(h=0,lty=2)

      # points(seq(150,9850,100),biomassCI[[i]][2,1:98]-biomassCI[[i]][2,2:99],cex=.5,pch=16)
      
      keep.dataset.id <- unique(x.meta[x.meta[,1]==site.id.list[i],4])
      
      rug(x.meta[x.meta[,1]== site.id.list[i],]$age_bacon,lwd=2)
      rug(control.pts[which(control.pts[,1]%in%keep.dataset.id),]$geo_age,lwd=3,col="red")

      ten.count.use = ten.count[which(x.meta$SiteID==site.id.list[i]),]
      prop.use <- prop.table(as.matrix(ten.count.use),margin=1)    
     
      for(p in rev(match(names(sort(colMeans(prop.use))),colnames(prop.use)))){
        prop.plot<- cbind(as.vector(x.meta[which(x.meta$site.id==site.id.list[i]),]$age_bacon),as.matrix(prop.use[,p]))      	
        prop.plot<-prop.plot[order(prop.plot[,1]),]
        plot(x=prop.plot[,1],y=prop.plot[,2],type="l",xlim=c(0,10000),ylim=c(0,max(prop.use[,p])),ylab=NA,yaxt='n', xaxt='n')
      	 #axis(2,at=signif(max(prop.use)/2),labels=signif(max(prop.use[,p])/2,digits=1))
      	 ciEnvelope(prop.plot[,1],rep(0,nrow(prop.use)),prop.plot[,2],col="darkblue")
      	 legend('topright',colnames(prop.use)[p],bty="n")
      	 legend('topleft',paste(signif(mean(prop.use[,p]),digits=2)),bty="n")
      } 

    }
dev.off()



name.vec <- seq(100,10000,100)
only.means<-laply(biomassCI,function(x){return(x[seq(2,299,3)])})
get.lat.long<-matrix(0,183,2)
for(i in 2:183){
	get.lat.long[i,1]<- unique(x.meta[x.meta[,1]==unique(x.meta[,1])[i],c(3)])
	get.lat.long[i,2]<- unique(x.meta[x.meta[,1]==unique(x.meta[,1])[i],c(2)])
}


only.means[only.means>200] <- NA

breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)))

pdf('pred.points.map.pdf')
for(r in 1:99){

data_binned <-  cut(only.means[,r], c(breaks), include.lowest = FALSE, labels = FALSE)

map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(get.lat.long[which(only.means[,r]>0),1],get.lat.long[which(only.means[,r]>0),2], pch=21,
		cex=1.1, bg=colors[data_binned[which(only.means[,r]>0)]],lwd=.2)
plotInset(-90,47,-82.5,50,
          expr={
          	keep.col<-unique(data_binned)
          	keep.col<-keep.col[!is.na(keep.col)]
          	keep.col<-sort(keep.col)
          	is.na.vec <- rep(NA,10)
          	is.na.vec[keep.col]<-colors[keep.col]
          	
          	hist(data_binned,col=is.na.vec,xaxt="n",xlab=NA,
          	ylab=NA,main=NA,cex.lab=.5,cex.axis=.5,
          	xlim=c(0,length(breaks)),ylim=c(0,20),breaks=seq(0,12,1))
          	
          	axis(side=1,breaks,at = seq(0,11,1),cex.axis = .5,las=2,line=0)
          	mtext(side = 1, "Biomass (Mg/ha)", line = 1.5,cex=.5)
         mtext(side = 2, "Frequency", line = 1.7,cex=.5)
          	})
title(paste("Biomass @",name.vec[r]))
}
dev.off()

only.means<-laply(biomassCI,function(x){return(x[seq(2,299,3)])})
only.means1<-only.means[-which(is.na(only.means)),1:99]
get.lat.long1<-get.lat.long[-which(is.na(only.means)),]


keep.max <- list()
for(i in plot.which.keep){    
    keep.max[[i]] <- which((biomassCI[[i]][2,1:98]/biomassCI[[i]][2,2:99]) == min(biomassCI[[i]][2,1:98]/biomassCI[[i]][2,2:99]))
    }

max.breaks <-  seq(0,10000,1000)
max.colors.ramp <- colorRampPalette(c("red","yellow","green","purple","blue"))
max.colors<-max.colors.ramp(length(max.breaks)-1)
data_binned <-  cut(unlist(keep.max)*100, c(max.breaks), include.lowest = FALSE, labels = FALSE)
map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(get.lat.long[plot.which.keep,1],get.lat.long[plot.which.keep,2], pch=21,
		cex=1.1, bg=max.colors[data_binned],lwd=.2)

plotInset(-90,47,-82.5,50,
          expr={
          	keep.col<-unique(data_binned)
          	keep.col<-keep.col[!is.na(keep.col)]
          	keep.col<-sort(keep.col)
          	is.na.vec <- rep(NA,10)
          	is.na.vec[keep.col]<-max.colors[keep.col]
          	
          	hist(data_binned,col=is.na.vec,xaxt="n",xlab=NA,
          	ylab=NA,main=NA,cex.lab=.5,cex.axis=.5,
          	xlim=c(0,length(breaks)),ylim=c(0,20),breaks=seq(0,12,1))
          	
          	axis(side=1,max.breaks,at = seq(0,10,1),cex.axis = .5,las=2,line=0)
          	mtext(side = 1, "Time of Big Change", line = 1.5,cex=.5)
         mtext(side = 2, "Frequency", line = 1.7,cex=.5)
          	})

hist(unlist(keep.max)*100,col=max.colors)

prop.all <- prop.table(as.matrix(ten.count),margin=1)   
par(mfrow=c(4,4))
for(s in 1:20){
	plot(x.meta$age_bacon,prop.all[,s])
	title(colnames(prop.all)[s])
}  


diff.breaks <- c(seq(0,.99,length.out=5),seq(1.01,2,length.out=5))
diff.colors.ramp <- colorRampPalette(c("darkblue",'lightblue',"white","pink","red"))
diff.colors <- diff.colors.ramp(length(diff.breaks)-1)

pdf('diff.points.map.ratio.pdf')
for(r in rev(seq(1,98,1))){

diff.create <- only.means1[,r]/only.means1[,r+1] #first one is the youngest
diff.create[is.na(diff.create)] <- 0 
diff.create[diff.create<diff.breaks[1]]<-diff.breaks[1]
diff.create[diff.create>diff.breaks[length(diff.breaks)]]<-diff.breaks[length(diff.breaks)]
data_binned <-  cut(diff.create, c(diff.breaks), include.lowest = TRUE, labels = FALSE)
#data_binned[is.na(data_binned)]<-5
#data_binned<-as.numeric(data_binned)

map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(get.lat.long1[,1],get.lat.long1[,2], pch=21,
		cex=1.1, bg=diff.colors[data_binned],lwd=.2)
plotInset(-90,47,-82.5,50,
          expr={
          	
          	keep.col<-unique(data_binned)
          	keep.col<-keep.col[!is.na(keep.col)]
          	keep.col<-sort(keep.col)
          	is.na.vec <- rep(NA,10)
          	is.na.vec[keep.col]<-diff.colors[keep.col]
          	
          	hist(data_binned,col=is.na.vec,xaxt="n",xlab=NA,
          	ylab=NA,main=NA,cex.lab=.5,cex.axis=.5,
          	xlim=c(0,length(diff.breaks)),ylim=c(0,40),
          	breaks=seq(0,10,1))
          	
          	axis(side=1,signif(diff.breaks),at = seq(0,9,1),cex.axis = .5,las=2,line=0)
          	mtext(side = 1, "Biomass Diff", line = 1.5,cex=.5)
            mtext(side = 2, "Frequency", line = 1.7,cex=.5)
         
          	})
title(paste("Biomass @",name.vec[r],"-",name.vec[r+1]))
}
dev.off()







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