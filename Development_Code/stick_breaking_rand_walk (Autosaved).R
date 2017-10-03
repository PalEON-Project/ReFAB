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
  
  sigma ~ dunif(0,5) #GELMAN PAPER
  
  b[1,1] ~ dunif(0,145)
  
  for(t in 2:T){
    b[1,t] ~ T(dlnorm(log(b[1,t-1]),1/sigma^2),0,145)
  }
  
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

# plot_sites = c(2550,2548,1991,7538,9953,10211,7542,1979,10132,7531,292,1701,1652,10149,269,29,824,1772,1901,9870,9876,9874,2235,8552,9826,238,2519,1771,1119,1814,1590,2933,10169,1891,10156,9873,8551,1729,9825,9867,9871,245,1548,1483,2524,1128,34,39,9768,10240,26,59,1815,318,2293,8572,23,1593,531,2231,1645,2534,1,8567,1598,1101,8573,24,2937,666,196,8574,25,10160,1679,680,10036,10034,1944,325,2956,1411,9738,9993,2864,10239,10242,10251,2239,661,518,10163,1541,8559,360)

samples.pred.save = list()
biomassCI=list()

do.sites<-c('Radtke Lake','Gass Lake','Kellners Lake','Cub Lake','Green Lake','Chippewa Bog')

for(s in 60:length(unique(x.meta[,1]))){ 
  print(paste('working on number',s,'of 183',(s/183)*100,'% complete'))
  #site_number = unique(x.meta[x.meta$sitename==do.sites[s],1])
  #site_number = unique(x.meta[x.meta$sitename=='GassLake',1])
  #site_number = unique(x.meta[x.meta$site.name=='Kellys Hollow',1])
  #site_number = unique(x.meta[x.meta$site.name=='Chippewa Bog',1])
  #site_number = unique(x.meta[x.meta$site.name==do.sites[p],1])
  #s = which(unique(x.meta[,1])==site_number)
  site_number = unique(x.meta[,1])[s]
  
  site_number = unique(x.meta[x.meta$site.name=='Kirchner Marsh',1])

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
    #beta.est = matrix(colMeans(samples1[100:nrow(samples1),]),K,I) 
    new.biomass = seq(1,200,1)
    Z.new = matrix(0,nrow=length(new.biomass),ncol=K)
    
    u<-c(rep(attr(Z.knots,"Boundary.knots")[1],1),attr(Z.knots,"knots"),rep(attr(Z.knots,"Boundary.knots")[2],1))
    
    data.pred = list(Y = Y)
    
    constants.pred = list(beta = beta1.est.real, beta1 = beta2.est.real, I = I, J = J,
                          T = T, n = n, u = u, N0 = rep(0, (length(u)-1)), 
                          N1 = rep(0, (length(u))), N2 = rep(0, (length(u)+1)),
                          N3 = rep(0, (length(u)+2)), age.index = age.index)
    
    inits.pred = list(b = matrix(100,1,T), sigma = 4.5) #
    
    dimensions.pred = list(exp.phi = c(T,I), exp.phi1 = c(T,I), phi.first = c(T,I),
                           phi.first1 = c(T,I), Zb = dim(Zb), Y = dim(Y),
                           p.rel = c(T,I), p.true = c(T,I)) #  b = dim(inits.pred$b),
    
    set.seed(0)
    
    model_pred <- nimbleModel(pred_code, inits = inits.pred, constants = constants.pred,
                              data = data.pred, dimensions = dimensions.pred)
    
    cm <- compileNimble(model_pred)
    
    spec.pred <- configureMCMC(model_pred, thin = 10, print = FALSE, useConjugacy = FALSE,
                               control = list(log=TRUE))
    spec.pred$getSamplers()
    spec.pred$addMonitors(c('p.true')) 
    Rmcmc.pred <- buildMCMC(spec.pred)
    
    Cmcmc.pred <- compileNimble(Rmcmc.pred, project = model_pred) 
    Cmcmc.pred$run(5000)
    
   calclik <- nimbleFunction(
    setup = function(model) {
    },
    run = function(biomasses=double(1), nSims=double()) {
        n <- length(biomasses)
        returnType(double(1))
        out = numeric(n)
        for(i in 1:n) {
           #set.seed(0)
           
           	model$b[1, 6] <- biomasses[i]
                     
            out[i] = 0
            for(j in 1:nSims) {
                model$simulate()  # simulate p.rel and fill in resulting p.true
                model$calculate('Y') # calculate logdensity of 'y'
                out[i] <- out[i] + exp(model$getLogProb('Y'))  # add density values across iterations
            }
            out[i] = out[i]/nSims
        }
        return(out)
    })

rcalclik <- calclik(model_pred)
ccalclik <- compileNimble(rcalclik, project = model_pred)
bvals <- seq(1,150, length = 50)
ccalclik$run(bvals, 5000) 
    
    samples.pred <-as.matrix(Cmcmc.pred$mvSamples)
    #two time points high and low
    #finer grid scale
    keep.prob <-matrix(NA,500,length(age.index))
    for(i in 1:500){
    	one.init <- matrix(as.matrix(Cmcmc.pred$mvSamples)[i,2101:4200],T,I)
    	for(a in 1:length(age.index)){
    		set.seed(0)
    		keep.prob[i,a] <- dmulti(x = Y[a,], prob = one.init[age.index[a],], size = n[a])
    	}
    }
    
    keep.prob125 <- keep.prob
    pdf('k.marsh.liks.pdf')
    par(mfrow=c(4,4))
    for(i in 1:length(age.index)){
    	plot(c(10,25,50,100,125),c(colMeans(keep.prob10)[i],colMeans(keep.prob25)[i],colMeans(keep.prob50)[i],colMeans(keep.prob100)[i],colMeans(keep.prob125)[i]),typ='l',main=age.index[i],xlim=c(0,150),ylab=NA)
    }
    dev.off()
    ##validate by looking at predictions and seeing if the plots match up with the same estimates.
    
   par(mfrow=c(3,3))
   for(j in 1:J){
        plot(vals,exp(outLik[,j] - max(outLik[,j]))/-sum(outLik[,j])
        ,typ='l',ylab=NA,main=w[j])
    }
   
    samples.pred.save[[s]]<-samples.pred
    
    # rug(sample.ages/100)
    dyn.unload(model_pred$nimbleProject$cppProjects[[1]]$getSOName())
    biomassCI[[s]] <- apply(samples.pred[100:500,1:99],2,quantile,c(0.025,0.5,0.975))
    #samples.pred.good.inits.save[[s]] <- samples.pred.good.inits
    
  }
  save(biomassCI, file="biomass.CI10.Rdata")
  save(samples.pred.save, file="samples.pred.save2.Rdata")
}

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

control.pts<-read.csv(text=getURL('https://raw.githubusercontent.com/PalEON-Project/stepps-baconizing/master/data/bacon_age_markers_v1.csv'))

plot.which=numeric(500)

for(i in list.buddy[order(lat.save)]){#
  if(length(biomassCI[[i]])>1){
    plot.which[i] <- i
  }
}

plot.which.keep<-plot.which[order(lat.save)]
plot.which.keep<-plot.which.keep[plot.which.keep!=0]


pdf(paste0('all.sites.',Sys.Date(),'.pdf'))
for(i in plot.which.keep){
  if(length(biomassCI[[i]])>1){
    
    map('state', xlim=c(-97.3,-83), ylim=c(41.5,50))
    site_number = unique(x.meta[,1])[i]
    points(x.meta[x.meta$site.id==site_number,3],      
           x.meta[x.meta$site.id==site_number,2], pch=19, cex=1)  
    #title(unique(x.meta[,1])[i])
    title(x.meta[x.meta$site.id==site_number,'site.name'][1])
    
    fig.mat <- matrix(1,27,1)
    fig.mat[1:6,]<-1
    fig.mat[7:27,]<-seq(2,22,1)
    
    layout(fig.mat)
    par(mar = c(0, 4, 0, 3), oma = c(5, 0, 4, 0),tcl=.5)
    
    breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)))


  data_binned <-  cut(rev(biomassCI[[i]][2,]), c(breaks), include.lowest = FALSE, labels = FALSE)
    
     plot(seq(100,9900,100),rev(biomassCI[[i]][2,]),
         cex=.1,ylim=c(0,150),xlim=c(-10,10000),ylab="Biomass (Mg / Ha)",
         xlab="Years Before Present",main=NA,, xaxt='n')
    
    title(x.meta[x.meta$site.id==site_number,'site.name'][1],outer=TRUE)
    axis(3,at=seq(0,10000,1000),labels=rev(seq(0,10000,1000)),padj=1)
    
    ciEnvelope(seq(100,9900,100),rev(biomassCI[[i]][1,]),rev(biomassCI[[i]][3,]),col="gray")
    
    points(seq(100,9900,100),rev(biomassCI[[i]][2,]),cex=.8,pch=16,col = colors[data_binned])
        
    keep.dataset.id <- unique(x.meta[x.meta$site.id==site_number,4])
    
    rug(x.meta[x.meta[,1]==site_number,]$age_bacon,lwd=2)
    rug(control.pts[which(control.pts[,1]%in%keep.dataset.id),]$geo_age,lwd=3,col="red")
    
 ten.count.use = ten.count[which(x.meta[,1]==site_number),]
    prop.use <- prop.table(as.matrix(ten.count.use),margin=1)    
    
    for(p in rev(match(names(sort(colMeans(prop.use))),colnames(prop.use)))){
      prop.plot<- cbind(as.vector(x.meta[which(x.meta$site.id==site_number),]$age_bacon),as.matrix(prop.use[,p]))      	
      prop.plot<-prop.plot[order(prop.plot[,1]),]
      plot(x=rev(prop.plot[,1]),y=prop.plot[,2],type="l",xlim=c(-10,10000),ylim=c(0,max(prop.use[,p])),ylab=NA,yaxt='n', xaxt='n')
      #axis(2,at=signif(max(prop.use)/2),labels=signif(max(prop.use[,p])/2,digits=1))
      ciEnvelope(rev(prop.plot[,1]),rep(0,nrow(prop.use)),prop.plot[,2],col="darkblue")
      legend('topleft',colnames(prop.use)[p])
      #legend('topleft',paste(signif(mean(prop.use[,p]),digits=2)),bty="n")
    } 
    
    
  }
  
  
}
dev.off()

####
#### PLOT ALL SITES ####
####



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
only.means.all<-lapply(biomassCI,function(x){return(x[seq(2,299,3)])})
get.lat.long<-matrix(0,183,2)
for(i in 2:183){
  get.lat.long[i,1]<- unique(x.meta[x.meta[,1]==unique(x.meta[,1])[i],c(3)])
  get.lat.long[i,2]<- unique(x.meta[x.meta[,1]==unique(x.meta[,1])[i],c(2)])
}


only.means[only.means>200] <- NA

breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)))

pdf(paste0('pred.points.map.fewer',Sys.Date(),'.pdf'))
for(r in seq(1,99,length.out=15)){
  
  only.means <- unlist(lapply(only.means.all,function(x){return(x[r])}))
  
  data_binned <-  cut(only.means, c(breaks), include.lowest = FALSE, labels = FALSE)
  
  long.keep <- list()
  lat.keep <- list()
  for(i in 1:length(only.means)){
  	long.keep[[i]] <- x.meta[x.meta[,1]==names(only.means)[i],'long'][1]
  	lat.keep[[i]] <- x.meta[x.meta[,1]==names(only.means)[i],'lat'][1]
  }
  
  
  map('state', xlim=c(-98,-81), ylim=c(41.5,50))
  points(unlist(long.keep),unlist(lat.keep), pch=21,
         cex=1.1, bg=colors[data_binned],lwd=.2)
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

only.means<-lapply(biomassCI,function(x){return(x[seq(2,299,3)])})
only.means1<-only.means[-which(is.na(only.means)),1:99]
get.lat.long1<-get.lat.long[-which(is.na(only.means)),]


keep.max <- list()
for(i in plot.which.keep){    
  keep.max[[i]] <- which((biomassCI[[i]][2,1:98]/biomassCI[[i]][2,2:99]) == min(biomassCI[[i]][2,1:98]/biomassCI[[i]][2,2:99]))
}

max.breaks <-  seq(0,10000,1000)
max.colors.ramp <- colorRampPalette(c("yellow","blue"))
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