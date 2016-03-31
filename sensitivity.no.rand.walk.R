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

   	if(validation == TRUE){
   		counts = round(plot_biomass_pollen[,4:(ncol(plot_biomass_pollen)-2)])
   		colnames(counts) <- colnames(plot_biomass_pollen[,4:(ncol(plot_biomass_pollen)-2)])
   		counts = counts[,-which(colnames(counts)==c("PINUSX"))]
   		trees <-c("ALNUSX","JUGLANSX","ACERX","CUPRESSA","FRAXINUX","FAGUS","CYPERACE","LARIXPSEU","TSUGAX","QUERCUS","TILIA","BETULA","PICEAX","OSTRYCAR","ULMUS","ABIES","POPULUS")
   		other.trees <- c("TAXUS","NYSSA","CASTANEA","PLATANUS","SALIX","LIQUIDAM")
   		ten.count = matrix(0,nrow(counts),length(trees)+3)
   		prairie <- c("ARTEMISIA","ASTERX","POACEAE","AMBROSIA","CHENOAMX","CORYLUS")
   		ten.count[,1] <- rowSums(counts[,prairie])
   		ten.count[,2] <- rowSums(counts[,other.trees])
   		ten.count[,3:(length(trees)+2)] <- counts[,trees]
   		ten.count[,(length(trees)+3)] <- rowSums(counts) - rowSums(ten.count)
   		colnames(ten.count)<-c("prairie","other trees",trees,"other herbs")
   		
   		biomass = plot_biomass_pollen[,1]
   		
   		Y=ten.count
   		counts = ten.count
   	}

ten.count.save = ten.count
ten.count = round(ten.count.save)
counts.save <- counts
ten.count.add.save = array(0,dim=c(61,20,nrow(counts.save)))
count.vec <- seq(2,61,3)
slope.save <- matrix(0,ncol(counts.save),nrow(counts.save))
biomassCI.save <- list()
site.number.save <- list()

load("slope.save.Rdata")

for(r in 84:nrow(counts.save)){
#site_number = sample(x = unique(x.meta[,1]), size = 1) #1128#plot_sites[i]
site.number.save[[r]] <- r#site_number

#ten.count.use = ten.count[which(x.meta[,1]==site_number),]
ten.count.use = counts.save[r,]
Y = as.matrix(ten.count.use)

ten.count.pick <- ten.count.use#[sample(size=1,x=1:nrow(ten.count.use)),]
ten.count.add <- ten.count.pick

for(i in 1:length(ten.count.pick)){
	mat.num <- ten.count.pick[i]
	if(mat.num >= 2){
		add.num = mat.num + c(-2,2,4)
	}
	if(mat.num == 1){
		add.num = mat.num + c(-1,2,4)
	}
	if(mat.num == 0){
		add.num = c(1,3,5)
	}
	add.vecs = rbind(ten.count.pick,ten.count.pick,ten.count.pick)
	add.vecs[,i] = add.num
	ten.count.add = rbind(ten.count.add,add.vecs)
}

ten.count.add.save[,,r] <- ten.count.add
Y = ten.count.add
counts = Y

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
Cmcmc.pred$run(10000)
samples.pred <- as.matrix(Cmcmc.pred$mvSamples)
proc.time() - ptm

ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}
samples.pred = samples.pred[2000:10000,]
biomassCI.save[[r]] = apply(samples.pred,2,quantile,c(0.025,0.5,0.975))
biomassCI = apply(samples.pred,2,quantile,c(0.025,0.5,0.975))

	for(i in 1:ncol(ten.count.add)){
		slope.save[i,r] = lm(c(biomassCI[2,1],biomassCI[2,count.vec[i]:(count.vec[i]+2)]) ~ c(ten.count.add[1,i],ten.count.add[count.vec[i]:(count.vec[i]+2),i]))$coefficients[[2]]
    }
}
save(slope.save,file="slope.save.Rdata")
pdf("pollen_count_sensitivity_valid_sites.pdf")
par(mfrow=c(2,2))
for(s in 1:20){
	breaks <-  c(min(slope.save),seq(-2,-1,1),-.25,.25)
	breaks1 <-  c(seq(1,4,1),max(slope.save))
	colors <- c("darkblue","blue", "lightblue","antiquewhite3")
	colors1 <- c("pink","salmon", "red", "darkred","maroon")
	data_binned1 <-  cut(slope.save[s,], c(breaks,breaks1),
	 include.lowest = FALSE, labels = FALSE)
	 colors.comb <- c(colors,colors1)
	 
	map('state', xlim=c(-98,-81), ylim=c(41,50))
	for(i in 1:nrow(final_coors)){
		#points(x.meta[x.meta$SiteID==site.number.save[[i]],3],
		#x.meta[x.meta$SiteID==site.number.save[[i]],2], pch=19,
		# cex=1, col=colors.comb[data_binned1[i]])	
		points(final_coors[i,2],final_coors[i,1], pch=19,
		cex=1, col=colors.comb[data_binned1[i]])
		}
	title(colnames(ten.count)[s])
	if(s == c(0,0,3,0,0,6,0,0,9,0,0,12,0,0,15,0,0,18,0,20)[s]){
	plot.new()
	legend("center",legend=c("-57 - -2.01","-2 - -1.01","-1 - -.26","-.25 - .24",".25 - .99","1 - 1.99","2 - 2.99","3 - 3.99","4 - 36"),col=c(colors,colors1),pch=19)
	}	
}
dev.off()

hist(cut(slope.save, c(breaks,breaks1),
	 include.lowest = TRUE, labels = FALSE), col=colors.comb, breaks=seq(0,10,1),
	 include.lowest=FALSE, main="hist of slopes", xaxt='n')
	 legend("topright",legend=c("-57 - -2.01","-2 - -1.01","-1 - -.26","-.25 - .24",".25 - .99","1 - 1.99","2 - 2.99","3 - 3.99","4 - 36"),col=c(colors,colors1),pch=19)
	

if(prairie.sensitivity ==TRUE){
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
}

if(forest.sensitivity ==TRUE){
par(mfrow=c(4,5))
for(i in 1:20){
	plot(Y[Y[,i]!=forest.vec[i],i],colMeans(samples.pred)[Y[,i]!= forest.vec[i]],ylim=c(0,max(biomassCI)+1),pch=16,cex=1,xlab=paste("counts",names(prairie.vec)[i]),ylab="biomass")
	ciEnvelope(Y[Y[,i]!= forest.vec[i],i],biomassCI[1,Y[,i]!= forest.vec[i]],biomassCI[3,Y[,i]!= forest.vec[i]],col="lightblue")
	points(Y[Y[,i]!= forest.vec[i],i],colMeans(samples.pred)[Y[,i]!= forest.vec[i]],pch=16,cex=1)
	abline(h=mean(samples.pred[,1]),col="red")
	abline(h=biomassCI[1,1],col="red",lty=2)
	abline(h=biomassCI[3,1],col="red",lty=2)
	abline(v=Y[1,i],col="red")
}
}

if(validation==TRUE){
	
	reg <- lm(colMeans(samples.pred)[-sites_rm]~biomass[-sites_rm]+0)
	reg1 <- lm(colMeans(samples.pred)[sites_rm]~biomass[sites_rm]+0)
	plot(biomass[-sites_rm],colMeans(samples.pred)[-sites_rm],ylim=c(0,200),xlim=c(0,200),pch=16,cex=1,main="Validation",
	xlab="True Biomass",ylab="Estimated Biomass")
	points(biomass[sites_rm],colMeans(samples.pred)[sites_rm],col="red",pch=16,cex=1)
	abline(reg)
	abline(reg1,col="red")
	text(35,175,paste("R^2 calib sites",signif(summary(reg)[[8]],digits=3)))
	text(35,165,paste("R^2 left out sites",signif(summary(reg1)[[8]],digits=3)))
	legend("bottomright",c("calib sites","left out sites"),pch=c(16,16),col=c("black","red"))
	
	hist(colMeans(samples.pred)[-sites_rm],breaks=20,freq=FALSE,col="gray",main="Validation")
	lines(density(colMeans(samples.pred)[sites_rm]),lwd=2)
	legend('topright',c("left out sites"),lty=c(1),lwd=2)
	
}

