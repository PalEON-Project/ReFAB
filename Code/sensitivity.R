####
#### Prediction Models
####

load("all.preds.min.list.RData")

x = pol.cal.count[pol.cal.count$Age>=200,]
x = x[x$Age<=2000,]

x = x[,-which(colnames(x)==c("PINUSX"))]
trees <- c("ACERX","CUPRESSA","FRAXINUX","FAGUS","CYPERACE","LARIXPSEU","TSUGAX","QUERCUS","TILIA","BETULA","PICEAX","OSTRYCAR","ULMUS","ABIES","POPULUS")
other.trees <- c("TAXUS","NYSSA","JUGLANSX","CASTANEA","PLATANUS","SALIX","LIQUIDAM","ALNUSX")
ten.count = matrix(0,nrow(x),length(trees)+3)
prairie <- c("CORYLUS","ARTEMISIA","ASTERX","POACEAE","AMBROSIA","CHENOAMX")
ten.count[,1] <- rowSums(x[,prairie])
ten.count[,2] <- rowSums(x[,other.trees])
ten.count[,3:(length(trees)+2)] <- as.matrix(x[,trees])
ten.count[,(length(trees)+3)] <- rowSums(x[,7:ncol(x)]) - rowSums(ten.count)
colnames(ten.count)<-c("prairie","other trees",trees,"other herbs")

ten.count <- round(ten.count)

#par(mfrow=c(3,3))
do_sites = c(269,1988) #1101,318
quartz()
par(mfrow=c(3,3))
for(r in 1:length(do_sites)){  
for(f in 1:18){
  rows_use = which(x[,1]==do_sites[r])

  counts <- t(as.matrix(ten.count[rows_use[c(3)],]))

  pick_spp = colnames(counts)[f]
  
  base_num = counts[,pick_spp]
  
  min.calc1 = min(ten.count[rows_use,f])
  max.calc1 = max(ten.count[rows_use,f])
  
  if(base_num>30){
    up.seq = seq(1,max.calc1 + 10,20)
    down.seq = seq(0,base_num,10)
    for(i in 1:length(up.seq)){
      counts = rbind(counts, c(counts[1,colnames(counts)==pick_spp]+up.seq[i],counts[1,colnames(counts)!=pick_spp]))
    }
    for(i in 1:length(down.seq)){
      counts = rbind(counts, c(counts[1,colnames(counts)==pick_spp]-down.seq[i],counts[1,colnames(counts)!=pick_spp]))
    }
  } 
  
  if(base_num>0 & base_num<30){
    up.seq = seq(1, max.calc1 + 5,2)
    down.seq = seq(0, base_num,2)
    for(i in 1:length(up.seq)){
      counts = rbind(counts, c(counts[1,colnames(counts)==pick_spp]+up.seq[i],counts[1,colnames(counts)!=pick_spp]))
    }
    for(i in 1:length(down.seq)){
      counts = rbind(counts, c(counts[1,colnames(counts)==pick_spp]-down.seq[i],counts[1,colnames(counts)!=pick_spp]))
    }
  } 
  
  if(base_num==0){
    up.seq1 = seq(1,10,1)
    for(i in 1:length(up.seq1)){
      counts = rbind(counts, c(counts[1,colnames(counts)==pick_spp]+up.seq1[i],counts[1,colnames(counts)!=pick_spp]))
    }
  }

  if(f == 18){
    counts = counts[-1,c(seq(2,f,1),1)]
  }
  if(f == 1){
    counts = counts
  }
  if(f != 1 & f!=18){
  counts = counts[-1,c(seq(2,f,1),1,seq((f+1),18,1))]
  }

J = nrow(counts)
Zb = matrix(NA,J,ncol(Z.knots))
DFS = ncol(Zb)
p = matrix(NA,J,ncol(counts)); phi.first = p; phi = p
beta.est.real = matrix(summary(csamp.real.cal)$statistics[,1],ncol(Z.knots),ncol(Y))
new.biomass = seq(1,400,1)
Z.new = bs(new.biomass,intercept=TRUE,knots = attr(Z.knots,"knots"),Boundary.knots = attr(Z.knots,"Boundary.knots"))

data.real.pred = list("I" = ncol(counts), "DFS" = DFS, "J" = J, "Zb" = Zb, "Y" = counts , "n" = rowSums(counts), "Z" = Z.new, "beta" = beta.est.real, "p" = p, "phi.first" = phi.first)

inits.pred = list(list(b = rep(20,J)),list(b=rep(300,J)))

mod.real.pred <- jags.model(paste0(model.dir,'biomass_pred_jags.R'),data = data.real.pred, n.chains = length(inits.pred), n.adapt = n.adapt, inits = inits.pred)
csamp.real.pred <- coda.samples(mod.real.pred,c("b"),n.iter=n.iter)

biomass_vec = summary(csamp.real.pred)$statistics[,1]

props_vec = counts/rowSums(counts)

to_save = cbind(props_vec,biomass_vec)

print(paste(pick_spp,do_sites[r]))

prop_all_age = ten.count[rows_use,]/rowSums(ten.count[rows_use,])

min.calc = min(prop_all_age[,f])
max.calc = max(prop_all_age[,f])

#pdf(paste0("site_",do_sites[r],"_",pick_spp,"_sensitivity.pdf") )

plot(props_vec[,f],biomass_vec,pch=19,cex=1,ylim=c(0,400),
     main = paste0(pick_spp,"- site ",do_sites[r]),
     xlab = "Pollen Prop",ylab="Biomass",xlim = c(0,max(props_vec[,f])))
abline(v = min.calc, col = "red",lwd = 2)
abline(v = max.calc, col = "red",lwd = 2)
axis(3,at = props_vec[,f], labels = counts[,f])
#dev.off()
}
}