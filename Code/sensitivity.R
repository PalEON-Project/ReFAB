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

par(mfrow=c(3,3))
do_sites = c(1988,269,1593,1820,2933,8574) #1101,318
quartz()
par(mfrow=c(3,3))
for(r in 1:length(do_sites)){  
for(f in 1:length(colnames(counts))){
  rows_use = which(x[,1]==do_sites[r])

  counts <- t(as.matrix(ten.count[rows_use[3],]))

  pick_spp = colnames(counts)[f]
  
  if(counts[,pick_spp]>0){
    up.seq = seq(1,100,5)
    max_neg = counts[,pick_spp]
    down.seq = seq(1,max_neg,5)
    for(i in 1:length(up.seq)){
      counts = rbind(counts, c(counts[1,colnames(counts)==pick_spp]+up.seq[i],counts[1,colnames(counts)!=pick_spp]))
    }
    for(i in 1:length(down.seq)){
      counts = rbind(counts, c(counts[1,colnames(counts)==pick_spp]-i,counts[1,colnames(counts)!=pick_spp]))
    }
  } else {
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

#pdf(paste0("site_",do_sites[r],"_",pick_spp,"_sensitivity.pdf") )

plot(props_vec[,1],biomass_vec,pch=19,cex=1,ylim=c(0,400),
     main = paste0(pick_spp,"- site ",do_sites[r]),
     xlab = "Pollen Prop",ylab="Biomass")
#points(props_vec[1,1],biomass_vec[1],pch=19,cex=1,col = "red")
#dev.off()

#save(to_save,file = paste0("site_",do_sites[r],pick_spp,"sens.RData"))
}
}

site=318
load("site_1988ABIESsens.RData")
#pdf(paste0("site_",do_sites[f],"_sensitivity.pdf") )

#pdf(paste0("site_",site,"_sensitivity.pdf") )
seq1 = seq(1,(nrow(props_vec)),4)
seq2 = seq(1,18,1)
#par(mfrow = c(3,3)) 
quartz()
par(mfrow=c(3,3))
    for(i in 1:18){
      plot(props_vec[seq1[i]:(seq1[i]+3),i],biomass_vec[seq1[i]:(seq1[i]+3)],pch = 19, cex = 1.5, main = colnames(ten.count)[i],xlim=c(0,1),ylim = c(0,400))
    } 
#dev.off()  


    }








save_list1 = cbind(props_vec[,1],biomass_vec)

save_list = cbind(save_list[,1:5],biomass_vec,props_vec[,pick_col])

colnames(save_list) = c("NA","biom_prairie","prop_prairie","biom_oak",
                        "prop_oak","biom_birch","prop_birch")

quartz()
par(mfrow=c(1,3))
colors = c(rep("black",3),rep("blue",3),"black","blue")
plot(save_list[,3],save_list[,2],pch = 19, cex = 1.5, main = "Prairie - Site 1988",xlim=c(0,1),ylim = c(0,400),col = colors)
plot(save_list[,5],save_list[,4],pch = 19, cex = 1.5, main = "Oak - Site 1988",xlim=c(0,1),ylim = c(0,400), col = colors)
plot(save_list[,7],save_list[,6],pch = 19, cex = 1.5, main = "Birch - Site 1988",xlim=c(0,1),ylim = c(0,400),col = colors)



