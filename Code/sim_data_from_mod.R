
sim <- matrix(0,nrow(counts),20)
for(i in 1:nrow(counts)){
	dirch.vec <- rdirichlet(1,exp(Z.knots%*%beta.est.real)[i,])
	sim[i,]<-rmultinom(1,prob = dirch.vec, size = rowSums(counts)[i])
}

sim.props <- prop.table(sim,1)

pdf('sim.props.pdf')
par(mfrow=c(3,3))
for(i in 1:20){
	plot(biomass,sim.props[,i],pch=19,main=colnames(counts)[i],ylim = range(c(sim.props[,i],counts[,i]/total_counts)))
	points(biomass,counts[,i]/total_counts,pch=19,col='red',cex= .7)
}
dev.off()
