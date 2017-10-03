





# log posterior
#pdf(paste0(do.sites[s],'logpostlik.withlognorm.pdf'))
par(mfrow=c(4,4))
for(age in age.index){
plot(vals,outPost[,age],main=paste('Posts, Age =',age))	
#abline(v=biomassCI[,age])
}

par(mfrow=c(4,4))
for(age in age.index){
plot(vals,outLik[,age],main=paste('Liks, Age =',age))
#abline(v=biomassCI[2,age],col='red')
#abline(v=biomassCI[c(1,3),age])
}

# dev.off()



inits.keep<-matrix(NA,1,T)

for(s in 1:length(do.sites)){
    site_number <- unique(x.meta[x.meta$sitename==do.sites[s],1])
    sample.ages <- x.meta[x.meta[,1]==site_number,]$age_bacon
    age.bins <- seq(0,10000,100)
    age.index <- as.matrix(as.numeric(cut(sample.ages,breaks = age.bins,labels=seq(1:(length(age.bins)-1)))))
for(age in age.index){
			normliks <- exp(outLik[,age] - max(outLik[,age]))/-sum(outLik[,age])
    inits.keep[,age]<-which(normliks==max(normliks))
	}
}

init.mean<-rowMeans(inits.keep,na.rm=TRUE)

for(s in 1:length(do.sites)){
	inits.keep[s,is.na(inits.keep[s,])]<-init.mean[s]
	
}



outLik.save[[s]]<-outPost
outPost.save[[s]]<-outLik
# par(mfrow=c(4,2))
# for(i in 15:25) {
	# plot(samples.pred.save[[5]][,i],typ='l',ylab=i)
	# #abline(h=157)
# }
 
 par(mfrow=c(3,3))
 prop.stop<-prop.table(Y,margin=1)
 for(i in 1:20) { 
 	plot(sample.ages,prop.stop[,i],main=colnames(Y)[i])
 }
 
 biomassCI <- apply(samples.pred[,1:99],2,quantile,c(0.025,0.5,0.975))
 	 testCI<-apply(test.pred[,1:99],2,quantile,c(0.025,0.5,0.975))
plot(biomassCI[2,],typ='l',ylim=c(0,200))
lines(apply(samples.pred[500:2000,1:99],2,quantile,c(0.025,0.5,0.975))[1,])
lines(apply(samples.pred[500:2000,1:99],2,quantile,c(0.025,0.5,0.975))[3,])
