#####
##### Trace Plots #####
#####

for(t in 1:100){
plot(samplesList[[1]][,t],typ='l',main=t,ylim=c(0,150),col=rainbow(3,alpha = 1)[1])
points(samplesList[[2]][,t],typ='l',col=rainbow(3,alpha = 0.6)[2])
points(samplesList[[3]][,t],typ='l',col=rainbow(3,alpha = 0.6)[3])

}

for(i in 101){
  plot(samplesList[[1]][,i],typ='l',main='SIGMA',
       ylim=c(0,max(c(samplesList[[1]][,i],samplesList[[2]][,i],samplesList[[3]][,i]))),col=rainbow(3,alpha = 1)[1])
  points(samplesList[[2]][,i],typ='l',col=rainbow(3,alpha = 0.6)[2])
  points(samplesList[[3]][,i],typ='l',col=rainbow(3,alpha = 0.6)[3])
}