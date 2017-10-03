model_pred <- nimbleModel(pred_code, inits = inits.pred, constants = constants.pred,
                          data = data.pred, dimensions = dimensions.pred)
spec.pred <- configureMCMC(model_pred, thin = 10, print = FALSE,
                           useConjugacy = FALSE,control = list(log=TRUE))#,control = list(log=TRUE)
spec.pred$removeSamplers('b')
for(i in 1:T){
  spec.pred$addSampler(paste0('b[1,',i,']'),'slice')
}
spec.pred$addMonitors(c("b")) 
Rmcmc.pred <- buildMCMC(spec.pred)
cm <- compileNimble(model_pred, Rmcmc.pred)

samplesList <- runMCMC(mcmc = cm$Rmcmc.pred, niter = 30000, nchains = 3,
        inits = list(list(b = matrix((25),1,T), sigma = 100),
                     list(b = matrix((100),1,T), sigma = 50),
                     list(b = matrix((140),1,T), sigma = 150)))

quartz()
pdf('cub.lake.long.run.trace.slice.pdf')
par(mfrow=c(3,3))
for(i in 1:100){
  plot(samplesList[[1]][,i],typ='l',main=i,ylim=c(0,150),col=rainbow(3,alpha = 1)[1])
  points(samplesList[[2]][,i],typ='l',col=rainbow(3,alpha = 0.6)[2])
  points(samplesList[[3]][,i],typ='l',col=rainbow(3,alpha = 0.6)[3])
}
for(i in 101){
  plot(samplesList[[1]][,i],typ='l',main='SIGMA',
       ylim=c(0,max(c(samplesList[[1]][,i],samplesList[[2]][,i],samplesList[[3]][,i]))),col=rainbow(3,alpha = 1)[1])
  points(samplesList[[2]][,i],typ='l',col=rainbow(3,alpha = 0.6)[2])
  points(samplesList[[3]][,i],typ='l',col=rainbow(3,alpha = 0.6)[3])
}
dev.off()

biomassCI1 <-  apply(as.data.frame(rbind(samplesList[[1]][,1:100],
                                         samplesList[[2]][,1:100],
                                         samplesList[[3]][,1:100])),
                     2,quantile,c(0.025,0.5,0.975))

quartz()
plot_biomass_ts(site_number = site_number, biomassCI = biomassCI1)

### Make plots of ind. cub lake ts

### Work on adding more knots



