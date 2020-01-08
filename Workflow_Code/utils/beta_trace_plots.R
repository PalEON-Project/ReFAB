
path_to_betas <- '~/Downloads/betas_FULL/'

samples.mixed.keep <- list()

for(i in 1:50){
  load(paste0(path_to_betas,list.files(path_to_betas)[i]))
  samples.mixed.keep[[i]] <- samples.mixed[49900:50000,]
  
}

samps.all <- do.call(rbind,samples.mixed.keep)

pdf('[10]beta_trace_plots.pdf')
par(mfrow=c(3,3))
set.seed(2)
for(i in sample(size=9,x=1:220)){
  plot(samps.all[,i],type = 'l',main=paste('Trace Plot for',colnames(samples.mixed)[i]),
       ylab = colnames(samples.mixed)[i],
       xlab = 'MCMC Iteration')
}
dev.off()
