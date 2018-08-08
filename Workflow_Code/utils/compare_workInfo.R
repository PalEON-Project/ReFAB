
compare_workInfo <- function(path_to_workInfo, locn){
  
  load('prediction.data_v2.Rdata')
  
  runnums <- dataID[dataID$name==locn,'ID']
  
  locnClean <- gsub(' ', '-', locn)
  
  out.keep <- list()
  for(i in 1:length(runnums)){
    load(paste0(path_to_workInfo,'workInfo_',runnums[i],'_',locnClean,'_Beta_',i,'.Rdata'))
    out.keep[[i]] <- out
  }
  
  minAge = 0
  maxAge = 10000
  ageInterval = 100
  
  site_number = unique(x.meta[x.meta$site.name == locn,1])
  x.meta.use <- x.meta[x.meta$site.name == locn,]
  
  source(file.path('Workflow_Code','utils','test_site.R'))
  test_site(x.meta.use)
  
  ten_count_use = ten.count[which(x.meta$site.name == locn), ]
  ten_count_use[which(is.na(ten_count_use))] <- 0
  Y = as.matrix(ten_count_use)
  
  sample_ages <- x.meta.use$age_bacon
  age_bins <- seq(minAge, maxAge, ageInterval)
  age_index <- as.matrix(as.numeric(
    cut(sample_ages, breaks = age_bins, labels=seq(1:(length(age_bins)-1)))
  ))
  
  tmp <- data.frame(cbind(age_index, Y))
  names(tmp)[1] <- 'age_index'
  
  Y2 <- aggregate(tmp, by = list(tmp$age_index), FUN = sum)
  dim(Y2)
  
  
  
  bMax <- 150
  par(mfrow=c(2,3))
  for(i in 1:(ncol(out)-1)){
    plot(as.numeric(Y2[i,3:24]/sum(Y2[i,3:24])),typ='o', ylab = 'Pollen Prop', xlab = NA, xaxt = 'n')
    points(as.numeric(Y2[i+1,3:24]/sum(Y2[i+1,3:24])),typ='o', ylab = 'Pollen Prop',
           xlab = NA,xaxt = 'n',col='red')
    axis(side = 1,at = 1:22,labels = colnames(Y),las=2,cex=1)
    
    plot(seq(5, bMax-5, by = 2),
         exp(out.keep[[1]][,i]-max(out.keep[[1]][,i]))/-sum(out.keep[[1]][,i]),
         typ='l',main=paste('age',Y2$Group.1[i]),ylab='Likelihood')
    for(b in 2:20){
      points(seq(5, bMax-5, by = 2),
             exp(out.keep[[b]][,i]-max(out.keep[[b]][,i]))/-sum(out.keep[[b]][,i]),
             typ='l')
    }
    
    plot(seq(5, bMax-5, by = 2),
         exp(out.keep[[1]][,i+1]-max(out.keep[[1]][,i+1]))/-sum(out.keep[[1]][,i+1]),
         typ='l',main=paste('age',Y2$Group.1[i+1]),col='red',ylab='Likelihood')
    for(b in 2:20){
      points(seq(5, bMax-5, by = 2),
             exp(out.keep[[b]][,i+1]-max(out.keep[[b]][,i+1]))/-sum(out.keep[[b]][,i+1]),
             typ='l',col='red')
    }
  }
  
}
