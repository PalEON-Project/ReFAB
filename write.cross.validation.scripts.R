
sigma.vals <- c(.01,.03,.09,.27,.81)


for(s in seq_along(sigma.vals)){
  for(g in 1:10){
    master.test <- readLines('~/ReFAB/master_cross_validation.R')
    master.test <- gsub("SIGMA", sigma.vals[s], master.test)
    master.test <- gsub("GROUP", g, master.test)
    locnClean <- gsub(' ', '-', SITE)
    locnClean <- gsub("'", '-', locnClean)
    writeLines(master.test, con=paste0('master.',locnClean,g,sigma.vals[s],'.R'))
    
    jobsh <- readLines('~/ReFAB/template.job.sh')
    jobsh <- gsub('@SITE@',paste0(locnClean,g,sigma.vals[s]),jobsh)
    writeLines(jobsh, con=paste0(locnClean,g,sigma.vals[s],'.sh'))
    
    system(paste('qsub',paste0(locnClean,g,sigma.vals[s],'.sh')))
  }

}


