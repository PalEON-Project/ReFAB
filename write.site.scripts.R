
# load in data for all sites, per Ann's original code
if(!file.exists('allPredData.Rda'))
  source('prep_data.R') 
load('allPredData.Rda')

site.names <- unique(x.meta$site.name)
how.many <- list()

for(i in 1:length(site.names)){ #44:length(site.names)
  locn <- site.names[i]
  site_number = unique(x.meta[x.meta$site.name == locn,1])
  ten_count_use = ten.count[which(x.meta$site.id == site_number), ]
  
  Y = as.matrix(ten_count_use)
  if(length(Y)>21 & nrow(Y) > 15 &
     max(x.meta[x.meta$site.name == locn,'age_bacon'])>8000 & 
     min(x.meta[x.meta$site.name == locn,'age_bacon'])<2000){
    
    how.many[[i]]<- locn
  }
}

names(how.many) <- site.names[1:177]
how.many <- unlist(how.many)


SITE <- names(how.many)[1]

master.test <- readLines('~/ReFAB/genPareto/master.R')
master.test <- gsub("SITE", SITE, master.test)
locnClean <- gsub(' ', '-', SITE)
writeLines(master.test, con=paste0('master.',locnClean,'.R'))

jobsh <- readLines('~/ReFAB/template.job.sh')
jobsh <- gsub('@SITE@',locnClean,jobsh)
writeLines(jobsh, con=paste0(locnClean,'.sh'))
