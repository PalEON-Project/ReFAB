
sigma.vals <- c(.01,.03,.09,.27,.81)

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

stop()

sigma.vals <- c(.01,.03,.09,.27,.81)

dataID <- data.frame(name = sort(rep(names(how.many),50)), ID = 1:3100,
                     sigma = rep(sigma.vals,620), group = rep(sort(rep(1:10,5)),62))

new.runs <- c('Tower Lake','Penegor Lake','Cub Lake','Wintergreen Lake','Kirchner Marsh','Gass Lake')
dataID.new <- data.frame(name = sort(rep(new.runs,5)), ID = 3101:3130,
                     sigma = rep(sigma.vals,6), group = rep(NA,30))

dataID <- rbind(dataID,dataID.new)

dataID.new1 <- data.frame(name = c(new.runs,rep('Lake Mendota',6)), ID = 3131:3142,
                         sigma = c(rep(2.43,7),sigma.vals), group = rep(NA,12))

dataID <- rbind(dataID,dataID.new1)

new.sig.vals <- seq(.03,.27,.03)
dataID.new2 <- data.frame(name = sort(rep(names(how.many),90)), ID = 3143:8722,#5580
                     sigma = rep(new.sig.vals,620), group = rep(sort(rep(1:10,9)),62))

dataID <- rbind(dataID,dataID.new2)

dataID.new3 <- data.frame(name = sort(names(how.many)), ID = 8723:8784,#5580
                          sigma = rep(.81,62), group = rep(NA,62))
dataID<-rbind(dataID,dataID.new3)

dataID.new4 <- data.frame(name = sort(names(how.many)), ID = 8785:8846,#5580
                          sigma = rep(.12,62), group = rep(NA,62))
dataID<-rbind(dataID,dataID.new4)

write.csv(dataID,file = 'dataID.csv')

# for(i in 1:length(how.many)){
#   site <- names(how.many)[i]
#   
# 
# for(s in seq_along(sigma.vals)){
#   for(g in 1:10){
#     master.test <- readLines('~/ReFAB/master_cross_validation.R')
#     master.test <- gsub("SIGMA", sigma.vals[s], master.test)
#     master.test <- gsub("GROUP", g, master.test)
#     master.test <- gsub("SITE", site, master.test)
#     locnClean <- gsub(' ', '-', site)
#     locnClean <- gsub("'", '-', locnClean)
#     writeLines(master.test, con=paste0('master.',locnClean,g,sigma.vals[s],'.R'))
#     
#     jobsh <- readLines('~/ReFAB/template.job.sh')
#     jobsh <- gsub('@SITE@',paste0(locnClean,g,sigma.vals[s]),jobsh)
#     writeLines(jobsh, con=paste0(locnClean,g,sigma.vals[s],'.sh'))
#     
#     system(paste('qsub',paste0(locnClean,g,sigma.vals[s],'.sh')))
#   }
# 
# }
# 
# }
