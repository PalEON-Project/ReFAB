#####
##### Get DOI for each site #####
#####


library(neotoma)

meta_look <- get_dataset(datasettype='pollen', gpid=c("Wisconsin", 
                                                 "Michigan",
                                                 "Minnesota",
                                                 "Illinois",
                                                 "Indiana"), ageyoung=0)

test <- get_site('Gass Lake')
meta_dl <- get_download(test)

meta_dl$`860`$dataset

