#####
##### Get DOI for each site #####
#####


library(neotoma)

meta_look <- get_dataset(datasettype='pollen', gpid=c("Wisconsin", 
                                                 "Michigan",
                                                 "Minnesota",
                                                 "Illinois",
                                                 "Indiana"), ageyoung=0)
meta_dl <- get_download(test)


test <- get_site('Gass Lake')
