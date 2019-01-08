test_site <- function(x.meta.use){
  if(length(unique(x.meta.use$lat)) != 1) {
    print(paste('lats are not the same for this site name',unique(x.meta.use$site.name)))
    stop()
  }
  if(length(unique(x.meta.use$long)) != 1) {
    print(paste('longs are not the same for ',unique(x.meta.use$site.name)))
    stop()
  }
  
  print(paste(unique(x.meta.use$site.name),' Passed Tests.'))
}