test_site <- function(x.meta.use){
  if(length(unique(x.meta.use$lat)) != 1) {
    print('lats are not the same for this site name')
    stop()
  }
  if(length(unique(x.meta.use$long)) != 1) {
    print('longs are not the same for this site name')
    stop()
  }
  
  print('Site Passed Tests.')
}