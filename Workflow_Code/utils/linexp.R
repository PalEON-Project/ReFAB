linexp <- nimbleFunction(
  run = function(x = double(2), J = double(0), I = double(0)) {
    returnType(double(2))
    for(j in 1:J){
      for(i in 1:I){
        if(x[j,i] < 0){
          x[j,i] <- exp(x[j,i])
        }else{
          x[j,i] <- x[j,i] + 1
        }
      }
    }
    return(x)
  })