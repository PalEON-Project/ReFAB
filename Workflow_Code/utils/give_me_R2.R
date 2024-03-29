

give_me_R2 <- function(preds,actual){
  rss <- sum(( preds - actual ) ^ 2)  ## residual sum of squares
  tss <- sum((actual - mean(actual)) ^ 2)  ## total sum of squares
  rsq <- 1 - rss/tss
  return(rsq)
}


# 1-(sum(true) - predict^2)/(sum(true)-mean(true)^2)

give_me_R2(1:10,1:10)
give_me_R2(c(rep(0,10),1,2),c(rep(0,10),1,2.5))


