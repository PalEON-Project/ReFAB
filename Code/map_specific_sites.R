load("all.preds.min.list.RData")
all.preds = summary(csamp.real.pred)
#follows from Run_Model.R ##all.preds = summary(csamp.real.pred)

site.factor = factor(x[,1],labels = seq(1,142,1))
all.preds1 = cbind(site.factor,x[,1:6],all.preds$quantiles[,1],all.preds$quantiles[,3],all.preds$quantiles[,5])
all.preds1 = cbind(all.preds1,ten.count)
all.preds1 = all.preds1[order(all.preds1$LatitudeNorth),]
all.preds1[,1] = factor(all.preds1[,2],labels = seq(1,142,1),
                        levels=unique(all.preds1[,2]),ordered=FALSE)

#ten.count = ten.count[order(all.preds1$LatitudeNorth),]

count_mat = count(df = as.data.frame(all.preds1),vars = "SiteID")

add1 = matrix(0,nrow(all.preds1),1)

for(i in 1:nrow(all.preds1)){
  add1[i,] = count_mat[count_mat[,1]==all.preds1$SiteID[i],2]
}

plot_sites = c(1988,269,1593,1820,2933,8574,318,1101)
quartz()
map('state', xlim=range(all.preds1[,4])+c(-2, 2), ylim=range(all.preds1[,3])+c(-1, 1))
for(i in 1:length(plot_sites)){
points(all.preds1[all.preds1[,2]==plot_sites[i],4],
       all.preds1[all.preds1[,2]==plot_sites[i],3], pch=19, cex=1)
text(unique(all.preds1[all.preds1[,2]==plot_sites[i],4]),
     unique(all.preds1[all.preds1[,2]==plot_sites[i],3]+.2),
     labels = unique(all.preds1[all.preds1[,2]==plot_sites[i],2]),cex=.6)
}
title("Sensitivity Analysis Sites")

