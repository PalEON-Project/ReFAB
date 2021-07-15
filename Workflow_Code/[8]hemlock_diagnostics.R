load('info.Rdata')
now_all <- do.call(rbind,info) #info comes from the average.biomass figure

info <- info[-c(22,27,73)]

#gut check #follows from clustering
plot(unlist(lapply(info, function(x)
  x$lon[1]))[order(agb.mat$name)],
  unlist(lapply(info, function(x)
    x$lat[1]))[order(agb.mat$name)],
  col = colors_tri[clusters@cluster],
  pch = 19)

lons_do <- unlist(lapply(info,FUN=function(x) x$lon[1]))

lons_cuts <- cut(lons_do,breaks=pretty(lons_do),labels=F,include.lowest=T)
lons_cols <- fields::tim.colors(length(pretty(lons_do))-1)

pdf('hemlock_ts_nohighlight.pdf',height=12,width=8)
par(mfrow=c(2,1))
plot(
  info[[47]]$age,
  info[[47]]$TSUGAX,
  xlim = c(100, 0),
  typ = 'b',
  ylim = c(0, .37),
  lwd = 3,
  ylab = 'Hemlock Pollen Proportion',
  xlab = 'Age YBP',
  xaxt = 'n',col='white'
)
abline(v = c(50,30),lty=2)
axis(1,
     at = seq(0, 100, 10),
     labels = seq(0, 10, length.out = 11))
for (i in 1:length(info)) {
  if (max(info[[i]]$TSUGAX) < .075)
    next()
  points(
    info[[i]]$age,
    info[[i]]$TSUGAX,
    xlim = c(100, 0),
    typ = 'l',
    col = lons_cols[lons_cuts[i]],
    lwd = 1.5
  )
}

legend('topleft',title = 'Longitude',legend = pretty(lons_do)[-9],col=lons_cols,lty =1,lwd=2)

plot(
  info[[47]]$age,
  info[[47]]$max_est,
  xlim = c(100, 0),
  typ = 'b',
  ylim = c(0, 250),
  lwd = 3,
  ylab = 'Biomass Maximum Likelihood Estimate',
  xlab = 'Age YBP',
  xaxt = 'n',col='white'
)
abline(v = c(50,30),lty=2)
axis(1,
     at = seq(0, 100, 10),
     labels = seq(0, 10, length.out = 11))
for (i in 1:length(info)) {
  if (max(info[[i]]$TSUGAX) < .075)
    next()
  points(
    info[[i]]$age,
    info[[i]]$max_est,
    xlim = c(100, 0),
    typ = 'l',
    col = lons_cols[lons_cuts[i]],
    lwd = 1.5
  )
}

dev.off()


hem_do <- unlist(lapply(info,FUN=function(x) x$TSUGAX))
biomass_do <- unlist(lapply(info,FUN=function(x) x$max_est))
age_do <- unlist(lapply(info,FUN=function(x) x$age))


info_hem <- info[which(unlist(lapply(info, FUN = function(x) max(x$TSUGAX)))>.075)]

hem_do <- unlist(lapply(info_hem,FUN=function(x) x$TSUGAX))
biomass_do <- unlist(lapply(info_hem,FUN=function(x) x$max_est))
age_do <- unlist(lapply(info_hem,FUN=function(x) x$age))
lon_do <- unlist(lapply(info_hem,FUN=function(x) x$lon))


mean(hem_do[which(age_do<45 & age_do>35)])
mean(hem_do[which(age_do<60 & age_do>50)])

mean(biomass_do[which(age_do<45 & age_do>35)])
mean(biomass_do[which(age_do<60 & age_do>50)])

df <- as.data.frame(scale(data.frame(hem_do,lon_do)))
summary(lm(log(biomass_do) ~ df$hem_do + df$lon_do + df$hem_do*df$lon_do))

age_breaks <- seq(0,100,1)
age_cuts <- cut(age_do,age_breaks,labels=F,include.lowest=T)

which.hem <- which(hem_do > .05)
which.hem2 <- which(hem_do > .05 & age_do > 30 & age_do < 50)
plot(hem_do[which.hem],
     biomass_do[which.hem],
     pch = 19,
     col = viridis::viridis(length(age_breaks)-1)[age_cuts[which.hem]])
points(hem_do[which.hem2],
       biomass_do[which.hem2],co='red')





#only central cluster
idx_sites <- order(agb.mat$name)[which(colors_tri[clusters@cluster]=='magenta3')]

#breaks and colors for age bins
breaks <- seq(0,100,10)
colors <- rev(viridis::viridis(length(breaks)-1))

pdf('hemlock_diagnositcs_bysite_east.pdf',width = 15,height=10)
for(i in idx_sites){
  
  #if(sum(info[[i]]$TSUGAX)==0) next()
  par(mfrow=c(2,3),oma=c(2,2,2,2))
    plot(
      info[[i]]$lon,
      info[[i]]$lat,
      pch = 19,
      col = 'black',
      main = info[[i]]$site[1]
    )
  maps::map('state', add = T)
  
  cuts <- cut(info[[i]]$age,breaks,include.lowest=T,labels=F)
  
  
  ###scatter plots
  plot(info[[i]]$TSUGAX,
       info[[i]]$max_est,
       ylab='max lik biomass estimate',
       xlab = 'hemlock pollen proportion',
       pch = 19,
       col = colors[cuts],
       main=info[[i]]$site[1],
       ylim=c(5,235),xlim=c(0,.25))
  
  age_get <- which(info[[i]]$age<50 &info[[i]]$age>30)
  points(info[[i]]$TSUGAX[age_get],
         info[[i]]$max_est[age_get],col = 'red',pch = 0,cex=1.5)
  
  plot(info[[i]]$TSUGAX + info[[i]]$FAGUS,
       info[[i]]$max_est,
       ylab='max lik biomass estimate',
       xlab = 'hem + beech pollen proportion',
       pch = 19,
       col = colors[cuts],
       main=info[[i]]$site[1],
       ylim=c(5,235),xlim=c(0,.25))
  
  points(info[[i]]$TSUGAX[age_get] + info[[i]]$FAGUS[age_get],
         info[[i]]$max_est[age_get],col = 'red',pch = 0,cex=1.5)
  
  plot.new()
  #### timeseries plots
  plot(
    info[[i]]$age,
    info[[i]]$TSUGAX,
    typ = 'b',
    col = 'limegreen',
    xlim = c(100, 0),
    ylim = c(0,.4),
    xlab = 'age',
    ylab = 'hemlock prop'
  )
  par(new = TRUE)
  plot(
    info[[i]]$age,
    info[[i]]$max_est,
    typ = 'b',
    xlim = c(100, 0),
    ylim = c(0,235),
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  axis(side = 4, at = pretty(range(info[[i]]$max_est)))
  mtext("biomass", side = 4, line = 3)
  legend(
    'topleft',
    c('hemlock proportion', 'max lik biomass'),
    lty = 1,
    col = c('limegreen', 'black')
  )
  
  
  plot(
    info[[i]]$age,
    info[[i]]$TSUGAX + info[[i]]$FAGUS,
    typ = 'b',
    col = 'limegreen',
    xlim = c(100, 0),
    ylim = c(0,.4),
    xlab = 'age',
    ylab = 'hemlock + beech prop'
  )
  par(new = TRUE)
  plot(
    info[[i]]$age,
    info[[i]]$max_est,
    typ = 'b',
    xlim = c(100, 0),
    ylim = c(0,235),
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  axis(side = 4, at = pretty(range(info[[i]]$max_est)))
  mtext("biomass", side = 4, line = 3)
  legend(
    'topleft',
    c('hem + beech proportion', 'max lik biomass'),
    lty = 1,
    col = c('limegreen', 'black')
  )
  
}
dev.off()
