library(ncdf4) #1.17.1
library(mgcv) #v1.8.38
library(maps) #v3.4.0

logged <- FALSE

nc <-
  nc_open('ReFAB_site_reconstruction_v2.0.nc')

load('albers_site_hull.Rdata')

agwb <- ncvar_get(nc, "AGWB")

dm <- dim(agwb)
nAges <- dm[4]
nSites <- dm[1]
nReps <- dm[3]


## nc$dim$lon$vals -> lon
## nc$dim$lat$vals -> lat

x <- albers.df_site$x/100000
y <- albers.df_site$y/100000


YBP <- ncvar_get(nc, "time")


## old code for lat/lon
##    load('pts_in.Rdata') #pts in convexhull
##    load('coors_dat.Rdata') #upper midwest
##    coors_dat_use <- coors_dat[pts_in,] #convexhull
coors_dat_use <- albers.df_hull/100000
    
expanded <- expand.grid(idx = 1:nrow(coors_dat_use), age = 1:nAges)
tmp <- data.frame(idx = 1:nrow(coors_dat_use), coors_dat_use)

fullnewdata <- merge(tmp, expanded)[,-1]
fullnewdata <- fullnewdata[order(fullnewdata[,'age']),]

nDraws <- 25

preds <- array(0, c(nrow(coors_dat_use), nAges, nDraws, nReps))

for(iter in seq_len(nReps)) {

    ## single draw
    agwb_time_slice_mean <- matrix(0, nSites, nAges)
    for(tt in 1:nAges) 
        agwb_time_slice_mean[ , tt] <- diag(agwb[ , , iter, tt])
    
    all.preds <- cbind(rep(x, nAges),
                       rep(y, nAges),
                       sort(rep(YBP, nSites)),
                       c(agwb_time_slice_mean))
        
    # head(all.preds)
    colnames(all.preds) <- c('x','y', 'age', 'agwb')
    
    if(logged)
        all.preds[,'agwb'] <- log(all.preds[,'agwb'])
    all.preds[,'age']  <- all.preds[,'age'] / 100
    
    ## uses ~3 GB memory, 80 secs with ~8 cores
    system.time(
        b <- gam(agwb ~ te(
                     x,
                     y,
                     age,
                     d = c(2, 1),
                     bs = c("tp", "cr"),
                     k = 40
                 ),
                 data = as.data.frame(all.preds))
    )

    ## Quasi-posterior draws following Wood GAM book.
    
    ## uses 36 GB memory
    if(iter == 1) {
        Xp <- predict(b, newdata = fullnewdata, type = 'lpmatrix')
    }
    br <- rmvn(n = nDraws, coef(b), vcov(b))
    
    ## could do this as full matrix-matrix multiply
    for(i in 1:nDraws) {
        preds[ , , i, iter] <- Xp %*% br[i, ]
    }
    print(c(date(), iter))
}

rm(Xp)
rm(br)
fn <- 'preds.Rda'
if(logged) fn <- 'preds_log.Rda'
save.image(file.path('/tmp', fn))

## Fitting on original scale means we get negative predictions.
## Per PLS work, fitting on log scale would probably bias estimates downwards.
preds[preds < 0] <- 0

## Fig 2a, 2c
predsSpaceTimePostMeans <- apply(preds, c(1,2), mean)
## Fig 2b
spatialAvgDrawsByTime <- apply(preds, c(2,3,4), mean)

if(FALSE){
#    preds[preds>350] <- 350
    summed <- apply(preds, c(2,3,4), sum)
    mn <- apply(summed, 1, mean)/6561 # try median too
    low <- apply(summed, 1, quantile, .025)/6561
    hi <- apply(summed, 1, quantile, .975)/6561

    low <- apply(summed, 1, quantile, .1)/6561
    hi <- apply(summed, 1, quantile, .9)/6561

    plot(1:100, rev(mn), type = 'l', ylim = c(45,110), xlab = 'time point', ylab = "avg Mg/ha")
    lines(1:100, rev(low), lty = 2)
    lines(1:100, rev(hi), lty = 2)

    plot(1:100, rev(mn), type = 'l', ylim = c(35, 120))
    tmp <- apply(summed, c(1,3), mean)/6561
    for(i in 101:240) 
        lines(1:100,rev(tmp[,i]),col='red')

    mn <- apply(summed, 1, mean)/6561 # try median too
    low <- apply(summed, 1, quantile, .025)/6561
    hi <- apply(summed, 1, quantile, .975)/6561
    plot(1:100, rev(mn), type = 'l', ylim = c(45,110), xlab = 'time point', ylab = "avg Mg/ha")
    lines(1:100, rev(low), lty = 2)
    lines(1:100, rev(hi), lty = 2)    
}

