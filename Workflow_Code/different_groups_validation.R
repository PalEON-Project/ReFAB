
######
###### Grouping species by traits
######

####
#### Shade Tolerance from LINKAGES Model Parameterization
####
spp.params <- read.csv('~/linkages_package/inst/spp_matrix.csv') #from https://github.com/araiho/linkages_package
shade_list <- list()
for(i in 1:length(trees_look)){
  shade_list[[i]] <- spp.params[grep(trees_look[i],spp.params$Spp_Name),'ITOL']
}
names(shade_list) <- trees_look ## seems low even for grams
sort(unlist(lapply(shade_list,mean,na.rm=T))) ## right order?


####
#### SSD from TRY database
####
SSD <- read.csv(file.path('Data','SSD-database.csv'))

trees_look <- c('Juglans','Fraxinu','Ostry','Carpin','Ulmus','Tilia','Carya','Fagus','Tsuga',
                'Quercus','Betula','Pinus','Acer','Alnus','Picea','Abies',
                'Populus','Larix','Pseu','Cupres','Castanea','Platanus','Salix','Liquidam',
                'Taxus','Nyssa')

WD <- list()
for(i in 1:length(trees_look)){
  spec <- SSD[grep(trees_look[i],SSD$Binomial,ignore.case = T),]
  if(trees_look[i]=='cupres'){
    WD[[i]] <- mean(spec[,4])
  }else{
    WD[[i]] <- spec[spec$Region=='NorthAmerica',4]
  }
}
names(WD) <- trees_look

####
#### Seeds from USDA
####
seeds_dat <- read.csv(file.path('Data','seeds_per_pound_USDA.csv'))

seeds_list <- list()
for(i in 1:length(trees_look)){
  seeds_list[[i]] <- udunits2::ud.convert(seeds_dat[grep(trees_look[i],seeds_dat$Scientific.Name),
                                                    'Seeds.per.Pound'],
                                          '1/lb','g')
}
names(seeds_list) <- trees_look ## seems low even for grams
sort(unlist(lapply(seeds_list,mean,na.rm=T))) ## right order?

####
#### LMA from BETY database
####

SLA <- traits::betydb_search('SLA')
save(SLA,file=file.path('Data','SLA.Rdata'))
load(file.path('Data','SLA.Rdata'))

sla_list <- list()
for(i in 1:length(trees_look)){
  
  sla_list[[i]] <- SLA$stat[grep(trees_look[i],SLA$genus)]
  
}
names(sla_list) <- trees_look
sort(unlist(lapply(sla_list,mean,na.rm=T))) ## right order?


####
#### Determine trait goups
####

ag <- matrix(NA,nrow=length(trees_look),ncol = 3)
ag[,1] <- unlist(lapply(WD,mean,na.rm=T))
ag[,2] <- unlist(lapply(seeds_list,mean,na.rm=T))
ag[,3] <- unlist(lapply(sla_list,mean,na.rm=T))

ag[is.na(ag)] <-0
rownames(ag) <- trees_look
colnames(ag) <- c('WoodDensity','SeedMass','SpecificLeafArea')

x <- prcomp(ag)
biplot(x,cex=1)
names(x)

sla_group <- which(x$x[,2]>.3)
seed_mass_group <- which(x$x[,1]< -2)
wood_density_group <- trees_look[-which((trees_look)%in%names(c(seed_mass_group,sla_group)))]

#######
####### Calibration / Validation for 5 different groupings
#######


load('two.thirds.cast.x.Rdata')
load('twothirds_v2.0.Rdata')
load('nimble_pull2018-10-31.Rdata')
trees <- c("JUGLANSX","FRAXINUX","OSTRYCAR","ULMUS","TILIA","CARYA",
           "FAGUS","TSUGAX","QUERCUS","BETULA",
           'PINUSX',"ACERX","ALNUSX","PICEAX","ABIES","POPULUS",
           "LARIXPSEU","CUPRESSA") #
other.trees <- c("CASTANEA","PLATANUS","LIQUIDAM","TAXUS","NYSSA")#NULL#c()
drop.taxa <- NA#c('other_herbs')
all.pollen.taxa.names <- colnames(pol_cal_count)[11:length(colnames(pol_cal_count))]

source(file.path('Workflow_Code','utils','taxa_selection.R'))
Y.pft <- taxa_selection(trees = trees, other.trees = other.trees,
                    cast.x = ag.two.thirds.cast.x, sites_rm = 0,
                    all.pollen.taxa.names = all.pollen.taxa.names,
                    prairie.include = F, bigwoods.include=F, other.herbs.include = F,
                    other.trees.include = F, drop.taxa = drop.taxa,
                    PFT.do = T)

Y.biome <- taxa_selection(trees = trees, other.trees = other.trees,
                        cast.x = ag.two.thirds.cast.x, sites_rm = 0,
                        all.pollen.taxa.names = all.pollen.taxa.names,
                        prairie.include = F,bigwoods.include=F, other.herbs.include = F,
                        other.trees.include = F, drop.taxa = drop.taxa,
                        biome.do = T)

Y.PFT.NEW <- taxa_selection(trees = trees, other.trees = other.trees,
                       cast.x = ag.two.thirds.cast.x, sites_rm = 0,
                       all.pollen.taxa.names = all.pollen.taxa.names,
                       prairie.include = F,bigwoods.include=F, other.herbs.include = F,
                       other.trees.include = F, drop.taxa = drop.taxa,
                       PFT.NEW.do = T)

Y.succ <- taxa_selection(trees = trees, other.trees = other.trees,
                        cast.x = ag.two.thirds.cast.x, sites_rm = 0,
                        all.pollen.taxa.names = all.pollen.taxa.names,
                        prairie.include = F,bigwoods.include=F, other.herbs.include = F,
                        other.trees.include = F, drop.taxa = drop.taxa,
                        succession.do = T)

Y.spp <- taxa_selection(trees = trees, other.trees = other.trees,
                        cast.x = ag.two.thirds.cast.x, sites_rm = 0,
                        all.pollen.taxa.names = all.pollen.taxa.names,
                        prairie.include = F,bigwoods.include=F, other.herbs.include = F,
                        other.trees.include = F, drop.taxa = drop.taxa,
                        spp.do = T)

plot_prop <- function(Y,biomass){
  
  props <- prop.table(as.matrix(Y),1)
  
  for(i in 1:ncol(props)){
    plot(biomass,props[,i],pch=19,main=colnames(Y)[i])
  }
  
}

pdf('new_groups_validation_scatters.pdf')
par(mfrow=c(5,3),mai=c(.2,.2,.2,.2))

plot_prop(biomass=biomass,Y=Y.pft)
plot_prop(biomass=biomass,Y=Y.PFT.NEW)
plot_prop(biomass=biomass,Y=Y.biome)
plot_prop(biomass=biomass,Y=Y.succ)
plot_prop(biomass=biomass,Y=Y.spp)

dev.off()

#####
##### Run Validation with the 3 different groups
#####
run_calib_valid <- function(Niters, u_middle, bMax, group_rm, Y.calib,
                            Y.pred,biomass.calib,biomass.pred){
  library(nimble)
  library(splines)
  library(maps)
  library(methods)
  
  ciEnvelope <- function(x,ylo,yhi,...){
    polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                        ylo[1])), border = NA,...) 
  }
  
  #### Making sure Z.knots and u are the same between calibration and validation
  u <- c(0,u_middle,bMax) #c(rep(attr(Z.knots,"Boundary.knots")[1],1),attr(Z.knots,"knots"),rep(attr(Z.knots,"Boundary.knots")[2],1))
  
  source("Workflow_Code/utils/bs_nimble.R")
  Z.test <- matrix(NA,length(biomass.calib),5)
  for(i in 1:length(biomass.calib)){
    Z.test[i,] <- bs_nimble(u_given = biomass.calib[i], u = u, N0 = rep(0, (length(u)-1)), N1 = rep(0, (length(u))),
                            N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))
  }
  
  Z.knots <- Z.test
  
  source(file.path('Workflow_Code','models','calibration.model.R'))
  samples.mixed <- calibration_model(Y = Y.calib, biomass = biomass.calib,
                                     Z.knots = Z.knots, u = u, Niters = Niters,
                                     group_rm = group_rm)
  
  burnin <- round(.2 * nrow(samples.mixed))
  new.biomass <- 1:bMax
  
  Z.new = matrix(0,nrow=length(new.biomass),ncol=ncol(Z.knots))
  for(i in 1:length(new.biomass)){
    u_given <- new.biomass[i]
    Z.new[i,] = bs_nimble(u_given, u=u, N0 = rep(0, (length(u)-1)),
                          N1 = rep(0, (length(u))), 
                          N2 = rep(0, (length(u)+1)), 
                          N3 = rep(0, (length(u)+2)))
  }
  source(file.path('Workflow_Code','utils','getLik.R'))
  outLik <- getLik(Z = Z.new, u = u, beta = (samples.mixed[nrow(samples.mixed),]),
                   bMax = bMax, Y = Y.pred, knots = length(u)+2)
  
  source(file.path('Workflow_Code','models','validation.R'))
  samples.pred <- validation_model(Y = Y.pred, Z.knots = Z.knots, 
                                   samples.mixed = samples.mixed, u = u,
                                   Niters = Niters, bMax = bMax, group_rm = group_rm,
                                   outLik = outLik)
  
}

plot_calilb_valid <- function(biomass,samples.pred,bMax,TITLE){
  #par(mfrow=c(1,1))
  plot(biomass, colMeans(samples.pred, na.rm = T),
       xlim=c(0,bMax), ylim=c(0,bMax), pch=19,
       xlab="True Biomass", ylab="Predicted Mean Biomass")
  abline(a=0,b=1)
  lm.mod <- lm(biomass~colMeans(samples.pred)+0)
  abline(lm.mod,lty=2)
  
  #points(biomass[sites_rm],colMeans(samples.pred[,sites_rm], na.rm = T),
  #       col='red',pch=19)
  title(main=TITLE)
  legend('bottomright',paste("R2 =",signif(summary(lm.mod)$r.squared,digits = 3)))
  
  arrows(x0 = biomass, y0 = apply(samples.pred,2,FUN = quantile,.05),
         x1 = biomass, y1 = apply(samples.pred,2,FUN = quantile,.975),
         code = 0, lwd=.1)
  points(biomass, colMeans(samples.pred, na.rm = T),pch=21,col='gray')
  #library(calibrate)
  #textxy(biomass,colMeans(samples.pred),1:length(biomass))
}

### Setup
source(file.path('Workflow_Code','utils','validation_args.R')) #file with constants that should be constant between validation exercises

### PFT OLD
run_calib_valid(Niters = Niters, u_middle = u_middle, bMax = bMax, group_rm = 'pft_old',
                Y.calib = Y.pft, Y.pred = Y.pft, biomass.calib = biomass,
                biomass.pred = biomass)
load("~/ReFAB/samples.pred.grouppft_oldbetaNA.Rdata")
samples.pred_old.pft <- samples.pred
plot_calilb_valid(biomass,samples.pred=samples.pred_old.pft,
                  bMax=bMax)

### PFT NEW
run_calib_valid(Niters = Niters, u_middle = u_middle, bMax = bMax, group_rm = 'pft_new',
                Y.calib = Y.PFT.NEW, Y.pred = Y.PFT.NEW, biomass.calib = biomass,
                biomass.pred = biomass)
load("~/ReFAB/samples.pred.grouppft_newbetaNA.Rdata")
samples.pred_new.pft <- samples.pred
plot_calilb_valid(biomass,samples.pred=samples.pred_new.pft,
                  bMax=bMax)

### Biome
run_calib_valid(Niters = Niters, u_middle = u_middle, bMax = bMax, group_rm = 'biome',
                Y.calib = Y.biome, Y.pred = Y.biome, biomass.calib = biomass,
                biomass.pred = biomass)
load("~/ReFAB/samples.pred.groupbiomebetaNA.Rdata")
samples.pred_biome <- samples.pred
plot_calilb_valid(biomass,samples.pred=samples.pred_biome,
                  bMax=bMax)

### SUCC
run_calib_valid(Niters = Niters, u_middle = u_middle, bMax = bMax, group_rm = 'succ',
                Y.calib = Y.succ, Y.pred = Y.succ, biomass.calib = biomass,
                biomass.pred = biomass)
load("~/ReFAB/samples.pred.groupsuccbetaNA.Rdata")
samples.pred_succ<- samples.pred
plot_calilb_valid(biomass,samples.pred=samples.pred_succ,
                  bMax=bMax)

### SPP
run_calib_valid(Niters = Niters, u_middle = u_middle, bMax = bMax, group_rm = 'spp',
                Y.calib = Y.spp, Y.pred = Y.spp, biomass.calib = biomass,
                biomass.pred = biomass)
load("~/ReFAB/samples.pred.groupsppbetaNA.Rdata")
samples.pred_spp<- samples.pred
plot_calilb_valid(biomass,samples.pred=samples.pred_spp,
                  bMax=bMax)


pdf('new_groups_validation_scatters_r2.pdf')
par(mfrow=c(5,4),mar=c(2,2,1,1), oma = c(1.5, 1.5, 1.5, 0))

plot_prop(biomass=biomass,Y=Y.pft)
plot_calilb_valid(biomass,samples.pred=samples.pred_old.pft,
                  bMax=bMax,TITLE = 'Old PFTs')

plot_prop(biomass=biomass,Y=Y.PFT.NEW)
plot_calilb_valid(biomass,samples.pred=samples.pred_new.pft,
                  bMax=bMax,TITLE = 'New PFTs')

plot_prop(biomass=biomass,Y=Y.biome)
plot_calilb_valid(biomass,samples.pred=samples.pred_biome,
                  bMax=bMax,TITLE = 'Biome Types')

plot_prop(biomass=biomass,Y=Y.succ)
plot_calilb_valid(biomass,samples.pred=samples.pred_succ,
                  bMax=bMax,TITLE = 'Successional Status')

plot_prop(biomass=biomass,Y=Y.spp)
plot_calilb_valid(biomass,samples.pred=samples.pred_spp,
                  bMax=bMax,TITLE = 'Species')

mtext(side = 2, text = 'Pollen Proportion', outer = T, line = 0)
mtext(side = 1, text = 'Settlement Biomass (Mg/ha)', outer = T, line = .5,at=.4)
dev.off()

