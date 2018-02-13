
#### trees                    the taxa that should be stand alone
#### other trees              are the arboreal taxa that should be grouped
#### cast.x                   the whole dataset from neotoma
#### sites_rm                 the rows that should be removed to maintain the 1/3 gold standard sites
#### all.pollen.taxa.names    the names of all the pollen taxa in the dataset
#### prairie.include          T/F should prairie be included?
#### other.herbs.include      T/F should other herbs be included?
#### drop.taxa                vector of taxa that should be dropped from the final data frame
#### pft.do                   making scenario case where pollen is put into pfts

taxa_selection <- function(trees, other.trees = NULL, cast.x, sites_rm = NULL,
                           all.pollen.taxa.names, prairie.include = TRUE,
                           other.herbs.include = TRUE, other.trees.include=TRUE,
                           drop.taxa = NA, PFT.do = FALSE){
  prairie <- c("ARTEMISIA","ASTERX","POACEAE","AMBROSIA","CHENOAMX","CORYLUS")
  all.veg.box <- cast.x[,all.pollen.taxa.names]
  
  if(!is.na(drop.taxa)){
    all.veg.box <- all.veg.box[,-which(colnames(all.veg.box)%in%drop.taxa)]
  }
  
  tree.box <- all.veg.box[,trees]
  
  if(PFT.do == TRUE){
    ever.box <- all.veg.box[,c("TSUGAX",'PINUSX',"PICEAX","LARIXPSEU","CUPRESSA")]
    decid.box <- all.veg.box[,c("FAGUS","QUERCUS","BETULA","JUGLANSX","ACERX",
                                "FRAXINUX","OSTRYCAR","ULMUS","TILIA","ALNUSX",
                                 "ABIES","POPULUS","CARYA","TAXUS","NYSSA",
                                "CASTANEA","PLATANUS","SALIX","LIQUIDAM")]
    cat.box <- data.frame(prairie = rowSums(all.veg.box[,prairie],na.rm = T),
                          evergreen = rowSums(ever.box, na.rm = T),
                          decidious = rowSums(decid.box, na.rm = T))
    
    test.check <- which(rowSums(cast.x[,all.pollen.taxa.names],na.rm = T) - rowSums(cat.box,na.rm = T)!=0)
    print(test.check)
    print('If any non-zeros returned, this means that not all pollen was included. Is this what you want?')
    
    counts <- round(cat.box)
    counts[is.na(counts)] <- 0
    
    if(any(sites_rm)){
      Y <- counts[-sites_rm,] #remove sites "sites_rm" defined above
      Y <- counts <- Y[,rev(order(colMeans(Y)))]
      print(paste(dim(Y)[1],'<- should be equal to length(biomass)')) # should be zero
    }else{
      Y <- counts <- counts[,rev(order(colMeans(counts)))]
    }
    
    return(Y)
  }

  if(other.trees.include==TRUE){
    other.tree.box <- data.frame(other_trees = rowSums(all.veg.box[,other.trees],na.rm = T))
    cat.box <- cbind(tree.box,other.tree.box)
  }else{
    cat.box <- tree.box
  }
  
  if(prairie.include == TRUE){
    prairie.box <- data.frame(prairie = rowSums(all.veg.box[,prairie],na.rm = T))
    cat.box <- cbind(cat.box,prairie.box)
  }else{
    all.veg.box <- all.veg.box[,-which(colnames(all.veg.box)%in%prairie)]
  }
  
  if(other.herbs.include == TRUE){
    herb.box <- data.frame(other_herbs = rowSums(all.veg.box,
                        na.rm = T) - rowSums(cat.box,na.rm = T))
    cat.box <- cbind(cat.box,herb.box)
  }
  
  test.check <- which(rowSums(cast.x[,all.pollen.taxa.names],na.rm = T) - rowSums(cat.box,na.rm = T)!=0)
  print(test.check)
  print('If any non-zeros returned, this means that not all pollen was included. Is this what you want?')

  counts <- round(cat.box)
  counts[is.na(counts)] <- 0

  if(any(sites_rm)){
    Y <- counts[-sites_rm,] #remove sites "sites_rm" defined above
    Y <- counts <- Y[,rev(order(colMeans(Y)))]
    print(paste(dim(Y)[1],'<- should be equal to length(biomass)')) # should be zero
  }else{
    Y <- counts <- counts[,rev(order(colMeans(counts)))]
  }
  
  return(Y)
  
}