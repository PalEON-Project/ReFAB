
library(plyr)

hk_counts = read.csv(paste0(data.dir,"hotchkiss_lynch_calcote_counts_v0.csv"))
hk_meta = read.csv(paste0(data.dir,"hotchkiss_lynch_calcote_meta_v0.csv"))

hk_counts[is.na(hk_counts)] <- 0

hk_counts1 = rename(hk_counts,c("Apiaceae.Umbell"="APIACEAE",
                   "Arceuthobium"="ARCEUTHOBI",
                   "Cyperaceae"="CYPERACE",
                   "Equisetum"="EQUISETU",
                   "Ericaceae"="ERICACEX",
                   "Euphorbia"="EUPHORB",
                   "Larix"="LARIXPSEU",
                   "Ostrya"="OSTRYCAR",
                   "Picea"="PICEAX",
                   "Polypodium"="POLYPOD",
                   "Ranunculus"="RANUNCUL",
                   "Rosaceae"="ROSACEAX",
                   "Rumex"="RUMEOXYR",
                   "Selaginella.rupestris"="SELAGINE",
                   "Tsuga"="TSUGAX",
                   "Urtica"="URTICACX"))  

ACERX = rowSums(hk_counts[,c("Acer.rubrum","Acer.saccharum","Acer.spicatum","Acer.undif.")]) #ACERX
ALNUSX = rowSums(hk_counts[,c("alnus.3.pore","Alnus.4","Alnus.5","AlnusUndif" )])#ALNUSX
BETULA = rowSums(hk_counts[,c("Betula","Betulaceae.undif")]) #"BETULA"     
FRAXINUX=rowSums(hk_counts[,c("Fraxinus3","Fraxinus4","FraxUnid")]) #"FRAXINUX"
IVA=rowSums(hk_counts[,c("Iva.ciliata","Iva.xanthifolia")]) #"IVA"   
JUGLANSX=rowSums(hk_counts[,c("Juglans.cinerea","Juglans.nigra","Juglans.undiff")])#"JUGLANSX"                             
LYCOPODX=rowSums(hk_counts[,c("Lyco.annotinum","Lyco.clavatum","Lyco.complanatum","Lyco.lucidulum","Lyco.obscurum","Lyco.selago","Lyco.undiff")])#"LYCOPODX" 
PINUSX=rowSums(hk_counts[,c("Pinus.subg..Pinus","Pinus.subg..Strobus","Pinus.undiff")])#"PINUSX"
POLYGONAX=rowSums(hk_counts[,c("Polygonum.amphibium.type","Polygonum.lapth..Type","Polygonum.sect..Persicaria")])#"POLYGONAX"    
 
cols_keep = which(toupper(colnames(hk_counts1))%in%colnames(pol.cal.count))
hk_counts2 = cbind(hk_counts[,1:7],hk_counts1[,cols_keep],ACERX,ALNUSX,
                   FRAXINUX,IVA,JUGLANSX,LYCOPODX,PINUSX,POLYGONAX)
hk_counts2[,"Betula"] = BETULA
colnames(hk_counts2) <- toupper(colnames(hk_counts2))
Other = rowSums(hk_counts[,8:ncol(hk_counts)]) - rowSums(hk_counts2[,8:ncol(hk_counts2)])

hk_counts3 = cbind(matrix(0,nrow(hk_counts2),10),hk_counts2,Other)
colnames(hk_counts3)

hk_counts3[hk_counts3$NAME=="Deep Pine",]$NAME <- c("Deep Pine")

hk_counts3$NAME = as.factor(as.character(hk_counts3$NAME))

for(i in 1:nrow(hk_counts3)){
  hk_counts3[i,1:10]<-hk_meta[which(hk_meta$name==hk_counts3$NAME[i]),]
}

colnames(hk_counts3)<-c(colnames(hk_meta),colnames(hk_counts3[,11:ncol(hk_counts3)]))

hk_counts3 = as.data.frame(hk_counts3)
save(hk_counts3,file="hk_counts3.csv")
