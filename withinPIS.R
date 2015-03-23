################################################################################
################################################################################
#withinPIS code
################################################################################
################################################################################

#loading the required packages
library(maptools)
library(rgdal)
library(OpenStreetMap)
library(mapplots)

#setting the path to the datasets
setwd("~/work/Rfichiers/Githuber/withinPIS_data")


################################################################################
#loading GIS layers
################################################################################

#choice of the patches
patch_info<-read.table("coord_all_patch12.txt",header=TRUE,sep="\t")
selecpatch<-read.table("selected_patchesPIS14.txt",sep="\t",header=TRUE)
Aland<-readShapePoly("RECNO.shp",proj4string=CRS("+init=epsg:2393"))
Aland<-spTransform(Aland,CRS("+init=epsg:3067"))
patchshape<-readShapePoly("All patches.shp",
                          proj4string=CRS("+init=epsg:3067"))

plot(Aland,col=grey(0.85),lty=0)
plot(patchshape,col="red",lty=0,add=TRUE)
plot(patchshape[patchshape[[3]] %in% selecpatch[,1],1],col="blue",lty=0,
     add=TRUE)

points(patch_info[patch_info$ID %in% selecpatch$patch_ID,c("Longitude","Latitude")],
       cex=1,bg="red",pch=21)
points(patch_info[patch_info$ID %in% selecpatch[selecpatch$survey=="Benoit",]$patch_ID,
             c("Longitude","Latitude")],
       cex=1,bg="blue",pch=21) #export pdf 7 x 6 inches

plot(Aland,col=grey(0.85),lty=0)
points(patch_info[patch_info$ID %in% selecpatch$patch_ID,c("Longitude","Latitude")],
       cex=3,bg="red",pch=21)
points(patch_info[patch_info$ID %in% selecpatch[selecpatch$survey=="Benoit",]$patch_ID,
             c("Longitude","Latitude")],
       cex=3,bg="blue",pch=21) #export png 1300 x 1167 pixels

#Mapping using openStreetMap and raster instead of shapefile
patchshape.osm<-spTransform(patchshape, osm()) #project the coordinates in the 
                                           #right system
launchMapHelper() #in order to easily aim at the area you want to download
rasterAlandmap<-openmap(c(60.4592491336415,19.45404052734375),
                        c(59.902024360718464,20.88775634765625),
                        type="stamen-toner",minNumTiles=16,zoom=11)
rasterAlandmap<-openmap(c(60.4592491336415,19.45404052734375),
                        c(59.902024360718464,20.88775634765625),
                        type="esri-topo",minNumTiles=16,zoom=11)
rasterAlandmap<-openmap(c(60.4592491336415,19.45404052734375),
                        c(59.902024360718464,20.88775634765625),
                        type="bing",minNumTiles=16,zoom=11)
plot(rasterAlandmap)
plot(patchshape.osm,col="red",lty=0,add=TRUE)
plot(patchshape.osm[patchshape.osm[[3]] %in% selecpatch[,1],1],col="blue",lty=0,
     add=TRUE)

#example for the patch 294
raster294<-openmap(c(60.15604096042903,19.704923629760742),
                   c(60.152912102694046,19.713034629821777),
                   type="esri-topo",minNumTiles=16,zoom=17)
plot(raster294)
plot(patchshape.osm[patchshape.osm[[3]] %in% selecpatch[,1],1],col="blue",lty=0,
     add=TRUE)

#plotting only one patch at a time 
plot(patchshape[patchshape[[3]] %in% selecpatch[1,1],1],col="white",lty=1)
plot(patchshape[patchshape[[3]] %in% selecpatch[2,1],1],col="red",lty=1)
plot(patchshape[patchshape[[3]] %in% selecpatch[3,1],1],col="purple",lty=1)


################################################################################
#loading and polishing the survey data
################################################################################

#we reorganize the coord data file
patch_info<-patch_info[,c(1:12,14:dim(patch_info)[2])]
plot(Aland,col=grey(0.85),lty=0)
points(patch_info$Longitude,patch_info$Latitude,cex=1,bg="white",pch=21)

#loading the sample information data
sample_info<-read.table("sample_info5.txt", header=TRUE, sep="\t", dec=".")

#loading genotypes data
#geno_hom<-read.table("GEN_14.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
geno_hom<-read.table("GEN_corrected.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
geno_hom<-merge(geno_hom,sample_info,by.x="UNIC_ID",by.y="FIMM_ID")
geno_hom<-merge(geno_hom,patch_info,by.x="patche_ID",by.y="ID",all.x=TRUE)
geno_hom<-replace(geno_hom,geno_hom=="CA","AC")
geno_hom<-replace(geno_hom,geno_hom=="GA","AG")
geno_hom<-replace(geno_hom,geno_hom=="TA","AT")
geno_hom<-replace(geno_hom,geno_hom=="GC","CG")
geno_hom<-replace(geno_hom,geno_hom=="TC","CT")
geno_hom<-replace(geno_hom,geno_hom=="TG","GT")
geno_hom<-droplevels(geno_hom)

#adding the MLG column
MLG<-geno_hom
MLG<-replace(MLG,MLG=="AA",1)
MLG<-replace(MLG,MLG=="CC",2)
MLG<-replace(MLG,MLG=="GG",3)
MLG<-replace(MLG,MLG=="TT",4)
MLGmat<-vector()
nb_SNP<- 19 #the number of markers used
name_SNP<-dimnames(geno_hom)[[2]][3:(nb_SNP+3-1)]
for (i in 3:(nb_SNP+3-1)) MLGmat<-paste(MLGmat,MLG[,i],sep="/")
geno_hom<-data.frame(geno_hom,"MLG"=MLGmat,stringsAsFactors=FALSE)

#restrict the dataset to the 15 selected patches
dataPIS<-geno_hom[geno_hom$patche_ID %in% selecpatch[,1],]
dataPIS<-droplevels(dataPIS)
dataPISclean<-dataPIS[dataPIS$nb_missing==0,]


################################################################################
#Do "coinfection genotypes" correspond to a mix between existing MLG?
################################################################################

#first we split potential "pure" genotype and mixed genotypes
singlePIS<-dataPIS[dataPIS$nb_snp_het==0 & dataPIS$nb_missing==0,]
mixedPIS<-dataPIS[dataPIS$nb_missing==0 & dataPIS$nb_snp_het>0,]

#function that combines genotypes, for internal use of listgenomix function. 
#geno1 and geno2 are the two genotypes to combine
mimicmix<-function(geno1,geno2){
  comb<-c()
  for (i in 3:(3+nb_SNP-1)) {
    temp<-paste(substr(geno1[i],1,1),substr(geno2[i],1,1),sep="")
    comb<-data.frame(cbind(comb,temp,stringsAsFactors = FALSE))
  }
  return(comb)
}

#function which creates the potential compound genotypes matrix, for each 
#patch. Input file is a "pure" genotypes dataset
listgenomix<-function(parent) {
  combi<-c()
  listPATCH<-table(parent$patche_ID)
  for (k in 1:length(listPATCH)) {
    parentsub<-parent[parent$patche_ID==dimnames(listPATCH)[[1]][k],]
    parentsub<-parentsub[!duplicated(parentsub$MLG),] 
    comb1<-c()
    if (dim(parentsub)[1]==1) {
      comb1<-cbind(parentsub$patche_ID[1],"single.parent","single.parent",
                   parentsub[1,3:(3+nb_SNP-1)])
    } else {
      temp<-c()
      for (i in 1:(dim(parentsub)[1]-1)) {
        for (j in ((i+1):(dim(parentsub)[1]))) {
          temp<-cbind(parentsub$patche_ID[1],parentsub$MLG[i],parentsub$MLG[j],
                      mimicmix(parentsub[i,],parentsub[j,]))
          comb1<-rbind(comb1,temp)
        }
      }
      comb1<-comb1[,-5]
    }
    dimnames(comb1)[[2]]<-c("patche_ID","parent1_ID","parent2_ID",name_SNP)
    combi<-rbind(combi,comb1)
    combi<-as.data.frame(combi,stringsAsFactors=FALSE)
    combi<-replace(combi,combi=="CA","AC")
    combi<-replace(combi,combi=="GA","AG")
    combi<-replace(combi,combi=="TA","AT")
    combi<-replace(combi,combi=="GC","CG")
    combi<-replace(combi,combi=="TC","CT")
    combi<-replace(combi,combi=="TG","GT")
  }
  combi[,4]<-as.character(combi[,4])
  return(combi)
}

mixpotPIS<-listgenomix(singlePIS)

MLG<-mixpotPIS
MLG<-replace(MLG,MLG=="AA",1)
MLG<-replace(MLG,MLG=="CC",2)
MLG<-replace(MLG,MLG=="GG",3)
MLG<-replace(MLG,MLG=="TT",4)
MLGmat<-vector()
nb_SNP<- 19 #the number of markers used
name_SNP<-dimnames(mixpotPIS)[[2]][4:(nb_SNP+4-1)]
for (i in 4:(nb_SNP+4-1)) MLGmat<-paste(MLGmat,MLG[,i],sep="/")
mixpotPIS<-data.frame(mixpotPIS,"MLG"=MLGmat,stringsAsFactors=FALSE)


#this function compare each observed mixed genotype to a set of potential mixed
#genotype resulting from the combination of the "pure" genotypes (see function 
#listgenomix). Input data are a dataframe of the observed mixed genotype and a 
#dataframe of potential mixed genotype (output of the function listgenomix). 
#The output is the same dataframe of observed mixed genotype with 3 
#supplementary columns
compmix<-function(mix,mixpot){
  genealo<-c()
  for (k in 1:dim(mix)[1]) {
    pach<-mix[k,1]
    mixpotsub<-mixpot[mixpot$patche_ID==pach,]
    if (dim(mixpotsub)[1]==0) {
      proba<-c("no_pot",NA,NA)
    } else {
      temp<-c(which(mixedPIS$MLG[k]==mixpotsub$MLG))
      if (length(temp)<1) {
        proba<-c(0,NA,NA)
      } else {
        proba<-c(1,as.character(mixpotsub$parent1_ID[temp]),
                 as.character(mixpotsub$parent2_ID[temp]))
      }
    }
    genealo<-rbind(genealo,proba)
  }
  dimnames(genealo)[[2]]<-c("knownparent","MLG_parent1","MLG_parent2")
  mix<-cbind(mix,genealo)
  return(mix)
}

identimixPIS<-compmix(mixedPIS,mixpotPIS)

#building the final file for plotting
dataplot<-identimixPIS[,c("patche_ID","UNIC_ID","Earthcape_Plant_Code","long_plant",
                      "lat_plant","MLG","knownparent","MLG_parent1",
                      "MLG_parent2")]
dataplot$MLG_parent1<-as.character(dataplot$MLG_parent1)
dataplot$MLG_parent2<-as.character(dataplot$MLG_parent2)
dataplot[is.na(dataplot$MLG_parent1),"MLG_parent1"]<-dataplot$MLG[is.na(dataplot$MLG_parent1)]
dataplot[is.na(dataplot$MLG_parent2),"MLG_parent2"]<-dataplot$MLG[is.na(dataplot$MLG_parent2)]

temp2<-singlePIS[,c("patche_ID","UNIC_ID","Earthcape_Plant_Code","long_plant",
                    "lat_plant","MLG")]
temp2<-cbind(temp2,"knownparent"=1,"MLG_parent1"=temp2$MLG,
             "MLG_parent2"=temp2$MLG,stringsAsFactors=FALSE)

dataplot<-rbind(dataplot,temp2)


#plot of all the MLG at the same time
fsele<-1047
temp<-dataplot[dataplot$patche_ID==fsele,]
temp2<-as.factor(c(temp$MLG_parent1,temp$MLG_parent2))
nb_distinc_MLG<-length(levels(temp2))
levels(temp2)<-rainbow(nb_distinc_MLG)
temp$MLG_parent1<-temp2[1:(length(temp2)/2)]
temp$MLG_parent2<-temp2[(1+length(temp2)/2):(length(temp2))]

plot(patchshape[patchshape[[3]]==fsele,1],col="white",lty=1)
points(temp$long_plant[temp$knownparent==1],
       temp$lat_plant[temp$knownparent==1],cex=0.9,pch=19,
       col=as.character(temp$MLG_parent1[temp$knownparent==1]))
points(temp$long_plant[temp$knownparent==1],
       temp$lat_plant[temp$knownparent==1],cex=0.4,pch=19,
       col=as.character(temp$MLG_parent2[temp$knownparent==1]))
points(temp$long_plant[temp$knownparent==0],
       temp$lat_plant[temp$knownparent==0],cex=0.9,pch=17,
       col=as.character(temp$MLG_parent1[temp$knownparent==0]))
points(temp$long_plant[temp$knownparent==0],
       temp$lat_plant[temp$knownparent==0],cex=0.4,pch=17,
       col=as.character(temp$MLG_parent2[temp$knownparent==0]))








