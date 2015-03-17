################################################################################
################################################################################
#withinPIS code
################################################################################
################################################################################

#loading the required packages
library(maptools)
library(rgdal)
library(OpenStreetMap)

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
plot(patchshape.osm[patchshape.osm[[3]] %in% selecpatch[,1],1],col="blue",lty=0,add=TRUE)

#example for the patch 294
raster294<-openmap(c(60.15604096042903,19.704923629760742),
                   c(60.152912102694046,19.713034629821777),
                   type="esri-topo",minNumTiles=16,zoom=17)
plot(raster294)
plot(patchshape.osm[patchshape.osm[[3]] %in% selecpatch[,1],1],col="blue",lty=0,add=TRUE)

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
geno_hom<-read.table("GEN_14.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
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


################################################################################
#Do "coinfection genotypes" correspond to a mix between existing MLG?
################################################################################

#first we split potential "mix partner" and mixed genotypes
singlePIS<-dataPIS[dataPIS$nb_snp_het==0 & dataPIS$nb_missing==0,]
mixedPIS<-dataPIS[dataPIS$nb_snp_het==0 & dataPIS$nb_snp_het>0,]

#function that combines genotypes, for internal use of listgenomix function. geno1 and geno2 are the two 
#genotype to combine
mimicmix<-function(geno1,geno2){
  comb<-c()
  for (i in 6:(6+nb_SNP-1)) {
    temp<-paste(substr(geno1[i],1,1),substr(geno2[i],1,1),sep="")
    comb<-data.frame(cbind(comb,temp,stringsAsFactors = FALSE))
  }
  return(comb)
}

#function which creates the electronical compound genotypes matrix, for each patch. Input file is a "parent" dataset
listgenomix<-function(parent) {
  combi<-c()
  listPATCH<-table(parent$PATCH)
  for (k in 1:length(listPATCH)) {
    parentsub<-parent[parent$PATCH==dimnames(listPATCH)[[1]][k],]
    parentsub<-parentsub[!duplicated(parentsub$MLG),] 
    comb1<-c()
    if (dim(parentsub)[1]==1) {
      comb1<-cbind(parentsub$PATCH[1],"single.parent",parentsub[1,6:(6+nb_SNP-1)])
    } else {
      temp<-c()
      for (i in 1:(dim(parentsub)[1]-1)) {
        for (j in ((i+1):(dim(parentsub)[1]))) {
          temp<-cbind(parentsub$PATCH[1],paste("combgeno",i,"_",j,sep=""),
                      mimicmix(parentsub[i,],parentsub[j,]))
          comb1<-rbind(comb1,temp)
        }
      }
      comb1<-comb1[,-4]
    }
    dimnames(comb1)[[2]]<-c("PATCH","GENOMIX",name_SNP)
    combi<-rbind(combi,comb1)
    combi<-as.data.frame(combi,stringsAsFactors=FALSE)
    combi<-replace(combi,combi=="CA","AC")
    combi<-replace(combi,combi=="GA","AG")
    combi<-replace(combi,combi=="TA","AT")
    combi<-replace(combi,combi=="GC","CG")
    combi<-replace(combi,combi=="TC","CT")
    combi<-replace(combi,combi=="TG","GT")
  }
  return(combi)
}

mixpot2010<-listgenomix(parent2010)
mixpot2011<-listgenomix(parent2011)


#this function compare each observed mixed genotype to a set of potential mixed genotype resulting from 
#the combination of (see function listgenomix). Input data are a dataframe of the observed mixed genotype 
#and a dataframe of potential mixed genotype (output of the function listgenomix). The output is the same 
#dataframe of observed mixed genotype with 3 supplementary columns which are the name of potential 
#combination of 
compmix<-function(mix,mixpot){
  probacol<-c()
  for (k in 1:dim(mix)[1]) {
    pach<-mix[k,1]
    mixpotsub<-mixpot[mixpot$PATCH==pach,]
    if (dim(mixpotsub)[1]==0) {
      proba<-c("no_pot",NA,NA)
    } else {
      qqq<-c()
      for (i in 1:dim(mixpotsub)[1]) {
        val<-c()
        for (j in 6:(6+nb_SNP-1)) {
          if (is.na(mix[k,j])) {
            qq<-0.9      #penalty value if no genotype call in the observed mixed individual
          } else { ifelse (mix[k,j]==mixpotsub[i,j-3],qq<-1,qq<-0) #no penalty (=1) if genotype call 
                   #of the observed mixed individual 
                   #is identical to the electonically 
                   #mixed genotype created by combining 
                   #potential "parent" genotype. If 
                   #there is a discrepancy between 
                   #observed and simulated mixed 
                   #genotype, then the penalty score 
                   #is set to 0
          }
          val<-c(val,qq)
        }
        qqq<-rbind(qqq,val)
      }
      proba<-cbind(apply(qqq,1,prod),apply(qqq,1,sum))
      proba<-c(as.character(mixpotsub[(which(proba==max(proba),arr.ind=TRUE)[1]),2]),
               proba[which(proba==max(proba),arr.ind=TRUE)[1],])
    }
    probacol<-rbind(probacol,proba)
  }
  dimnames(probacol)[[2]]<-c("Combine_ID","proba_prod","proba_sum")
  mix<-cbind(mix,probacol)
  return(mix)
}

compmix2010<-compmix(mix2010,mixpot2010)
compmix2011<-compmix(mix2011,mixpot2011)















