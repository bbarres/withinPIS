################################################################################
################################################################################
#withinPIS code
################################################################################
################################################################################

#loading the required packages
library(maptools)
library(rgdal)
library(OpenStreetMap)


################################################################################
#loading GIS layers
################################################################################

#setting the path to the datasets
setwd("~/work/Rfichiers/Githuber/withinPIS_data")

#choice of the patches
patch<-read.table("coord_all_patch9.txt",header=TRUE,sep="\t")
summary(patch$PatchArea)
selecpatch<-read.table("selected_patchesPIS14.txt",sep="\t",header=TRUE)
Aland<-readShapePoly("RECNO.shp",proj4string=CRS("+init=epsg:2393"))
Aland<-spTransform(Aland,CRS("+init=epsg:3067"))
patchshape<-readShapePoly("All patches.shp",
                          proj4string=CRS("+init=epsg:3067"))

plot(Aland,col=grey(0.85),lty=0)
plot(patchshape,col="red",lty=0,add=TRUE)
plot(patchshape[patchshape[[3]] %in% selecpatch[,1],1],col="blue",lty=0,
     add=TRUE)

points(patch[patch$ID %in% selecpatch$patch_ID,c("Longitude","Latitude")],
       cex=1,bg="red",pch=21)
points(patch[patch$ID %in% selecpatch[selecpatch$survey=="Benoit",]$patch_ID,
             c("Longitude","Latitude")],
       cex=1,bg="blue",pch=21) #export pdf 7 x 6 inches

plot(Aland,col=grey(0.85),lty=0)
points(patch[patch$ID %in% selecpatch$patch_ID,c("Longitude","Latitude")],
       cex=3,bg="red",pch=21)
points(patch[patch$ID %in% selecpatch[selecpatch$survey=="Benoit",]$patch_ID,
             c("Longitude","Latitude")],
       cex=3,bg="blue",pch=21) #export png 1300 x 1167 pixels

#Mapping using openStreetMap and raster instead of shapefile
patchshape<-spTransform(patchshape, osm()) #project the coordinates in the 
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
plot(patchshape,col="red",lty=0,add=TRUE)
plot(patchshape[patchshape[[3]] %in% selecpatch[,1],1],col="blue",lty=0,add=TRUE)

#example for the patch 294
raster294<-openmap(c(60.15604096042903,19.704923629760742),
                   c(60.152912102694046,19.713034629821777),
                   type="esri-topo",minNumTiles=16,zoom=17)
plot(raster294)
plot(patchshape[patchshape[[3]] %in% selecpatch[,1],1],col="blue",lty=0,add=TRUE)


################################################################################
#loading survey data
################################################################################

#loading patches informations
patche_info<-read.table("coord_all_patch12.txt", header=TRUE, sep="\t", dec=".")
#we reorganize the coord data file
patche_info<-patche_info[,c(1:12,14:dim(patche_info)[2])]
points(patche_info$Longitude,patche_info$Latitude,cex=1,bg="white",pch=21)

#loading the sample information data
sample_info<-read.table("sample_info4.txt", header=TRUE, sep="\t", dec=".")

#loading genotypes data
geno_hom<-read.table("geno_hom10_13.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
geno_hom<-merge(geno_hom,sample_info,by.x="UNIC_ID",by.y="FIMM_ID")
geno_hom<-merge(geno_hom,patche_info,by.x="patche_ID",by.y="ID",all.x=TRUE)
geno_hom<-replace(geno_hom,geno_hom=="CA","AC")
geno_hom<-replace(geno_hom,geno_hom=="GA","AG")
geno_hom<-replace(geno_hom,geno_hom=="TA","AT")
geno_hom<-replace(geno_hom,geno_hom=="GC","CG")
geno_hom<-replace(geno_hom,geno_hom=="TC","CT")
geno_hom<-replace(geno_hom,geno_hom=="TG","GT")
geno_hom<-droplevels(geno_hom)


################################################################################
#adding the MLG column
################################################################################

#we add a column with the multilocus genotype (MLG)
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