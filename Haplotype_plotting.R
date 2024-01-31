#LD analysis and haplotype plotting
#This workflow refers to https://github.com/BGI-shenzhen/LDBlockShow
#This workflow refers to https://github.com/royfrancis/pophelper
#This workflow refers to https://github.com/raivokolde/pheatmap
#This workflow refers to https://github.com/fabianlindfors/reshape
#This workflow refers to https://github.com/andysouth/rworldxtra
#This workflow refers to https://github.com/topics/cluster
#This workflow refers to https://github.com/kassambara/factoextra
#This workflow refers to https://github.com/tidyverse/ggplot2

#!/bin/bash
#SBATCH -p G1Part_sce
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH--mem=4g
unset I_MPI_PMI_LIBRARY 
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0
export PATH=/es01/paratera/sce3000/software/LDBlockShow-master/bin:$PATH
#export PATH=/es01/paratera/sce3000/software/LDBlockShow-master/install/bin:$PATH
#export PATH=/es01/paratera/sce3000/software/plink-ng-master:$PATH
#LDBlockShow –InVCF Big_minipig18chr.vcf.gz -OutPut out -InGWAS gwas.pvalue -Region chr1:288027825:288278178 -OutPng  -SeleVar 2
LDBlockShow -InVCF 25_indiv_Big_minipig18chr.vcf -OutPut 25_zoom_indiv_block -InGFF /es01/paratera/sce3000/ANNOTATION2/111/02.gene/final.homo_RNA_specific.gff -Region chr1:288084133:288084668 -OutPng  -SeleVar 2

library(pophelper)
install.packages('')
install.packages('pheatmap')
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
annotation_col <- read.table("./annotation4snpworld-1.txt",header=T,stringsAsFactors = T)
rownames(annotation_col) = c(1:)
labels_col = rep(c(""), )
ann_colors = list(
  Region = c(Large="red",Mini="blue" , Medium="yellow"),
#Enter color code
  Breed = c()
)
GeneX <- read.table("./SNP.txt", header=F, stringsAsFactors = F)
colnames(GeneX) <- c(1:)
#out <- pheatmap(GeneX, show_colnames=FALSE, show_rownames=FALSE，legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cluster_row = FALSE, cluster_cols = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="")
#out <- pheatmap(GeneX, show_colnames=TRUE, show_rownames=FALSE, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=FALSE, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="")
out <- pheatmap(GeneX, show_colnames=TRUE, cutree_cols=3, show_rownames =FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cluster_cols = TRUE, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="")
out
install.packages("reshape")
head(GeneX)

install.packages("cluster")
install.packages("factoextra")
latlon <- read.table("laton-1.txt",header=T,stringsAsFactors = F, fill=TRUE)
library(cluster)
library(factoextra)
library(ggplot2)
# install.packages("vctrs")
data2 <- latlon[!is.na(latlon$Lat),]
df = scale(data2)
result <- dist(df, method = "euclidean")

result_hc <- hclust(d = result, method = "ward.D2")
data2$type <- cutree(result_hc, k=24)
lat_mean <- tapply(data2[,1],data2$type,mean,na.rm = TRUE)
lon_mean <- tapply(data2[,2],data2$type,mean,na.rm = TRUE)
data2$cluster1 <- NA
data2$cluster2 <- NA
for(i in 1:) {
  for(j in 1:){
    if(data2[i,3] == j ){
      data2[i,4] <- as.numeric(lat_mean[j])
      data2[i,5] <- as.numeric(lon_mean[j])
    } 
  }
}
new_latlon <- data2[,c(4,5)]
#input the number of samples
rownames(new_latlon) <- c(1:)



#out <- pheatmap(GeneX, show_colnames=TRUE, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="")

names <- as.numeric(colnames(GeneX[,out$tree_col[["order"]]]))


dataset <- new_latlon[names,]
new_latlon$id <- rownames(new_latlon)
dataset$sub <- NA
colnames(dataset)[1:3] <- c("Lat", "Lon", "id")
#head (GeneX)
#head (out)


dataset$id <- rownames(dataset)
head(dataset)
#which(dataset$id=70,)
#which(dataset$id==70,)
#which(dataset$id==95,)
#which(dataset$id==114,)
#which(dataset$id==10,)
#which(dataset$id==142,)
#which(dataset$id==153,)
#which(dataset$id==123,)
#which(dataset$id==133,)
#which(dataset$id==75,)
#which(dataset$id==61,)
#which(dataset$id==133,)
#which(dataset$id==75,)
#which(dataset$id==61,)
#which(dataset$id==84,)
#which(dataset$id==89,)
dataset$sub <- NA
#dataset$sub[1] <- 1
#dataset$sub[2] <- 2
#Input cluster grouping
dataset$sub[:] <- 1
dataset$sub[:] <- 2
#dataset$sub[119:131] <- 3

#dataset$sub[132:145] <- 4
#dataset$sub[146:160] <- 5
#dataset$sub[46:121] <- 6
#dataset$sub[237:248] <- 3
colnames(dataset)[4] <- "Sub.population.3"
dataset$value <- 1
reshape <- cast(dataset,Lat+Lon~Sub.population.3) 
reshape2 <- as.data.frame(reshape)
mapPies(reshape2,xlim=c(60,150),ylim=c(22,50),nameX="Lon",nameY="Lat",nameZs=c('1','2','3'),symbolSize=1.2,   zColours=c('blue','red','yellow'),barOrient='vert',oceanCol="#a7c5c7",landCol="#ede9e0",main="29")
out