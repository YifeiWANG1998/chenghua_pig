#GWAS file preparation and analysis workflow
#This workflow refers to https://github.com/chrchang/plink-ng
#This workflow refers to https://github.com/xiaolei-lab/rMVP
#This workflow refers to https://github.com/kaneplusplus/bigmemory
#!/bin/bash
#SBATCH -p G1Part_sce
#SBATCH -N 1
#SBATCH -n 3
#SBATCH -c 1
#SBATCH--mem=12g

/es01/paratera/sce3000/software/plink-ng-master/plink  --vcf bodylength.vcf --make-bed

new("big.matrix.descriptor", description = list(sharedType = "FileBacked", 
    filename = "mvp.plink.geno.bin", dirname = "./", totalRows = 25646934, 
    totalCols = 203L, rowOffset = c(0, 25646934), colOffset = c(0, 
    203), nrow = 25646934, ncol = 203, rowNames = NULL, colNames = NULL, 
    type = "char", separated = FALSE))




library("rMVP")
library('bigmemory')




MVP.Data(fileBed="plink",
         filePhe=NULL,
         fileKin=FALSE,
         filePC=FALSE,       
         #priority="speed",
         #maxLine=10000,
         out="mvp.plink"
)



phePath <- 'myphenotype.txt'
phenotype <- read.table(phePath, header=TRUE)
#print(dim(phenotype))
genoPath <- 'mvp.plink.geno.desc'
genotype <- attach.big.matrix(genoPath)
#print(dim(genotype))
mapPath <- 'plink.bim'
map <- read.table(mapPath , head = FALSE)
map <- map[,c(2,1,4,5,6)]
colnames(map) <- c('SNP','CHROM','POS','REF','ALT')

res <- MVP(phenotype, genotype, map)