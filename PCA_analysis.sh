#PCA analysis
#This workflow refers to https://github.com/chrchang/plink-ng
#This workflow refers to https://github.com/jianyangqt/gcta
#This workflow refers to https://github.com/Wanyi-Huang/PCA2normal_format.R


#!/bin/bash
#SBATCH -p G1Part_sce
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH--mem=4g
export PATH=/es01/paratera/sce3000/software/plink-ng-master:$PATH
plink --vcf with_id_272_0.9final_version0.9_nosex.vcf --make-bed --out World_snp --chr-set 18 no-xy

#!/bin/bash
#SBATCH -p G1Part_sce 
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH -c 1
#SBATCH--mem=4g
./gcta64 --make-grm --out World_snp.gcta --bfile World_snp --autosome-num 18

#!/bin/bash
#SBATCH -p G1Part_sce 
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH -c 1
#SBATCH--mem=4g
./gcta64 --grm World_snp.gcta --pca 20 --out World_snp.gcta

eigvec <- read.table("World_snp.gcta.eigenvec", header = F, stringsAsFactors = F)
write.table(eigvec[2:ncol(eigvec)], file = "finalWorld_gcta.eigenvector.xls", sep = "\t", row.names = F, col.names = T, quote = F)

eigval <- read.table("World_snp.gcta.eigenval", header = F)
pcs <- paste0("PC", 1:nrow(eigval))
eigval[nrow(eigval),1] <- 0
percentage <- eigval$V1/sum(eigval$V1)*100
eigval_df <- as.data.frame(cbind(pcs, eigval[,1], percentage), stringsAsFactors = F)
names(eigval_df) <- c("PCs", "variance", "proportion")
eigval_df$variance <- as.numeric(eigval_df$variance)
eigval_df$proportion <- as.numeric(eigval_df$proportion)
write.table(eigval_df, file = "finalworld_gcta.eigenvalue.xls", sep = "\t", quote = F, row.names = F, col.names = T)