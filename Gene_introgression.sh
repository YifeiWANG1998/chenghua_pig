
#Gene introgression file preparation and analysis process
#This workflow refers to https://github.com/millanek/Dsuite
#This workflow refers to https://github.com/tidyverse/ggplot2
#This workflow refers to https://github.com/tidyverse

#!/bin/bash
#SBATCH -p G1Part_sce
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -c 1
#SBATCH--mem=8g
export PATH=$PATH:/es01/paratera/sce3000/software/Dsuite/Build
Dsuite Dtrios ./D_ststics-final.vcf sets.txt -t treemix.0.treeout -o sample

#!/bin/bash
#SBATCH -p G1Part_sce
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -c 1
#SBATCH--mem=8g
export PATH=$PATH:/es01/paratera/sce3000/software/Dsuite/Build
Dsuite Dinvestigate ./D_ststics-final.vcf sets.txt test_trios.txt -w 10,10

!/bin/bash
#SBATCH -p G1Part_sce
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -c 1
#SBATCH--mem=8g
export PATH=$PATH:/es01/paratera/sce3000/software/Dsuite/Build
Dsuite Fbranch treemix.0.treeout sample_tree.txt > fbranch.out

library(ggplot2)
library(tidyverse)



my_dat <- tibble::tibble(read.table('./DWZ.txt'))
my_dat <- my_dat |> dplyr::mutate(statu=ifelse(V4>0,'A','B'))
# my_dat <- my_dat |> dplyr::mutate(V4=abs(V4))
my_dat <- my_dat |> dplyr::filter(V2 >=& V2<)
p<-ggplot(my_dat,aes(V2,V4,color=statu,fill=statu))+geom_col()+
  scale_y_continuous(breaks = seq(-1, 1, 0.5), labels = as.character(abs(seq(-1, 1, 0.5))) , limits = c(-1, 1))+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'))+
  scale_color_manual(values=c("B"="#1406eb","A"="#f61005"))+
  scale_fill_manual(values=c("B"="#1406eb","A"="#f61005"))
p
ggsave("DWZ.pdf",p,width = 12,height = 8)