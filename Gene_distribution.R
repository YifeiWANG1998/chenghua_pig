#Gene distribution
#This work refers to https://github.com/TickingClock1992/RIdeogram
import os
import sys
import re


def main(gff_fname,genelist_fname,outfname):
    gene_list = open(genelist_fname).readlines()
    gene_list = [i.strip() for i in gene_list]
    f = open(gff_fname)
    out = open(outfname,'w')
    for line in f:
        for gene in gene_list:
            if gene in line:
                out.write(line)
                continue
    f.close()
    out.close()




if __name__ == '__main__':
    gff_fname = "final.homo_RNA_specific.gff"
    genelist_fname = "gene.txt"
    outfname = "filtered_gff.gff"
    main(gff_fname,genelist_fname,outfname)


library(RIdeogram)

gene_density <- GFFex(input = "filtered_gff.gff.gz",
      karyotype = "pig.fai",
      feature = "mRNA",window = 1000000)


# gene_density <- GFFex(input = "final.homo_RNA_specific.gff.gz",
#                       karyotype = "pig.fai",
#                       feature = "gene",window = 1000000)
# 


pig_karyotype <- read.table("pig.fai", 
                            sep = "\t", header = T, 
                            stringsAsFactors = F)


ideogram(karyotype = pig_karyotype, overlaid = gene_density)
convertSVG("chromosome.svg", device = "png")
