#The analysis between the two groups includes two data calculations: (i) FST and (ii)ROD. 
#This workflow refers to https://github.com/vcftools/vcftools.
#This workflow refers to https://cran.r-project.org/package=argparse/vignettes/argparse.html
#This workflow refers to https://github.com/hardingnj/xpclr
#!/bin/bash
#SBATCH -p G1Part_sce
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH--mem=4g
export PATH=/es01/paratera/sce3000/software/vcftools-master/install/bin:$PATH
vcftools --vcf /es01/paratera/sce3000/population_filter/AAAAno_xyautosome_sort_popfilter.vcf --fst-window-size 20000 --window-pi-step 10000     \
        --weir-fst-pop /es01/paratera/sce3000/population_filter/Selection/big1.txt --weir-fst-pop /es01/paratera/sce3000/population_filter/Selection/mini.txt --out AFst.big.mini


#!/bin/bash
#SBATCH -p G1Part_sce
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH--mem=4g
export PATH=/es01/paratera/sce3000/software/vcftools-master/install/bin:$PATH
vcftools --vcf /es01/paratera/sce3000/population_filter/AAAAno_xyautosome_sort_popfilter.vcf   \
        --window-pi 20000 --window-pi-step 10000        \
        --keep /es01/paratera/sce3000/population_filter/Selection/mini.txt --out Api.mini




library("argparse")
parser <- ArgumentParser(description='The reduction of diversity (ROD), defined as: ROD = 1 - π(domesticated)/π(wild)')

parser$add_argument( "-a", "--domesticated", type="character",required=T,
                help="input popA pi result  file for domesticated pop[required]",
                metavar="domesticated")

parser$add_argument( "-b", "--wild", type="character",required=T,
                help="input popB pi result   file for wild pop[required]",
                metavar="wild")

parser$add_argument( "-o", "--outdir", type="character", default=getwd(),
                help="output file directory [default %(default)s]",
                metavar="path")
parser$add_argument("-p", "--prefix", type="character", default="ROD",
                help="out file name prefix [default %(default)s]",
                metavar="prefix")
opt <- parser$parse_args()

if( !file.exists(opt$outdir) ){
        if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
                stop(paste("dir.create failed: outdir=",opt$outdir,sep=""))
        }
}

pop1<-read.table(opt$domesticated,sep="\t",header=T,comment.char="")
pop1<-data.frame(CHROM=pop1$CHROM,BIN_START=pop1$BIN_START, BIN_END=pop1$BIN_END,PIA=pop1$PI)

#pop1$ID=paste(pop1$CHROM,pop1$BIN_START, pop1$BIN_END)
pop2<-read.table(opt$wild,sep="\t",header=T,comment.char="")
pop2<-data.frame(CHROM=pop2$CHROM,BIN_START=pop2$BIN_START, BIN_END=pop2$BIN_END,PIB=pop2$PI)

data<-merge(pop1,pop2,by=c( "CHROM" ,"BIN_START",  "BIN_END"))
#data<-subset(data,data$PIA/data$PIB>1)

data$ROD<-1-(data$PIA/data$PIB)


colnames(data)<-c( "CHROM" ,"BIN_START",  "BIN_END",opt$domesticated,opt$wild,"ROD")
head(data)
write.table(data,
                file=paste0(opt$outdir,"/",opt$prefix,".txt"), 
                sep = "\t",row.names = F,quote =F)


for i in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18;do
        xpclr --format vcf --input xpclr.vcf       \
                --samplesA phenotype1.txt --samplesB phenotype2.txt    \
                --out ./$i.xpclr --chr $i       \
                --size 200000 --step 100000 --maxsnps 200 --minsnps 10
done