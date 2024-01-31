#genome assembly
#This workflow refers to https://github.com/chhylp123/hifiasm
#This workflow refers to https://github.com/lh3/bwa
#This workflow refers to https://www.htslib.org/
# contig assembly
hifiasm  -t 60 --primary  -o test m64154_211103_151527.clean.fasta.gz m64236_211031_214554.clean.fasta.gz m64236_211110_130835.clean.fasta.gz m64236_211112_001301.clean.fasta.gz
awk '/^S/{print ">"$2;print $3}' test.p_ctg.gfa > test.p_ctg.fa

## HiC anchoring
#### step1  
python generate_site_positions.py enzyme test.p_ctg.fa test.p_ctg.fa
bwa index test.p_ctg.fa
/share/app/bwa-0.7.12/bwa mem -t 15 test.p_ctg.fa ganHC_1.fq.gz | samtools sort -@ 10 -n -O bam |python changReadName.py /1 ./ganHC_1.fq.gz.R1.sort.bam
/share/app/bwa-0.7.12/bwa mem -t 15 test.p_ctg.fa ganHC_2.fq.gz | samtools sort -@ 10 -n -O bam |python changReadName.py /2 ./ganHC_1.fq.gz.R2.sort.bam

#### step2 
samtools merge -@ 10 -n -O sam ganHC_1.fq.gz.sam ganHC_1.fq.gz.R1.sort.bam ganHC_1.fq.gz.R2.sort.bam

#### step3
awk -f dups.awk -v name=aligned/ merged_sort.txt
#### step4 运行3d-dna
run-asm-pipeline.sh -m haploid -r 2 test.p_ctg.fa merged_nodups.txt