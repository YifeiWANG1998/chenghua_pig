#Workflow of neighbor-joining (NJ) tree construction
#This workflow refers to https://github.com/BGI-shenzhen/VCF2Dis
#!/bin/bash
#SBATCH -p G1Part_sce
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -c 1
#SBATCH--mem=25g
export PATH=/es01/paratera/sce3000/software/VCF2Dis-master/bin:$PATH
VCF2Dis -InPut ./nr6a1.vcf.recode.vcf -OutPut nr6a1.vcf.recode.vcf.mat

f1 = open('./nr6a1.vcf.recode.vcf.mat')
res_f = open('./re_name_81nr6a1.vcf.recode.vcf.mat','w')
count = 0 
for line in f1:
    if count == 0:
        res_f.write(line)
    else:
        parts = line.split()
        parts[0] = 'A' + str(count).zfill(3)
        res_f.write("\t".join(parts) + '\n')
    count += 1
res_f.close()