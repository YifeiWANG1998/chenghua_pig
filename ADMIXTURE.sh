#ADMIXTURE file preparation and analysis workflow
#This workflow refers to https://github.com/chrchang/plink-ng
#This workflow refers to https://github.com/vcftools/vcftools
#This workflow refers to https://github.com/stevenliuyi/admix
#!/bin/bash
#SBATCH -p G1Part_sce
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH--mem=4g
export PATH=/es01/paratera/sce3000/software/plink-ng-master:$PATH
plink --vcf no_xyautosome_sort_popfilter.vcf --indep-pairwise 50 10 0.2 --out A_ld \
        --allow-extra-chr --set-missing-var-ids @:#

 #!/bin/bash
#SBATCH -p G1Part_sce
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH--mem=4g
export PATH=/es01/paratera/sce3000/software/plink-ng-master:$PATH
plink --vcf no_xyautosome_sort_popfilter.vcf --make-bed --extract A_ld.prune.in    \
        --out A_LDfiltered --recode vcf-iid --keep-allele-order --allow-extra-chr --set-missing-var-ids @:# 

#!/bin/bash
#SBATCH -p G1Part_sce
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -c 1
#SBATCH--mem=25g
export PATH=/es01/paratera/sce3000/software/vcftools-master/install/bin:$PATH
vcftools --vcf A_LDfiltered.vcf --plink \
        --out A_plink

#!/bin/bash
#SBATCH -p G1Part_sce
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH--mem=4g
export PATH=/es01/paratera/sce3000/software/plink-ng-master:$PATH
plink --noweb --file A_plink --recode12 --out A_admixture       \
        --allow-extra-chr --keep-allele-order

for k in {2..6};do
        ./admixture -j2 -C 0.01 --cv A_admixture.ped $k >A_admixture.log$k.out
done

