#Treemix file preparation and analysis process
#This workflow refers to https://github.com/chrchang/plink-ng
#This workflow refers to https://github.com/carolindahms/TreeMix
#!/bin/bash
#SBATCH -p G1Part_sce
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH--mem=4g
export PATH=/es01/paratera/sce3000/software/plink-ng-master:$PATH
plink --vcf Treemix.vcf --indep-pairwise 50 10 0.2 --out tep_ld        \
        --allow-extra-chr --set-missing-var-ids @:# --keep-allele-order

#!/bin/bash
#SBATCH -p G1Part_sce
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH--mem=4g
export PATH=/es01/paratera/sce3000/software/plink-ng-master:$PATH
plink --vcf Treemix.vcf --extract tep_ld.prune.in --freq --missing --within pop.cluster1.txt --out input --allow-extra-chr --set-missing-var-ids @:# --keep-allele-orde

gzip input.frq.strat

for i in {0..11}
do
        treemix -i input_treemix.frq.gz -m $i -o treemix.$i -bootstrap -k 500 > treemix_${i}_log &
done


source("plotting_funcs.R")
getwd()
poporder="pop.order.txt"
outstem="treemix.0"
plot_tree(outstem)
plot_resid(outstem,poporder)