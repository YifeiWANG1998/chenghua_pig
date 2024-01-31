#Files preparation and calculation of runs of homozygosity (ROH)
#This workflow refers to https://github.com/chrchang/plink-ng
#` -*- coding: utf-8 -*-
import os


def writeandrun_bwa_work(dir_path,dirs):
    for dir in dirs:
        path = os.path.join(dir_path,dir)
        di_path = '/es01/paratera/sce3000'
        dir_path = '/es01/paratera/sce3000'
        work_file = di_path + '/ROH_nopop/' +dir+'/'+'bed'+'.sh'
        f1 = open(work_file,'w')
        vcf = di_path + '/ROH_nofilterpop/' +dir+'/'+'filter'+'.vcf'
        f1.write('#!/bin/bash\n#SBATCH -p G1Part_sce\n#SBATCH -N 1\n#SBATCH -n 1\n#SBATCH -c 1\n#SBATCH--mem=4g\nexport PATH=/es01/paratera/s
ce3000/software/plink-ng-master:$PATH\n')
        f1.write('plink --vcf ' + vcf +  ' --make-bed --out new_snp')
        f1.close()
        now_path = os.getcwd()
        os.chdir(di_path + '/ROH_nopop/' +dir)
        os.system('sbatch ' + work_file)
        os.chdir(now_path)
if __name__ == '__main__':
    dir_path = '/es01/paratera/sce3000'
    dirs= ['BMX','WZS','LC','CH','CW','DWZ','NX','HT','JH','LWH','Min','MS','RC','TT','WJ']
    writeandrun_bwa_work(dir_path,dirs)

#!/bin/bash
#SBATCH -p G1Part_sce
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -c 1
#SBATCH--mem=25g
export PATH=/es01/paratera/sce3000/software/plink-ng-master:$PATH
plink \
        --bfile CH_snp \
        --homozyg \
        --homozyg-density 50 \
        --homozyg-gap 100 \
        --homozyg-kb 500 \
        --homozyg-snp 50 \
        --homozyg-window-het 1 \
        --homozyg-window-snp 50 \
        --homozyg-window-threshold 0.05 \
        --out CH_new