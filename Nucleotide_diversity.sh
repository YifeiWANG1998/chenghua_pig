#Nucleotide diversity calculation
#This workflow refers to https://github.com/vcftools/vcftools
#` -*- coding: utf-8 -*-
import os


def writeandrun_bwa_work(dir_path,dirs):
    for dir in dirs:
        path = os.path.join(dir_path,dir)
        di_path = '/es01/paratera/sce3000'
        dir_path = '/es01/paratera/sce3000'
        work_file = di_path + '/pai_nopo_filter/' +dir+'/'+'pai'+'.sh'
        f1 = open(work_file,'w')
        vcf = '/es01/paratera/sce3000/filter/noxysort.no_repeat_snps.filter.vcf'
        f1.write('#!/bin/bash\n#SBATCH -p G1Part_sce\n#SBATCH -N 1\n#SBATCH -n 1\n#SBATCH -c 1\n#SBATCH--mem=4g\nexport PATH=/es01/paratera/s
ce3000/software/vcftools-master/install/bin:$PATH\n')
        f1.write('vcftools --vcf ' + vcf +  ' --keep div.txt ' + '--window-pi 100000' + ' --out ' + os.path.join(os.getcwd(),dir,'nopopfilter' + '.pi')
)
        f1.close()
        now_path = os.getcwd()
        os.chdir(di_path + '/pai_nopo_filter/' +dir)
        os.system('sbatch ' + work_file)
        os.chdir(now_path)
if __name__ == '__main__':
    dir_path = '/es01/paratera/sce3000'
    dirs= ['Berkshire','BMX','CH','CW','Duroc','DWZ','HT','JH','Landrace','Largewhite','LC','LWH','MGL','Min','MS','NX','RC','TT','WJ','WZS']
    writeandrun_bwa_work(dir_path,dirs)