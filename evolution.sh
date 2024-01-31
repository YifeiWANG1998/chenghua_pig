#Evolution workflow
#This workflow refers to https://github.com/topics/blastp
#This workflow refers to https://github.com/amkozlov/raxml-ng
#This workflow refers to https://github.com/rcedgar/muscle
###step1  
/share/app/blast-2.2.26/bin/blastall -p blastp -m 8 -a 1 -e 1e-5 -F F -z 30 -d all.pep -i all.pep -o all.pep.m8

###step2
perl orthomcl.pl --pv_cutoff 1e-5 --pi_cutoff 0 --pmatch_cutoff 0 --mode 3 --blast_file all.pep.m8 --gg_file all.gg --outdir s2.orthomcl

###step3
perl orthomcl2stat.pl s2.orthomcl/all_orthomcl.out  Outdir 
perl get-fish.pl Outdir/all_orthomcl.out.single-copy.stat Outdir/all_orthomcl.out.all.stat  

###step4
/share/app/muscle-3.8.31/bin/muscle -quiet -in s3.singlecopy_genefamily/1.pep -out s3.singlecopy_genefamily/1.pep.muscle

####step5
perl supergene_from_famlies-category.pl -cds  s3.singlecopy_genefamily/ > Outdir/single-copy.cds.phy 

####step6
/share/app/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -T 4  -p 12345 -m PROTGAMMAGTR -s single-copy.pep.phy -n pep -x 12345 -f a -o CAPTE -# 100

