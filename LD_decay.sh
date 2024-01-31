#Calculation of linkage disequilibrium (LD) decay
#This workflow refers to https://github.com/BGI-shenzhen/PopLDdecay
` -*- coding: utf-8 -*-
import os


def writeandrun_bwa_work(dir_path,dirs):
    for dir in dirs:
        path = os.path.join(dir_path,dir)
        di_path = '/es01/paratera/sce3000'
        dir_path = '/es01/paratera/sce3000'
        work_file = di_path + '/LD_decay_nofiltr_pop/' +dir+'/'+'new_LD_modify'+'.sh'
        f1 = open(work_file,'w')
        vcf = '/es01/paratera/sce3000/filter/nosexsort.no_repeat_snps.filter'
        f1.write('#!/bin/bash\n#SBATCH -p G1Part_sce\n#SBATCH -N 1\n#SBATCH -n 1\n#SBATCH -c 1\n#SBATCH--mem=4g\nexport PATH=/es01/paratera/sce3000/software/PopLDdecay-master/bin:$PATH\n')
        f1.write('PopLDdecay -InVCF ' + vcf +  ' -SubPop div.txt -MaxDist 500 -OutStat ' + os.path.join(os.getcwd(),dir,'newstat')
)
        f1.close()
        now_path = os.getcwd()
        os.chdir(di_path + '/LD_decay_nofiltr_pop/' +dir)
        os.system('sbatch ' + work_file)
        os.chdir(now_path)
if __name__ == '__main__':
    dir_path = '/es01/paratera/sce3000'
    dirs= ['Duroc','Eur_wild','Landrace','Largewhite','MGL']
    writeandrun_bwa_work(dir_path,dirs)


#
` -*- coding: utf-8 -*-
import os


def writeandrun_bwa_work(dir_path,dirs):
    for dir in dirs:
        path = os.path.join(dir_path,dir)
        di_path = '/es01/paratera/sce3000'
        dir_path = '/es01/paratera/sce3000'
        work_file = di_path + '/LD_decay_nofiltr_pop/' +dir+'/'+'new_LD_modify'+'.sh'
        f1 = open(work_file,'w')
        vcf = '/es01/paratera/sce3000/filter/nosexsort.no_repeat_snps.filter'
        f1.write('#!/bin/bash\n#SBATCH -p G1Part_sce\n#SBATCH -N 1\n#SBATCH -n 1\n#SBATCH -c 1\n#SBATCH--mem=4g\nexport PATH=/es01/paratera/sce3000/software/PopLDdecay-master/bin:$PATH\n')
        f1.write('PopLDdecay -InVCF ' + vcf +  ' -SubPop div.txt -MaxDist 500 -OutStat ' + os.path.join(os.getcwd(),dir,'newstat'+ '--OutType 4')
)
        f1.close()
        now_path = os.getcwd()
        os.chdir(di_path + '/LD_decay_nofiltr_pop/' +dir)
        os.system('sbatch ' + work_file)
        os.chdir(now_path)
if __name__ == '__main__':
    dir_path = '/es01/paratera/sce3000'
    dirs= ['Duroc','Eur_wild','Landrace','Largewhite','MGL']
    writeandrun_bwa_work(dir_path,dirs)



###############################################################################
###         For more detail, please see following reference
###############################################################################
## Remington, D. L., Thornsberry, J. M., Matsuoka, Y.,
## Wilson, L. M., Whitt, S. R., Doebley, J., ... & Buckler, E. S. (2001). Structure of linkage disequilibrium
## and phenotypic associations in the maize genome. 
## Proceedings of the national academy of sciences, 98(20), 11479-11484. 
## https://doi.org/10.1073/pnas.201394398
###############################################################################


# import TASSEL LD output file
ld <- read.delim(file.choose(),stringsAsFactors = FALSE,header=TRUE, sep = "\t") 
str(ld)
ld$dist <- as.numeric(ld$X.Dist)
ld$rsq <- ld$R.2
str(ld)
file <- ld[,c(3:5)]
str(file)
Cstart <- c(C=0.1)
modelC <- nls(rsq ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * 
                ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( 2*N*(2+C*dist) * (11+C*dist) ) ), 
              data=file, start=Cstart, control=nls.control(maxiter=100))
rho <- summary(modelC)$parameters[1]
newrsq <- ( (10+rho*file$dist) / ( (2+rho*file$dist) * (11+rho*file$dist) ) ) *
  ( 1 + ( (3+rho * file$dist) * (12+12*rho*file$dist + (rho*file$dist)^2) ) / 
      (2*file$N*(2+rho*file$dist) * (11+rho*file$dist) ) )

newfile <- data.frame(file$dist, newrsq)
str(newfile)
#maxld <- max(newfile$newrsq,na.rm=TRUE) #using max LD value from adjusted data
#halfdecay = maxld*0.5
#halfdecaydist <- newfile$file.dist[which.min(abs(newfile$newrsq-halfdecay))]
newfile <- newfile[order(newfile$file.dist),]
str(newfile)



ld <- read.delim(file.choose(),stringsAsFactors = FALSE,header=TRUE, sep = "\t") 
str(ld)
ld$dist <- as.numeric(ld$X.Dist)
ld$rsq <- ld$R.2
str(ld)
file <- ld[,c(3:5)]
str(file)
Cstart <- c(C=0.1)
modelC <- nls(rsq ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * 
                ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( 2*N*(2+C*dist) * (11+C*dist) ) ), 
              data=file, start=Cstart, control=nls.control(maxiter=100))
rho <- summary(modelC)$parameters[1]
newrsq <- ( (10+rho*file$dist) / ( (2+rho*file$dist) * (11+rho*file$dist) ) ) *
  ( 1 + ( (3+rho * file$dist) * (12+12*rho*file$dist + (rho*file$dist)^2) ) / 
      (2*file$N*(2+rho*file$dist) * (11+rho*file$dist) ) )

newfile2 <- data.frame(file$dist, newrsq)
str(newfile2)
#maxld <- max(newfile$newrsq,na.rm=TRUE) #using max LD value from adjusted data
#halfdecay = maxld*0.5
#halfdecaydist <- newfile$file.dist[which.min(abs(newfile$newrsq-halfdecay))]
newfile2 <- newfile2[order(newfile2$file.dist),]


ld <- read.delim(file.choose(),stringsAsFactors = FALSE,header=TRUE, sep = "\t") 
str(ld)
ld$dist <- as.numeric(ld$X.Dist)
ld$rsq <- ld$R.2
str(ld)
file <- ld[,c(3:5)]
str(file)
Cstart <- c(C=0.1)
modelC <- nls(rsq ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * 
                ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( 2*N*(2+C*dist) * (11+C*dist) ) ), 
              data=file, start=Cstart, control=nls.control(maxiter=100))
rho <- summary(modelC)$parameters[1]
newrsq <- ( (10+rho*file$dist) / ( (2+rho*file$dist) * (11+rho*file$dist) ) ) *
  ( 1 + ( (3+rho * file$dist) * (12+12*rho*file$dist + (rho*file$dist)^2) ) / 
      (2*file$N*(2+rho*file$dist) * (11+rho*file$dist) ) )

newfile3 <- data.frame(file$dist, newrsq)
str(newfile3)
newfile3 <- newfile3[order(newfile3$file.dist),]


ld <- read.delim(file.choose(),stringsAsFactors = FALSE,header=TRUE, sep = "\t") 
str(ld)
ld$dist <- as.numeric(ld$X.Dist)
ld$rsq <- ld$R.2
str(ld)
file <- ld[,c(3:5)]
str(file)
Cstart <- c(C=0.1)
modelC <- nls(rsq ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * 
                ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( 2*N*(2+C*dist) * (11+C*dist) ) ), 
              data=file, start=Cstart, control=nls.control(maxiter=100))
rho <- summary(modelC)$parameters[1]
newrsq <- ( (10+rho*file$dist) / ( (2+rho*file$dist) * (11+rho*file$dist) ) ) *
  ( 1 + ( (3+rho * file$dist) * (12+12*rho*file$dist + (rho*file$dist)^2) ) / 
      (2*file$N*(2+rho*file$dist) * (11+rho*file$dist) ) )

newfile4 <- data.frame(file$dist, newrsq)
str(newfile4)
newfile4 <- newfile4[order(newfile4$file.dist),]


ld <- read.delim(file.choose(),stringsAsFactors = FALSE,header=TRUE, sep = "\t") 
str(ld)
ld$dist <- as.numeric(ld$X.Dist)
ld$rsq <- ld$R.2
str(ld)
file <- ld[,c(3:5)]
str(file)
Cstart <- c(C=0.1)
modelC <- nls(rsq ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * 
                ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( 2*N*(2+C*dist) * (11+C*dist) ) ), 
              data=file, start=Cstart, control=nls.control(maxiter=100))
rho <- summary(modelC)$parameters[1]
newrsq <- ( (10+rho*file$dist) / ( (2+rho*file$dist) * (11+rho*file$dist) ) ) *
  ( 1 + ( (3+rho * file$dist) * (12+12*rho*file$dist + (rho*file$dist)^2) ) / 
      (2*file$N*(2+rho*file$dist) * (11+rho*file$dist) ) )

newfile5 <- data.frame(file$dist, newrsq)
str(newfile5)
newfile5 <- newfile5[order(newfile5$file.dist),]

ld <- read.delim(file.choose(),stringsAsFactors = FALSE,header=TRUE, sep = "\t") 
str(ld)
ld$dist <- as.numeric(ld$X.Dist)
ld$rsq <- ld$R.2
str(ld)
file <- ld[,c(3:5)]
str(file)
Cstart <- c(C=0.1)
modelC <- nls(rsq ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * 
                ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( 2*N*(2+C*dist) * (11+C*dist) ) ), 
              data=file, start=Cstart, control=nls.control(maxiter=100))
rho <- summary(modelC)$parameters[1]
newrsq <- ( (10+rho*file$dist) / ( (2+rho*file$dist) * (11+rho*file$dist) ) ) *
  ( 1 + ( (3+rho * file$dist) * (12+12*rho*file$dist + (rho*file$dist)^2) ) / 
      (2*file$N*(2+rho*file$dist) * (11+rho*file$dist) ) )

newfile6 <- data.frame(file$dist, newrsq)
str(newfile6)
newfile6 <- newfile6[order(newfile6$file.dist),]


ld <- read.delim(file.choose(),stringsAsFactors = FALSE,header=TRUE, sep = "\t") 
str(ld)
ld$dist <- as.numeric(ld$X.Dist)
ld$rsq <- ld$R.2
str(ld)
file <- ld[,c(3:5)]
str(file)
Cstart <- c(C=0.1)
modelC <- nls(rsq ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * 
                ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( 2*N*(2+C*dist) * (11+C*dist) ) ), 
              data=file, start=Cstart, control=nls.control(maxiter=100))
rho <- summary(modelC)$parameters[1]
newrsq <- ( (10+rho*file$dist) / ( (2+rho*file$dist) * (11+rho*file$dist) ) ) *
  ( 1 + ( (3+rho * file$dist) * (12+12*rho*file$dist + (rho*file$dist)^2) ) / 
      (2*file$N*(2+rho*file$dist) * (11+rho*file$dist) ) )
newfile7 <- data.frame(file$dist, newrsq)
str(newfile7)
newfile7 <- newfile7[order(newfile7$file.dist),]


ld <- read.delim(file.choose(),stringsAsFactors = FALSE,header=TRUE, sep = "\t") 
str(ld)
ld$dist <- as.numeric(ld$X.Dist)
ld$rsq <- ld$R.2
str(ld)
file <- ld[,c(3:5)]
str(file)
Cstart <- c(C=0.1)
modelC <- nls(rsq ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * 
                ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( 2*N*(2+C*dist) * (11+C*dist) ) ), 
              data=file, start=Cstart, control=nls.control(maxiter=100))
rho <- summary(modelC)$parameters[1]
newrsq <- ( (10+rho*file$dist) / ( (2+rho*file$dist) * (11+rho*file$dist) ) ) *
  ( 1 + ( (3+rho * file$dist) * (12+12*rho*file$dist + (rho*file$dist)^2) ) / 
      (2*file$N*(2+rho*file$dist) * (11+rho*file$dist) ) )
newfile8 <- data.frame(file$dist, newrsq)

str(newfile8)
newfile8 <- newfile8[order(newfile8$file.dist),]

ld <- read.delim(file.choose(),stringsAsFactors = FALSE,header=TRUE, sep = "\t") 
str(ld)
ld$dist <- as.numeric(ld$X.Dist)
ld$rsq <- ld$R.2
str(ld)
file <- ld[,c(3:5)]
str(file)
Cstart <- c(C=0.1)
modelC <- nls(rsq ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * 
                ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( 2*N*(2+C*dist) * (11+C*dist) ) ), 
              data=file, start=Cstart, control=nls.control(maxiter=100))
rho <- summary(modelC)$parameters[1]
newrsq <- ( (10+rho*file$dist) / ( (2+rho*file$dist) * (11+rho*file$dist) ) ) *
  ( 1 + ( (3+rho * file$dist) * (12+12*rho*file$dist + (rho*file$dist)^2) ) / 
      (2*file$N*(2+rho*file$dist) * (11+rho*file$dist) ) )
newfile9 <- data.frame(file$dist, newrsq)

str(newfile9)
newfile9 <- newfile9[order(newfile9$file.dist),]

ld <- read.delim(file.choose(),stringsAsFactors = FALSE,header=TRUE, sep = "\t") 
str(ld)
ld$dist <- as.numeric(ld$X.Dist)
ld$rsq <- ld$R.2
str(ld)
file <- ld[,c(3:5)]
str(file)
Cstart <- c(C=0.1)
modelC <- nls(rsq ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * 
                ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( 2*N*(2+C*dist) * (11+C*dist) ) ), 
              data=file, start=Cstart, control=nls.control(maxiter=100))
rho <- summary(modelC)$parameters[1]
newrsq <- ( (10+rho*file$dist) / ( (2+rho*file$dist) * (11+rho*file$dist) ) ) *
  ( 1 + ( (3+rho * file$dist) * (12+12*rho*file$dist + (rho*file$dist)^2) ) / 
      (2*file$N*(2+rho*file$dist) * (11+rho*file$dist) ) )
newfile10 <- data.frame(file$dist, newrsq)

str(newfile10)
newfile10 <- newfile10[order(newfile10$file.dist),]




ld <- read.delim(file.choose(),stringsAsFactors = FALSE,header=TRUE, sep = "\t") 
str(ld)
ld$dist <- as.numeric(ld$X.Dist)
ld$rsq <- ld$R.2
str(ld)
file <- ld[,c(3:5)]
str(file)
Cstart <- c(C=0.1)
modelC <- nls(rsq ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * 
                ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( 2*N*(2+C*dist) * (11+C*dist) ) ), 
              data=file, start=Cstart, control=nls.control(maxiter=100))
rho <- summary(modelC)$parameters[1]
newrsq <- ( (10+rho*file$dist) / ( (2+rho*file$dist) * (11+rho*file$dist) ) ) *
  ( 1 + ( (3+rho * file$dist) * (12+12*rho*file$dist + (rho*file$dist)^2) ) / 
      (2*file$N*(2+rho*file$dist) * (11+rho*file$dist) ) )
newfile11 <- data.frame(file$dist, newrsq)

str(newfile11)
newfile11 <- newfile11[order(newfile11$file.dist),]


ld <- read.delim(file.choose(),stringsAsFactors = FALSE,header=TRUE, sep = "\t") 
str(ld)
ld$dist <- as.numeric(ld$X.Dist)
ld$rsq <- ld$R.2
str(ld)
file <- ld[,c(3:5)]
str(file)
Cstart <- c(C=0.1)
modelC <- nls(rsq ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * 
                ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( 2*N*(2+C*dist) * (11+C*dist) ) ), 
              data=file, start=Cstart, control=nls.control(maxiter=100))
rho <- summary(modelC)$parameters[1]
newrsq <- ( (10+rho*file$dist) / ( (2+rho*file$dist) * (11+rho*file$dist) ) ) *
  ( 1 + ( (3+rho * file$dist) * (12+12*rho*file$dist + (rho*file$dist)^2) ) / 
      (2*file$N*(2+rho*file$dist) * (11+rho*file$dist) ) )
newfile12 <- data.frame(file$dist, newrsq)

str(newfile12)
newfile12 <- newfile12[order(newfile12$file.dist),]




ld <- read.delim(file.choose(),stringsAsFactors = FALSE,header=TRUE, sep = "\t") 
str(ld)
ld$dist <- as.numeric(ld$X.Dist)
ld$rsq <- ld$R.2
str(ld)
file <- ld[,c(3:5)]
str(file)
Cstart <- c(C=0.1)
modelC <- nls(rsq ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * 
                ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( 2*N*(2+C*dist) * (11+C*dist) ) ), 
              data=file, start=Cstart, control=nls.control(maxiter=100))
rho <- summary(modelC)$parameters[1]
newrsq <- ( (10+rho*file$dist) / ( (2+rho*file$dist) * (11+rho*file$dist) ) ) *
  ( 1 + ( (3+rho * file$dist) * (12+12*rho*file$dist + (rho*file$dist)^2) ) / 
      (2*file$N*(2+rho*file$dist) * (11+rho*file$dist) ) )
newfile13 <- data.frame(file$dist, newrsq)

str(newfile13)
newfile13 <- newfile13[order(newfile13$file.dist),]



ld <- read.delim(file.choose(),stringsAsFactors = FALSE,header=TRUE, sep = "\t") 
str(ld)
ld$dist <- as.numeric(ld$X.Dist)
ld$rsq <- ld$R.2
str(ld)
file <- ld[,c(3:5)]
str(file)
Cstart <- c(C=0.1)
modelC <- nls(rsq ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * 
                ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( 2*N*(2+C*dist) * (11+C*dist) ) ), 
              data=file, start=Cstart, control=nls.control(maxiter=100))
rho <- summary(modelC)$parameters[1]
newrsq <- ( (10+rho*file$dist) / ( (2+rho*file$dist) * (11+rho*file$dist) ) ) *
  ( 1 + ( (3+rho * file$dist) * (12+12*rho*file$dist + (rho*file$dist)^2) ) / 
      (2*file$N*(2+rho*file$dist) * (11+rho*file$dist) ) )
newfile14 <- data.frame(file$dist, newrsq)

str(newfile14)
newfile14 <- newfile14[order(newfile14$file.dist),]


ld <- read.delim(file.choose(),stringsAsFactors = FALSE,header=TRUE, sep = "\t") 
str(ld)
ld$dist <- as.numeric(ld$X.Dist)
ld$rsq <- ld$R.2
str(ld)
file <- ld[,c(3:5)]
str(file)
Cstart <- c(C=0.1)
modelC <- nls(rsq ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * 
                ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( 2*N*(2+C*dist) * (11+C*dist) ) ), 
              data=file, start=Cstart, control=nls.control(maxiter=100))
rho <- summary(modelC)$parameters[1]
newrsq <- ( (10+rho*file$dist) / ( (2+rho*file$dist) * (11+rho*file$dist) ) ) *
  ( 1 + ( (3+rho * file$dist) * (12+12*rho*file$dist + (rho*file$dist)^2) ) / 
      (2*file$N*(2+rho*file$dist) * (11+rho*file$dist) ) )
newfile15 <- data.frame(file$dist, newrsq)

str(newfile15)
newfile15 <- newfile15[order(newfile15$file.dist),]




# plotting the values
pdf("LD_decay-1.pdf", height=5, width = 5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 
plot(newfile$file.dist, newfile$newrsq, ylim = c(0, 1), pch=".", cex=2, xlab="Distance (bp)", ylab=expression(LD ~ (r^2)), col="#ff0000")

lines(newfile2$file.dist, newfile2$newrsq, col="#000000", lwd=2)
lines(newfile3$file.dist, newfile3$newrsq, col="#2424ff", lwd=2)
lines(newfile4$file.dist, newfile4$newrsq, col="#a226f0", lwd=2)
lines(newfile5$file.dist, newfile5$newrsq, col="#00ff00", lwd=2)
lines(newfile6$file.dist, newfile6$newrsq, col="#ffa300", lwd=2)
lines(newfile7$file.dist, newfile7$newrsq, col="#ffff00", lwd=2)
lines(newfile8$file.dist, newfile8$newrsq, col="#a52a2a", lwd=2)
lines(newfile9$file.dist, newfile9$newrsq, col="#00ffbc", lwd=2)
lines(newfile10$file.dist, newfile10$newrsq, col="#00ffff", lwd=2)
lines(newfile11$file.dist, newfile11$newrsq, col="#bdbdbd", lwd=2)
lines(newfile12$file.dist, newfile12$newrsq, col="#83ccfa", lwd=2)
lines(newfile13$file.dist, newfile13$newrsq, col="#ffd700", lwd=2)
lines(newfile14$file.dist, newfile14$newrsq, col="#ee6363", lwd=2)

lines(newfile15$file.dist, newfile15$newrsq, col="#00b0ee", lwd=2)

dev.off()

#while (!is.null(dev.list()))  dev.off()
#if (dev.cur () > 1) dev.off ()

