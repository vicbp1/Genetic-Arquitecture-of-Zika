# Congenital-Zika-Syndrome
This repository include several scripts used on the article "The genetic architecture of the susceptibility to Congenital Zika syndrome in Brazilians".
The Running_ELR_ver5 code is an updated version to run a genome-wide exact logistic regression.
Currently is running for just one covariate. 

To run the code use the following command line

Run using: R < ELRrun.R <raw.plink.filename> <chromosome number> <plink bim file> <output.filename> <total number of MCMC> <Burnin> <covar file> --no-save > screen_output
  
To run this program you need al least six parameters:
-raw plink file       : Raw plink file obtaned with --recode A on plink, include genotype information\n
-chromosome number    : [0,..,22] a number that select a subset of genotypes\n
-plink bim file       : Bim file of the PLINK binary format of the same dataset\n
-output.filename      : String corresponding to output name\n
-total number of MCMC : Total of iterations to run including burnin\n
-Burnin               : Total number of iterations for burnin iterations\n
-cov file             : Covariate file (with header and TAB delimited), if included, the script will automatically adjust the regression\n
-IMPORTANT 1          : YOU NEED TO INSTALL THE ELRM PACKAGE!!\n
-IMPORTANT 2          : Currently, this version can adjust just for one covariate !!\n
-IMPORTANT 3          : The covariate need to be a discrete variable!!\n
