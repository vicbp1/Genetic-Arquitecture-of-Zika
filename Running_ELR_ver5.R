#==========================================
# (C) Copyright 2020, by LNCC and Contributors.
#
# -----------------
#  Script to run EXACT LOGISTIC REGRESSION by Chr ELRrun
#  -----------------
#
# Author: Victor Borda
# Contributions: Meddly Santolalla, Cinthya Cardoso
## 
## This script will win an exact logistic regression with an output including description of the variant, p-value, OR and OR ci.
## !!!! IMPORTANT: MUST HAVE PACKAGE "elrm" INSTALLED !!!!
## For R versions > 3.5 you have to install "elrm" via devtools using the following line:
##
##  require(devtools)
##  devtools::install_version("elrm", "1.2.2")
##
##  To run ELRrun type
##
##  R < ELRrun.R <raw.plink.filename> <chromosome number> <plink bim file> <output.filename> <total number of MCMC> <Burnin> <covar file> <covar number> --no-save
##
## !!!! IMPORTANT 2: Currently, this version can adjust just for one covariate !!!!
## !!!! IMPORTANT 3: Covariate format (3 columns with header and TAB separated):
## !!!!	IID	Phenotype	Covariate

usage=function()
{
	print(noquote("Run using: R < ELRrun.R <raw.plink.filename> <chromosome number> <plink bim file> <output.filename> <total number of MCMC> <Burnin> <covar 				file> --no-save > screen_output"))
	print(noquote("To run this program you need al least six parameters:"))
	print(noquote("raw plink file       : Raw plink file obtaned with --recode A on plink, include genotype information"))
	print(noquote("chromosome number    : [0,..,22] a number that select a subset of genotypes"))
	print(noquote("plink bim file       : Bim file of the PLINK binary format of the same dataset" ))
	print(noquote("output.filename      : String corresponding to output name"))
	print(noquote("total number of MCMC : Total of iterations to run including burnin"))
	print(noquote("Burnin               : Total number of iterations for burnin iterations"))
	print(noquote("cov file             : Covariate file (with header and TAB delimited), if included, the script will automatically adjust the regression"))
	print(noquote("IMPORTANT 1          : YOU NEED TO INSTALL THE ELRM PACKAGE!!!!"))
	print(noquote("IMPORTANT 2          : Currently, this version can adjust just for one covariate !!!!"))
	print(noquote("IMPORTANT 3          : The covariate need to be a discrete variable!!!"))

}

temp=commandArgs()

raw.infile=as.character(temp[2])
if (raw.infile=="help"){usage();q(save='no')}

chr.number=as.numeric(temp[3])
name.map=as.character(temp[4])
output.filename=as.character(temp[5])
total.iterations=as.numeric(temp[6])
burnin.its=as.numeric(temp[7])
cov.filename=as.character(temp[8])

library(plyr)
library(elrm)

### Reading files

exome.raw <- read.table(raw.infile,header=T)
snps.info <- read.table(name.map,header=FALSE,na.strings=c(""),colClasses = c("numeric","character","numeric","numeric","character","character"))

### Subsetting information, selecting data corresponding to the chromosome of command line

column.2.extract<-c()
column.2.extract.bim<-c()

for (allele in c(1:nrow(snps.info))){
	m<-allele+6
	bim.line<-as.numeric(snps.info[allele,1])
	if(isTRUE(chr.number == bim.line)){
	column.2.extract<-c(column.2.extract,m)
	column.2.extract.bim<-c(column.2.extract.bim,allele)
	}
}

snps.info.subset<-snps.info[column.2.extract.bim,]

if (isTRUE(cov.filename=="--no-save")){
	message("######		Running chain Logistic regression without adjusment		######")
	covar_exomes<-exome.raw[,c(2,6,column.2.extract)]
	snps_columns <- 3:dim(covar_exomes)[2]
	alpha=0.05

	## Creating table
	col.names.nonad<-c("TestedAllele","p-value","OR",paste("lower", 100 - 100 *alpha, "ci", sep = ""), paste("upper", 100 - 100 *alpha, "ci", sep = ""))
	p_snps_exact_nonadjust <- read.table(text = "",col.names = col.names.nonad)
	for (i in snps_columns){ 
	  
	  cross.info<-xtabs(~covar_exomes$PHENOTYPE+ interaction(covar_exomes[,i]))
	  number.alleles<-ncol(cross.info)
	  
	  if(number.alleles>2){
		cdat <- cdat <- data.frame(variant = c(0:2),
		Phenotype = cross.info[2, ], ntrials = colSums(cross.info))
	    
	  }else{
		cdat <- cdat <- data.frame(variant = c(0:1),
		Phenotype = cross.info[2, ], ntrials = colSums(cross.info))
	    
	  }
	  
	  #### IMPORTANT : If any interaction has 0 trials, the exact logistic regression will not run. You have to remove that row before 
	  
	  cdat<-cdat[cdat$ntrials!=0,]  
	  counting<-i-2
	  lastcol<-length(covar_exomes)-2

	  message(paste0("Running test ",counting," of ",lastcol," runs for non adjusted regression"))
	  
	  exact.log.test <- elrm(formula = Phenotype/ntrials ~ variant, interest = ~variant, iter = total.iterations, 
		                 dataset = cdat, burnIn = burnin.its,r=2)
	  
	  p_val.exact.czs<- as.numeric(exact.log.test$p.values)
	  coef.variant<-as.numeric(exact.log.test$coeffs)
	  
	  p_snps_exact_nonadjust[i-2,1]<-colnames(covar_exomes[i])
	  p_snps_exact_nonadjust[i-2,2]<-p_val.exact.czs
	  p_snps_exact_nonadjust[i-2,3]<-exp(coef.variant[1]) #calculating OR
	  p_snps_exact_nonadjust[i-2,4]<-exp(as.numeric(exact.log.test$coeffs.ci[1])) #calculating lower 
	  p_snps_exact_nonadjust[i-2,5]<-exp(as.numeric(exact.log.test$coeffs.ci[2]))#calculating upper
	  
	}

	new_table<-snps.info.subset[,-3]
	colnames(new_table)<-c("CHR","SNP","BP","A1","A2")
	final.result.exact_nonadjust<-cbind(new_table,p_snps_exact_nonadjust)

	#### OUTPUTTING

	write.table(final.result.exact_nonadjust,output.filename,sep="\t", row.names = F, quote = F)

}else{
	message("######		Running chain Logistic regression with adjusment		######")
	
	exome.genotypes<-exome.raw[,c(2,column.2.extract)]
	cov_list <- read.table(cov.filename,header=T,as.is=T)
	covar_exomes <-merge(cov_list, exome.genotypes, by.x = "IID", by.y = "IID",all.y = T)
	snps_columns <- 4:dim(covar_exomes)[2]
	alpha=0.05

	## Creating table
	col.names.ad<-c("TestedAllele","p-value","OR",paste("lower", 100 - 100 *alpha, "ci", sep = ""), paste("upper", 100 - 100 *alpha, "ci", sep = ""))
	p_snps_exact_adjust <- read.table(text = "", col.names = col.names.ad)

	##  running nusing the raw information file

	for (i in snps_columns){ 
		cross.info<-xtabs(~covar_exomes$Phenotype+ interaction(covar_exomes[,i],covar_exomes$Covariate))
		number.alleles<-ncol(cross.info)
		if(number.alleles>6){
			cdatJ <- cdatJ <- data.frame(variant = rep(c(0:2),3),covar = rep(0:2, each = 3),
			Phenotype = cross.info[2, ], ntrials = colSums(cross.info))
		}else{
			cdatJ <- cdatJ <- data.frame(variant = rep(c(0:1),3),covar = rep(0:2, each = 2),
			Phenotype = cross.info[2, ], ntrials = colSums(cross.info))
	}

	#### IMPORTANT : If any interaction has 0 trials, the exact logistic regression will not run. You have to remove that row before 

	cdatJ<-cdatJ[cdatJ$ntrials!=0,]  
	counting<-i-3
	lastcol<-length(covar_exomes)-3

	message(paste0("Running test ",counting," of ",lastcol," runs for adjusted regression"))

	diff.variants<-length(unique(cdatJ$variant))

	if (diff.variants>1){
		exact.log.test.adjust <- elrm(formula = Phenotype/ntrials ~ variant+covar, interest = ~variant, iter = total.iterations, 
		dataset = cdatJ, burnIn = burnin.its,r=4)
		p_val.exact.czs.adjust<- as.numeric(exact.log.test.adjust$p.values)
		coef.variant.adjust<-as.numeric(exact.log.test.adjust$coeffs)
		p_snps_exact_adjust[i-3,1]<-colnames(covar_exomes[i])
		p_snps_exact_adjust[i-3,2]<-p_val.exact.czs.adjust
		p_snps_exact_adjust[i-3,3]<-exp(coef.variant.adjust[1])				#calculating OR
		p_snps_exact_adjust[i-3,4]<-exp(as.numeric(exact.log.test.adjust$coeffs.ci[1]))	#calculating lower 
		p_snps_exact_adjust[i-3,5]<-exp(as.numeric(exact.log.test.adjust$coeffs.ci[2]))	#calculating upper

	}else{

	message(paste0("The variant number ",counting," do not have enough variability"))

	p_snps_exact_adjust[i-3,1]<-"empty"
	p_snps_exact_adjust[i-3,2]<-"empty"
	p_snps_exact_adjust[i-3,3]<-"empty"
	p_snps_exact_adjust[i-3,4]<-"empty"
	p_snps_exact_adjust[i-3,5]<-"empty"

		}
	}

	new_table<-snps.info.subset[,-3]
	colnames(new_table)<-c("CHR","SNP","BP","A1","A2")

	final.result.exact_adjust<-cbind(new_table,p_snps_exact_adjust)

	#### OUTPUTTING

	write.table(final.result.exact_adjust,output.filename,sep="\t", row.names = F, quote = F)
	#write.table(final.result.exact_ajust,paste0(output.path,"microcephaly_exact_model_ajust_january2020.txt"),sep="\t", row.names = F, quote = F)


}


