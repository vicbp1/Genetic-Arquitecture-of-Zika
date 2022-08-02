#==========================================
# (C) Copyright 2020, by LNCC and Contributors.
#
# -----------------
#  Script to run FIRTH's LOGISTIC REGRESSION
#  -----------------
#
# Author: Victor Borda
## 
## With this script you can run an Firth logistic regression. As output you will have the description of the variant, p-value, OR and OR ci
## !!!! IMPORTANT: MUST HAVE PACKAGE "logistf" INSTALLED !!!!
##
##  To run Running_LOGISTF_ver2.R
##  R < ELRrun.R <raw.plink.filename> <plink bim file>  <model A | AD>  <output.filename> <covar file> --no-save
## example: R < ELRrun.R zikadata.raw 22 <output.filename> --no-save

temp=commandArgs()

raw.infile=as.character(temp[2])
name.map=as.character(temp[3])
model=as.character(temp[4])
output.filename=as.character(temp[5])
outcome=as.character(temp[6])
covars=as.character(temp[7])
cov.filename=as.character(temp[8])
adjornot=as.character(temp[9])

print(covars)
covars=c(unlist(strsplit(covars,",")))
print(covars)

library("SKAT")
library("SNPRelate")
library("plyr")
library("logistf")

### Reading files
#name.map<-"/home/victor/LNCC/Exomes_zika/final_calling_March2020/Logistic/Firth_Logistic/by_chrs/dataset_AD_chr22.bim"
snps.info <- read.table(name.map,header=FALSE,na.strings=c(""),colClasses = c("numeric","character","numeric","numeric","character","character"))

######## Subsetting input depending of the model ########


#raw.infile<-"/home/victor/LNCC/Exomes_zika/final_calling_March2020/Logistic/Firth_Logistic/by_chrs/dataset_AD_chr22.raw"
if (isTRUE(model=="A")){
  exome.raw <- read.table(raw.infile,header=T)
}else{
  exome.AD <- read.table(raw.infile,header=T)
  columnsAD <- grep("HET", colnames(exome.AD))
  exome.AD<-exome.AD[,c(1:6,columnsAD)]
  exome.2.remove.novar<-exome.AD[,-c(1:6)]
  columns_non_unique<-exome.2.remove.novar[vapply(exome.2.remove.novar, function(x) length(unique(x)) > 1, logical(1L))]
  names2extract<-colnames(columns_non_unique)
  positions<-c()
  for (i in c(1:length(names2extract))){
    x<-which(colnames(exome.AD)==names2extract[i] )
    positions<-c(positions,x)
  }
  exome.raw<-exome.AD[,c(1:6,positions)]
  
}

if (isTRUE(adjornot=="NO")){
  
  exome.genotypes<-exome.raw[,-c(1,3,4,5,6)]
  cov_list <- read.table(cov.filename,header=T,as.is=T)
  covar_exomes <-merge(cov_list, exome.genotypes, by.x = "IID", by.y = "IID",all.y = T)
  snps_columns <- 20:dim(covar_exomes)[2]                     ### Depending on the number of columns corresponding to covariates
  alpha=0.05
  
  col.names.nonad<-c("TestedAllele","p-value","OR",paste("lower", 100 - 100 *alpha, "ci", sep = ""), paste("upper", 100 - 100 *alpha, "ci", sep = ""))
  p_snps_nonadjust <- read.table(text = "",col.names = col.names.nonad)
  
  ####################################### REGRESSION #################################################
  
  for (i in snps_columns){
    
    counting<-i-19
    lastcol<-length(covar_exomes)-19
    
    message(paste0("running variant ",counting," of ",lastcol," runs"))
    zika_nonadj<-logistf(data=covar_exomes,covar_exomes$Phenotype_CZS~covar_exomes[,i],firth = TRUE, pl=TRUE)
    
    p_snps_nonadjust[i-19,1]<-colnames(covar_exomes[i])
    p_snps_nonadjust[i-19,2]<-as.numeric(zika_nonadj$prob[2])
    p_snps_nonadjust[i-19,3]<-exp(zika_nonadj$coefficients[2])                #calculating OR
    p_snps_nonadjust[i-19,4]<-exp(as.numeric(zika_nonadj$ci.lower[2]))        #calculating lower 
    p_snps_nonadjust[i-19,5]<-exp(as.numeric(zika_nonadj$ci.upper[2]))        #calculating upper
    
  }
    rows2keep<-positions-6
    new_table<-snps.info[rows2keep,-3] ### do not copy the thrid column
    colnames(new_table)<-c("CHR","SNP","BP","A1","A2")
    
    final.result_nonadjust<-cbind(new_table,p_snps_nonadjust)

    ####### OUTPUTTING
    
    write.table(final.result_nonadjust,output.filename,sep="\t", row.names = F, quote = F)
    
}else{
  exome.genotypes<-exome.raw[,-c(1,3,4,5,6)]
  cov_list <- read.table(cov.filename,header=T,as.is=T)
  covar_exomes <-merge(cov_list, exome.genotypes, by.x = "IID", by.y = "IID",all.y = T)
  snps_columns <- 20:dim(covar_exomes)[2]                     ### Depending on the number of columns corresponding to covariates
  alpha=0.05
  
  col.names.table.adjust<-c("TestedAllele","p-value","OR",paste("lower", 100 - 100 *alpha, "ci", sep = ""), paste("upper", 100 - 100 *alpha, "ci", sep = ""))
  p_snps_adjust <- read.table(text = "", col.names = col.names.table.adjust)
  
  ####################################### REGRESSION #################################################

  
  for (i in snps_columns){
    
    counting<-i-19
    lastcol<-length(covar_exomes)-19
    
    message(paste0("running variant ",counting," of ",lastcol," runs"))
    
    elements<-paste("covar_exomes[,i]",paste(covars,collapse = "+"),sep = "+")
    f <- as.formula(
      paste(outcome, 
            elements, 
            sep = " ~ "))
    zika_adj   <-logistf(data=covar_exomes,f,control=logistf.control(maxit=150),firth = TRUE, pl=TRUE)
    
    p_snps_adjust[i-19,1]<-colnames(covar_exomes[i])
    p_snps_adjust[i-19,2]<-as.numeric(zika_adj$prob[2])
    p_snps_adjust[i-19,3]<-exp(as.numeric(zika_adj$coefficients[2]))                #calculating OR
    p_snps_adjust[i-19,4]<-exp(as.numeric(zika_adj$ci.lower[2])) #calculating lower 
    p_snps_adjust[i-19,5]<-exp(as.numeric(zika_adj$ci.upper[2]))#calculating upper
    
    
  }
  rows2keep<-positions-6
  new_table<-snps.info[rows2keep,]
  colnames(new_table)<-c("CHR","SNP","BP","A1","A2")
  final.result_adjust<-cbind(new_table,p_snps_adjust)
  
  #### OUTPUTTING

  write.table(final.result_adjust,output.filename,sep="\t", row.names = F, quote = F)

}

