#!/usr/bin/Rscript
"
GOAL:    Function to run Naive and Full GWAS in TASSEL using rtassel R library 
AUTHOR:  Luis Garreta
VERSION: r1.0: Naive GWAS
"

options(java.parameters = c("-Xmx1g", "-Xms1g"))
library(rTASSEL)

rTASSEL::startLogger(fullPath = getwd(), fileName = "out-TASSEL.log")

geno <- "tassel-geno.vcf"
pheno="phenot.tbl"

#----------------------------------------------------------
#----------------------------------------------------------
main <- function () {
	args = commandArgs (trailingOnly=T)
	geno  = args [1]
	pheno = args [2]
	runNaiveTASSEL (geno, pheno)
}

#----------------------------------------------------------
# Perform association analysis using General Linear Model 
#----------------------------------------------------------
runNaiveTASSEL <- function (genotype, phenotype, outFile) {
	msg ("Reading TASSEL Genotype/Phenotype ")
	tasGenoPheno <- rTASSEL::readGenotypePhenotype(
		genoPathOrObj    = genotype,
		phenoPathDFOrObj = phenotype
	)

	msg ("    >>>> Running TASSEL Naive (GLM )...")
	tasGLM <- rTASSEL::assocModelFitter(tasObj=tasGenoPheno, formula=.~., fitMarkers=T, kinship=NULL, fastAssociation=F)
	#
	print (tasGLM)
	write.table (tasGLM$GLM_Stats, file=outFile, quote=F, col.names=T, row.names=F, sep="\t")
}

#-------------------------------------------------------------
# Print a log message with the parameter
#-------------------------------------------------------------
msg <- function (...) 
{
  messages = unlist (list (...))
  cat ("\t>>>>", messages, "\n")
}

#----------------------------------------------------------
#----------------------------------------------------------
#main()
