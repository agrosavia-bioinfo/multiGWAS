#!/usr/bin/Rscript

options(java.parameters = c("-Xmx1g", "-Xms1g"))
library(rTASSEL)

rTASSEL::startLogger(fullPath = NULL, fileName = NULL)

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
#----------------------------------------------------------
runNaiveTASSEL <- function (genotype, phenotype) {
	# Load in hapmap file
	message (">>>>> Genotype")
	tasGenoHMP <- rTASSEL::readGenotypeTableFromPath(path=geno)

	# Load into rTASSEL "GenotypePhenotype" object
	message (">>>>> Phenotype")
	tasPheno <- rTASSEL::readPhenotypeFromPath(path=pheno)

	message (">>>> Geno/Pheno")
	tasGenoPheno <- rTASSEL::readGenotypePhenotype(
		genoPathOrObj = tasGenoHMP,
		phenoPathDFOrObj = tasPheno
	)

	print (tasGenoPheno)

	message ("Calculating GLM...")
	tasGLM <- rTASSEL::assocModelFitter(tasObj=tasGenoPheno, formula=.~., fitMarkers=T, kinship=NULL, fastAssociation=F)
	#
	print (tasGLM)
	write.table (tasGLM$GLM_Stats, file="tasse_glm.scores", quote=F, col.names=T, row.names=F, sep="\t")
}

#----------------------------------------------------------
#----------------------------------------------------------
main()
