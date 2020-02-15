#!/usr/bin/Rscript
# Recode a tetraploid genotype (0,1,2,3) to diploid (0,1,2)
# also recode to ACGT (diplo) using info in solcap scafolds

# Remove duplicated SNPs

options (width=300)
library(parallel)

options(stringsAsFactors = FALSE, width=300)
#----------------------------------------------------------
# Convert numeric (gwaspoly) genotype to ACGT using solcap ref/alt alleles
#----------------------------------------------------------
convertGenotypeNumericToACGTFormat <- function (genotypeFile, SNPsFile) {
	genotype = read.csv (genotypeFile, header=T)
	SNPs     = read.table (SNPsFile, header=T, row.names=1)

	dossages = genotype [,-c(1)]
	markers  = genotype [,1]
	samples  = colnames (dossages)

	print (dossages [1:10,1:10])
	rownames (dossages) = genotype [,1] 

	setAllele <- function (id) {
		dossages [id,which (dossages [id,]==0)] = SNPs[id, "AAAA"]
		dossages [id,which (dossages [id,]==1)] = SNPs[id, "AAAB"]
		dossages [id,which (dossages [id,]==2)] = SNPs[id, "AABB"]
		dossages [id,which (dossages [id,]==3)] = SNPs[id, "ABBB"]
		dossages [id,which (dossages [id,]==4)] = SNPs[id, "BBBB"]
		return (dossages[id,])
	}

	#alleles = matrix (sapply (markers, setAllele), ncol=ncol(dossages), nrow=nrow(dossages), byrow=T)
	alleles = t(sapply (markers, setAllele))
	print (alleles[1:10,1:10])

	newGenotype = cbind (Markers=genotype [,1], alleles)

	outFile = paste0 (strsplit (genotypeFile, split="[.]")[[1]][1],"-tetra-ACGT.tbl")
	message (">>> Writing ACGT genotype to ", outFile, "...")
	write.table (file=outFile, newGenotype, quote=F,row.names=F, sep="\t")
}
#----------------------------------------------------------
#----------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)

args = c(genotype="agrosavia-genotype-NUM-CLEANED.tbl", snps="solcap-SNPs-refs-alts-AGs.tbl")
#args = c(genotype="agNUM.tbl", snps="solcap-SNPs-refs-alts-AGs.tbl")
genotypeFile = args [1]
SNPsFile     = args [2]

convertGenotypeNumericToACGTFormat (genotypeFile, SNPsFile)


