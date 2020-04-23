#!/usr/bin/Rscript

# Recode a numeric clusterCall tetraploid genotype to GWASpoly ACGT format
# INPUT: Two files: ClusterCall (SNP|---Markers---) genotype and SNPs info (SNP|REF|ALT|XXX|CHROM|XXX|POS)
# Remove duplicated SNPs
USAGE="
Convert clusterCall (numeric) file to GWASpoly (ACGT) file
USAGE: recode-clusterCall-to-GWASpoly.R <ClusterCall file> <SNPs Map File>"

library(parallel)
options(stringsAsFactors = FALSE, width=300)

#----------------------------------------------------------
# main
#----------------------------------------------------------
main <- function () {
	args = commandArgs(trailingOnly = TRUE)
	if (length (args) < 2)
		stop (USAGE)


	#args = c(genotype="agrosavia-genotype-NUM-CLEANED.tbl", snps="solcap-SNPs-refs-alts-AGs.tbl")
	#args = c(genotype="agNUM.tbl", snps="solcap-SNPs-refs-alts-AGs.tbl")
	genotypeFile = args [1]
	SNPsFile     = args [2]

	#convertGenotypeNumericToACGTFormat (genotypeFile, SNPsFile)
	numericToACGTFormatGenotype (genotypeFile, SNPsFile)
}

#----------------------------------------------------------
# Convert ClusterCall numeric format to GWASpoly format
# Basic genotype: [SNPs,...Alleles...] to [SNPs,CHR,POS,...Alleles...]
# Numeric tetra genotype (0,1,2,3,4) to ACGT with chromosome and pos
# Uses solcap ref/alt alleles
# Warning!!! It takes too long
#----------------------------------------------------------
numericToACGTFormatGenotype <- function (genotypeFile, SNPsFile)
{
	genotype     <- read.csv (genotypeFile, header=T, check.names=F)
	SNPsMap      <- read.table (SNPsFile, header=T, check.names=F)
	rownames (SNPsMap) <- SNPsMap [,1]

	genotypeMap = getChromosomeInfoFromMapFile (genotype, SNPsMap)
	alleles     = genotypeMap [,-c(2,3)]  # For GWASpoly numeric

	allelesACGT  <- numericToACGTFormatAlleles (alleles, SNPsMap)
	genoACGT     = data.frame (genotypeMap [,1:3], allelesACGT [,-1])

	outFile = paste0 (strsplit (genotypeFile, split="[.]")[[1]][1],"-ACGT.tbl")
	#msgmsg ("Writing ACGT genotype to ", outFile, "...")
	write.table (file=outFile, genoACGT, quote=F,row.names=F, sep=",")
}

getChromosomeInfoFromMapFile  <- function (genotype, SNPsMap) {
	SNPs         = genotype [,1]
	chromosomes  = SNPsMap [SNPs, 5]
	positions    = SNPsMap [SNPs, 6]
	newGenotype  = cbind (SNP=genotype [,1], CHROM=chromosomes, POS=positions, genotype [,-c(1)])
	return (newGenotype)
}
	

numericToACGTFormatAlleles <- function (alleles, SNPs)
{
	#>>>>>>>> local functions <<<<<<<<<<<<<<<<<
	setA <- function (allelesVec, refs, alts) {
		id  = allelesVec [1]
		gnt <- as.numeric (allelesVec [-1])
		ref = refs [id,2]
		alt = alts [id,2]
		gnt [gnt==4] = strrep(ref,4)
		gnt [gnt==3] = paste0 (strrep(alt,1),strrep(ref,3))
		gnt [gnt==2] = paste0 (strrep(alt,2),strrep(ref,2))
		gnt [gnt==1] = paste0 (strrep(alt,3),strrep(ref,1))
		gnt [gnt==0] = strrep(alt,4)
		return (gnt)
	}
	#>>>>>>>> local functions <<<<<<<<<<<<<<<<<
	refs <- data.frame (SNPs [,c(1,2)])
	rownames (refs) <- SNPs [,1]
	alts <- data.frame (SNPs [,c(1,3)])
	rownames (alts) <- SNPs [,1]
	alls <- alleles

	allelesNUM <- t(apply (alleles,  1, setA, refs, alts ))
	colnames (allelesNUM) = colnames (alleles [-1])
	rownames (allelesNUM) = rownames (alleles)

	newAlleles <- data.frame (Markers=alleles [,1], allelesNUM)
	return (newAlleles)
}

#----------------------------------------------------------
# Call to main
#----------------------------------------------------------
main ()

