#!/usr/bin/Rscript

# Convert gwaspoly genotype format to plink transposed format (tped, tfam)

library (parallel)

options(stringsAsFactors = FALSE)

#------------------------------------------------------------------------------
## Format and write tassel phenotype
#------------------------------------------------------------------------------
writeTasselPheno <- function (gwaspolyPhenotypeFile) {
	gwaspolyPhenotype = read.table (file=gwaspolyPhenotypeFile, header=T)
	idNames = as.character (gwaspolyPhenotype [,1])
	Taxa = gwaspolyPhenotype [,2]
	BLIGHT = gwaspolyPhenotype [,3]
	plinkPheno = cbind (Taxa, BLIGHT)

	message (">>> Writing tassel phenotype...")
	plinkPhenoFilename = paste0 (strsplit (gwaspolyPhenotypeFile, split="[.]")[[1]][1], "-tassel.tbl")
	sink (plinkPhenoFilename)
	cat ("<Phenotype>\n")
	cat ("taxa\tdata\n")
	write.table (file="", plinkPheno, col.names=T, row.names=F, quote=F, sep="\t")
	sink()
}

#------------------------------------------------------------------------------
## Write and format tassel structure
#------------------------------------------------------------------------------
writeTasselStruct <- function (structureFile) {
	structData = read.csv (file=structureFile, header=T)
	idNames = as.character (structData [,1])
	Trait = sapply (strsplit (idNames, split="And_",fixed=T), function (x) x[[2]][1])
	struct = structData [,-1]
	newData = cbind ("<Trait>"=Trait, struct)

	message (">>> Writing tassel structure")
	outFilename = paste0 (strsplit (structureFile, split="[.]")[[1]][1], "-tassel.tbl")
	sink (outFilename)
	cat ("<Covariate>\n")
	write.table (file="", newData, col.names=T, row.names=F, quote=F, sep="\t")
	sink()
}

#------------------------------------------------------------------------------
## Format and write tassel phenotype
#------------------------------------------------------------------------------
createPlinkPhenotype <- function (gwaspolyPhenotypeFile) {
	phenotype = read.csv (file=gwaspolyPhenotypeFile, header=T)
	idNames = as.character (phenotype [,1])
	# Remove "And_" prefix
	Samples = sapply (strsplit (idNames, split="And_",fixed=T), function (x) x[[2]][1])
	BLIGHT = phenotype [,2]

	message (">>> Writing plink phenotype...")
	plinkPheno = cbind (FID=0,IID=Samples, BLIGHT= BLIGHT)
	outFilename = paste0 (strsplit (gwaspolyPhenotypeFile, split="[.]")[[1]][1], "-plink.tbl")
	write.table (file=outFilename, plinkPheno, col.names=T, row.names=F, quote=F, sep="\t")

	# Normalized phenotype
	message (">>> Writing normalized plink phenotype...")
	MinMaxScaling <- function(x){return((x-min(x))/(max(x)-min(x)))}
	BLIGHT_normalized = round (MinMaxScaling (BLIGHT), 8)
	plinkPheno = cbind (FID=0,IID=Samples, BLIGHT= BLIGHT_normalized)
	outFilename = paste0 (strsplit (gwaspolyPhenotypeFile, split="[.]")[[1]][1], "-plink-norm.tbl")
	write.table (file=outFilename, plinkPheno, col.names=T, row.names=F, quote=F, sep="\t")
}

#----------------------------------------------------------
#----------------------------------------------------------
createFileMAP <- function (genotype, plinkFilename) {
	message (">>> Writing MAP files...")
	markersIds <- gsub ("solcap_snp_", "", genotype [,1])
	chrm    <- genotype [,2]+1
	pos     <- genotype [3]

	genoMAP    = cbind (chr=chrm, iid=markersIds, dist=0, pos=pos)
	filenameMAP = paste0 (plinkFilename, "-plink.map")
	write.table (file=filenameMAP, genoMAP, col.names=F, row.names=F, quote=F, sep="\t")
}

#----------------------------------------------------------
#----------------------------------------------------------
createFilePED <- function (genotype, plinkFilename) {
	message (">>> Writing PED files...")
	#namesGeno  = c("fid", "iid", "pid", "mid", "sex", "phe", rep ("X", ncol (talleles)))
	alleles    <<- genotype [,-c(1,2,3)]

	message (">>> Creating transposed genotype...")
	markersIds <<- gsub ("solcap_snp_", "", genotype [,1] )
	samplesIds <<- gsub ("And_", "", colnames (alleles))
	allelesTransposed <<- t(alleles)
	rownames (allelesTransposed) = samplesIds
	colnames (allelesTransposed) = markersIds

	filenamePEDTransp = paste0 (plinkFilename, "-transposed.tbl")
	write.table (file =filenamePEDTransp, allelesTransposed, col.names=T, row.names=T, quote=F, sep="\t")

	message (">>> Creating PED genotype...")
	talleles   = tabAlleles (alleles)
	message (">>> ...Created PED genotype")
	#talleles   = t (alleles)

	colnames (talleles) = markersIds
	genoPED    = cbind (0, samplesIds, 0,0,0,-9, talleles)

	filenamePED = paste0 (plinkFilename, "-plink.ped")
	write.table (file=filenamePED, genoPED, col.names=F, row.names=F, quote=F, sep="\t")

}

#----------------------------------------------------------
# Add tabs to alleels changign AG --> A	G
#----------------------------------------------------------
tabAlleles <- function (alleles) {
	alleles  [is.na (alleles)] = "00"
	ncols = ncol (alleles)
	nrows = nrow (alleles)
	alleles = t (alleles)
	tabs <- function (x) {return (sprintf("%s %s", substr(x,2,2),substr(x,1,1)))}
	allelesTabbed =matrix (mclapply (alleles,tabs,mc.cores=4), ncol=nrows, nrow=ncols, byrow=F)
	return (allelesTabbed)
}

#----------------------------------------------------------
# Main
#----------------------------------------------------------

args = commandArgs(trailingOnly = TRUE)

args = c("agrosavia-genotype-diplo-AGs.tbl","agrosavia-phenotype.tbl", "structure-checked.tbl")
gwaspolyGenotypeFile  = args [1]
gwaspolyPhenotypeFile = args [2]
structureFile         = args [3]

plinkFilename = strsplit (gwaspolyGenotypeFile, split="[.]")[[1]][1]

# Read the genotype
gwaspolyGenotype = read.table (file=gwaspolyGenotypeFile, header=T)
genotype = gwaspolyGenotype [order (gwaspolyGenotype$Chrom,	gwaspolyGenotype$Position),]

## Write and format tassel phenotype
createPlinkPhenotype (gwaspolyPhenotypeFile)

# Write MAP
createFileMAP (genotype, plinkFilename)

# Write PED
createFilePED (genotype, plinkFilename)

## Write and format tassel structure
writeTasselStruct (structureFile) {


