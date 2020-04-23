#!/usr/bin/Rscript

# Get from Genome information file (.gff3) the SNP, CHROM, POS
USAGE="getSNPsInfo.R <genome file>  ---> SNPs File"

options (width=300)
args = commandArgs (trailingOnly=T)

main <- function () {
	genomeFile    = args [1]
	getSnpChromosomePositionFromGenome (genomeFile)
}

#---------------------------------------------------------
# Get from Genome information file (.gff3) the SNP, CHROM, POS
# Write results to table "genome-SNPs-Chrom-Pos.tbl"
#---------------------------------------------------------
getSnpChromosomePositionFromGenome <- function (genomeFile) {
	#>>>>>>>>>>> local function <<<<<<<<<<<<<<<
	# Get SNP id from text string
	getIdMarker <- function (s) {
		return (strsplit (strsplit (s, split="=")[[1]][3], split=";")[[1]][1])
	}
	# Remove prefix "chr" from "chrXX"
	getChromMarker <- function (s) {
		return (strsplit (s, split="chr")[[1]][2])
	}
	#>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<

	# Get important columns from genome: CHRM, POS, SNP Info
	genomeAll    = read.table (genomeFile, header=F)
	genomeTmp = genomeAll [,c(1, 4, 9)]

	# Parsing of markers ids and chromosomes value
	genomeMarkers     = unlist (lapply (as.vector (genomeAll [,c(9)]), getIdMarker))
	#genomeMarkers     = gsub ("solcap_snp_","", genomeMarkers)  # To remove long prefix
	genomeChromosomes = as.numeric (unlist (lapply (as.vector (genomeAll [,c(1)]), getChromMarker)))
	genomePositions   = unlist (genomeAll [,c(4)])

	# Create new genome and write
	genome       = data.frame (SNP=genomeMarkers, CHROM=genomeChromosomes+1, POS=genomePositions)

	#gf=genotypeFiltered = gc [!duplicated (gc[,1],),]
	genomeNoDups = genome [!duplicated (genome$SNP),]
	genomeNoDups = genomeNoDups [order (genomeNoDups$SNP),]
	write.table (file="genome-SNPs-Chrom-Pos.tbl", genomeNoDups, quote=F, sep="\t", row.names=F)
}

#----------------------------------------------------------
# Call to main function
#----------------------------------------------------------
main ()


