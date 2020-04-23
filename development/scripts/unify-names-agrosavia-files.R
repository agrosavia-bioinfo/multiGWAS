#!/usr/bin/Rscript

# Unify genotype and phenotype to same samples and markers, 
# It completes the genotype with chromosomes info from genome

#-------------------------------------------------------------
# Print a log message with the parameter
#-------------------------------------------------------------
msg <- function (...) {
  messages = unlist (list (...))
  cat (">>>>", messages, "\n")
}

options (width=300)
args = commandArgs(trailingOnly = TRUE)
args = c ("original/agrosavia-phenotype-Gota-SanJorge1997-CLEANED.tbl", 
		  "original/agrosavia-genotype-CCC-CLEANED.tbl", 
		  "original/potato_8303SNPs_potato_dm_v4.03.gff3")

phenotypeFile = args [1]
genotypeFile  = args [2]
genomeFile    = args [3]
#-------------------------------------------------------------

	msg ("Reading files and checking colnames and rownames...")
	pa=phenotypeAll = read.csv (phenotypeFile, header=T, sep=",")
	ga=genotypeAll  = read.csv (genotypeFile, header=T, sep=",")
	genomeAll    = read.table (genomeFile, header=F)

	# Get samples names for all
	sp=samplesPheno  = phenotypeAll$Parcela
	sg=samplesGeno   = colnames (genotypeAll)

	#---------------------------------------------------------
	# Unify samples
	#---------------------------------------------------------
	cs=commonSamples = Reduce (intersect, list (samplesGeno, samplesPheno))

	phenotype  = phenotypeAll [phenotypeAll$Parcela %in% commonSamples, c("Parcela", "Gota")]
	genotype   = genotypeAll  [,colnames(genotypeAll) %in% c("Markers", commonSamples)]
	colnames (phenotype) = c ("Parcela", "Gota")

	# Write checked files checked sample names
	write.table (file="agrosavia-phenotype-checked.tbl", phenotype, row.names=F, quote=F, sep=",")
	
	#---------------------------------------------------------
	# Add Chomosome and position from genome to genotype #####
	#---------------------------------------------------------
	# Gen markers from genome, writing their markers names as the id at first col
	genomeTmp = genomeAll [,c(1, 4, 9)]
	getIdMarker <- function (s) {
		return (strsplit (strsplit (s, split="=")[[1]][3], split=";")[[1]][1])
	}
	getChromMarker <- function (s) {
		return (strsplit (s, split="chr")[[1]][2])
	}

	# Parsing of markers ids and chromosomes value
	genomeMarkers     = unlist (lapply (as.vector (genomeAll [,c(9)]), getIdMarker))
	genomeChromosomes = as.numeric (unlist (lapply (as.vector (genomeAll [,c(1)]), getChromMarker)))
	genomePositions   = unlist (genomeAll [,c(4)])

	# Create new genome
	genome    = data.frame (gsub ("solcap_snp_","", genomeMarkers), genomeChromosomes+1, genomePositions)
	colnames (genome) = c("Markers","Chrom","Position")

	# Merge genome with genotype
	gc=genotypeChrom = merge (genome, genotype, by="Markers")
	gf=genotypeFiltered = gc [!duplicated (gc[,1],),]
	write.table (file="agrosavia-genotype-checked.tbl", genotypeFiltered, row.names=F, quote=F, sep=",")

