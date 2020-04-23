#!/usr/bin/Rscript

#-------------------------------------------------------------
# Print a log message with the parameter
#-------------------------------------------------------------
msg <- function (...) {
  messages = unlist (list (...))
  cat (">>>>", messages, "\n")
}

options (width=300)
args = commandArgs(trailingOnly = TRUE)
args = c ("original/papa-agrosavia-fenotipoGota.txt", 
		  "original/ClusterCall_prediction_CCC-fixed.txt", 
		  "original/K5_estructura_tetraploides_2017.txt",
		  "original/potato_8303SNPs_potato_dm_v4.03.gff3")

phenotypeFile = args [1]
genotypeFile  = args [2]
structureFile = args [3]
genomeFile    = args [4]
#-------------------------------------------------------------

	msg ("Reading files and checking colnames and rownames...")
	phenotypeAll = read.table (phenotypeFile, header=T, sep="\t")
	genotypeAll  = read.table (genotypeFile, header=T, sep=",")
	structureAll = read.table (structureFile, header=T)
	genomeAll    = read.table (genomeFile, header=F)

	# Get samples names for all
	samplesPheno  = phenotypeAll$Parcela
	samplesGeno   = colnames (genotypeAll)
	samplesStruct = structureAll$sample

	#---------------------------------------------------------
	# Unify samples
	#---------------------------------------------------------
	commonSamples = Reduce (intersect, list (samplesGeno, samplesPheno, samplesStruct))

	phenotype  = phenotypeAll [phenotypeAll$Parcela %in% commonSamples, c("Parcela", "GE")]
	genotype   = genotypeAll  [,colnames(genotypeAll) %in% c("Markers", commonSamples)]
	structure  = structureAll [structureAll[,1] %in% commonSamples,]
	colnames (phenotype) = c ("sample", "gota")
	genotypeStructure = merge (phenotype, structure, by="sample")

	# Write checked files checked sample names
	write.table (file="phenotype-checked.tbl", phenotype, row.names=F, quote=F, sep=",")
	write.table (file="structure-checked.tbl", structure, row.names=F, quote=F, sep=",")
	write.table (file="genostruct-checked.tbl", genotypeStructure, row.names=F, quote=F, sep=",")
	
	#---------------------------------------------------------
	# Add Structure to phenotype
	#---------------------------------------------------------

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
	genomeChromosomes = unlist (lapply (as.vector (genomeAll [,c(1)]), getChromMarker))
	genomePositions   = unlist (genomeAll [,c(4)])

	# Create new genome
	genome    = cbind (genomeMarkers, genomeChromosomes, genomePositions)
	colnames (genome) = c("Markers","Chrom","Position")

	# Merge genome with genotype
	genotypeChrom = merge (genome, genotype, by="Markers")
	write.table (file="genotype-checked.tbl", genotypeChrom, row.names=F, quote=F, sep=",")

