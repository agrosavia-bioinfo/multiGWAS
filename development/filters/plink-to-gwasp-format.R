#!/usr/bin/Rscript

## Create gwaspoly files (geno and pheno) from plink files (.ped,.map, pheno)

library (dplyr)

options (width=300, stringsAsFactors=T)
args = commandArgs(trailingOnly=T)

#----------------------------------------------------------
# Plink separated alleles (e.g. A C A G A A) are joined (AC AG AA)
#----------------------------------------------------------
joinAlleles <- function (alleles) {
	nCols = ncol (alleles)
	matAGs = matrix (nrow = nrow (alleles), ncol=nCols/2)
	matGAs = matrix (nrow = nrow (alleles), ncol=nCols/2)
	j = 0
	for (i in seq(1,nCols,2)){
		j=j+1
 		matAGs [,j] = cbind (paste0 (alleles [,i], alleles [,i+1]))
 		matGAs [,j] = cbind (paste0 (alleles [,i+1],"\t", alleles [,i]))
	}
	return (list(matAGs=matAGs, matGAs=matGAs))
}

#----------------------------------------------------------
# Create gwaspoly genotype from plink  .ped and .map files
#----------------------------------------------------------
createGwaspolyGenotype <- function (pedFile, mapFile) {
	outName = paste0 (strsplit (pedFile, split="[.]")[[1]][1], "-FLT-gwasp.tbl")

	# Process genotype
	ped = read.table (pedFile, header=F,stringsAsFactors=F, colClasses=c("character"))
	samples = ped [,2]
	alleles = ped [,-c(1,2,3,4,5,6)]
	allelesJoined = joinAlleles (alleles)

	allelesAGs = t(allelesJoined$matAGs)

	map     = read.table (mapFile, header=F)
	chrs    = map [,1]
	markers = as.character (map [,2])
	traits  = map [,3]
	pos     = map [,4]

	message (">>> PED files:")
	print (ped[1:10,1:20])

	message ("\n>>> MAP file:")
	print (map[1:10,])

	#rownames (allelesAGs) = markers
	colnames (allelesAGs) = samples

	genotype = cbind (Markers=markers,Chrom=chrs,Position=pos,allelesAGs)
	gn = genotype 
	#allelesDF = cbind (map [,c(1,2,3,4)], allelesAGs)

	message ("\n>>> Gwaspoly genotype file:")
	print (genotype [1:10,1:10])
	write.table (file=outName, genotype, row.names=F,quote=F, sep=",")

	return (list(samples=as.numeric(samples), markers=markers, alleles=allelesJoined$matGAs, traits=traits))
}

#----------------------------------------------------------
#----------------------------------------------------------
createGwaspolyPhenotype <- function (phenoFile, sampleNames) {
	phenoGwaspFile = gsub (".tbl", "-FLT-gwasp.tbl", phenoFile)

	pheno      <- read.table (phenoFile, header=T)
	#phTmp   <<-  pheno %>% filter (sampleNames %in% IID)
	#phenoDF = pheno %>% arrange (sampleNames %in% IID) %>% as.data.frame()
	#pheno %>% filter (sampleNames %in% IID)

	phenoDF = pheno %>% filter (IID %in% sampleNames) %>% arrange (match (IID, sampleNames)) %>% data.frame()

	phenoGwasp <- cbind (Samples=phenoDF[,2],BLIGHT=phenoDF[,3])

	message ("\n>>> Gwaspoly penotype file:")
	print (head (phenoGwasp))
	write.table (file=phenoGwaspFile, phenoGwasp, row.names=F,quote=F,sep=",")

	return (phenoDF[,3])

}

#----------------------------------------------------------
# Filter the tetra gwaspoly genotype using the plink filtered genotype
#----------------------------------------------------------
filterTetraGwaspolyWithPlinkGeno <- function (gwaspGenoTetraFile, samples, markers) {
	gwp = read.csv (file=gwaspGenoTetraFile, header=T, check.names=F)
	gwp [,2] = gwp [,2]+1

	gwpDF = gwp %>% select (c(1,2,3), as.character(samples)) %>% filter (Markers %in% markers) %>% 
			arrange (match(Markers, markers)) %>% data.frame(check.names=F)

	outName = gsub (".tbl", "-FLT.tbl", gwaspGenoTetraFile)
	write.csv (file=outName, gwpDF, quote=F,row.names=F)
}

#----------------------------------------------------------
#----------------------------------------------------------
filterStructureGwaspoly <- function  (structFile) {
	struct = read.csv (file=structFile, header=T, check.names=F)

	structDF = struct %>% filter (Samples %in% samples) %>% arrange (match (Samples, samples)) %>% data.frame (check.names=F)
	outName = gsub (".tbl", "-FLT.tbl", structFile)
	write.csv (file=outName, structDF, quote=F,row.names=F)
}


#----------------------------------------------------------
# Main
#----------------------------------------------------------
args = c("geno-6.ped", "geno-6.map", "pheno.tbl", "agrosavia-genotype-tetra-NUM.tbl", "agrosavia-structure.tbl")

plinkPEDFile  = args [1]
plinkMAPFile  = args [2]
phenoFile  = args [3]
gwaspGenoTetraFile = args [4]
structFile = args [5]


genotype = createGwaspolyGenotype (plinkPEDFile, plinkMAPFile)
samples  = genotype$samples
markers  = genotype$markers
alleles  = genotype$alleles
traits   = genotype$traits
#
traits = createGwaspolyPhenotype (phenoFile, genotype$samples)
#
filterTetraGwaspolyWithPlinkGeno (gwaspGenoTetraFile, samples, markers)
#
filterStructureGwaspoly (structFile)

##genoFile  = "agrosavia-genotype-FLT-plink.ped"
##phenoFile = "agrosavia-phenotype-FLT.plink"
##
##geno  = read.table (file=genoFile, header=F,colClasses= c("character"))
##pheno = read.table (file=phenoFile, header=T)
##
##ng = cbind (IID=paste0("A",pheno$IID), BLIGHT=pheno$BLIGHT, geno [,-c(1,2,3,4,5,6)])
##outName = paste0 (strsplit (genoFile, split="-FLT-")[[1]][1], "-FLT-shesis.tbl")
##write.table (file=outName, ng, row.names=F, col.names=F, quote=F, sep="\t")
##
##snpsNames  = geno [,-c(1,2,3,4,5,6)]
##outName = paste0 (strsplit (genoFile, split="-FLT-")[[1]][1], "-FLT-shesis-names.txt")
##write.csv (file=outName, snpsNames, row.names=F, quote=F)


#---------------------
outFile = paste0 (strsplit (plinkPEDFile, split="[.]")[[1]][1], "-FLT-shesis.tbl")
genoShesis = cbind (paste0("A", samples),traits,alleles)
write.table (file=outFile, genoShesis, row.names=F, col.names=F, quote=F, sep="\t")

outFile       = paste0 (strsplit (plinkPEDFile, split="-FLT-")[[1]][1], "-FLT-shesis-names.txt")
con = file (outFile)
writeLines (markers, con)
close (con)

