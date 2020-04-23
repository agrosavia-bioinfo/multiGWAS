#!/usr/bin/Rscript

# LOG: r1.0: Format for plink and gwas using a default "out" filenames

#options (width=300, stringsAsFactors=F)
#args = commandArgs(trailingOnly=T)
library (parallel)
library (stringi)
suppressMessages (library (dplyr))
formatsLogFile="log-formats.log"

#------------------------------------------------------------------------------
## Format gwaspoly phenotype to plink format
#------------------------------------------------------------------------------
gwasp2plinkPhenotype <- function (gwaspPhenoFile, outFile="") 
{
	msg ("    >>>> Creating plink phenotype...")
	phenotype = read.csv (file=gwaspPhenoFile, header=T)
	idNames = as.character (phenotype [,1])

	samples = phenotype [,1]
	traits  = phenotype [,2]

	msg ("    >>>> Writing plink phenotype...")
	plinkPheno = data.frame (FID=0,IID=samples, TRAIT= traits)
	if (outFile=="")
		outFile = paste0 ("out/",strsplit (gwaspPhenoFile, split="[.]")[[1]][1], "-plink.tbl")
	write.table (file=outFile, plinkPheno, col.names=T, row.names=F, quote=F, sep="\t")
}

#------------------------------------------------------------------------------
## Format and write gwasp to tassel phenotype
#------------------------------------------------------------------------------
gwasp2tasselPhenotype<- function (gwaspPhenotypeFile, outFilename="") 
{
	gwaspPhenotype <- read.csv (file=gwaspPhenotypeFile, header=T)
	Taxa <- as.character (gwaspPhenotype [,1])
	TRAIT <- gwaspPhenotype [,2]
	tasselPheno <- cbind (Taxa, TRAIT)

	msg ("    >>>> Writing gwasp to tassel phenotype...")
	if (outFilename=="")
		outFilename = paste0 (strsplit ("out/",gwaspPhenotypeFile, split="[.]")[[1]][1], "-tassel.tbl")

	sink (outFilename)
	cat ("<Phenotype>\n")
	cat ("taxa\tdata\n")
	write.table (file="", tasselPheno, col.names=T, row.names=F, quote=F, sep="\t")
	sink()
}

#------------------------------------------------------------------------------
## Create plink MAP file from gwaspoly genotype 
#------------------------------------------------------------------------------
gwasp2plinkGenoMap <- function (gwaspGenoFile) 
{
	msg ("    >>>> Creating plink MAP file...")
	genotype    <- read.table (file=gwaspGenoFile, header=T,stringsAsFactors=T,sep=",")
	markers     <- as.character(genotype [,1])
	chromosomes <- genotype [,2]
	positions   <- genotype [,3]

	plinkMap     <- data.frame (chr=chromosomes, iid=markers, dist=0, positions=positions)
	plinkMapSorted <- plinkMap %>% arrange (chr, positions)
	outFile   = paste0 ("out/",strsplit (basename (gwaspGenoFile), split="[.]")[[1]][1], "-plink.map")
	write.table (file=outFile, plinkMapSorted, col.names=F, row.names=F, quote=F, sep="\t")
	return (plinkMapSorted$iid)
}

#----------------------------------------------------------
# Getting ref/alt alleles for SNPs
#----------------------------------------------------------
bases <- c("A","C","G","T")
get.ref <- function(x) 
{
	y <- paste(na.omit(x),collapse="")
	ans <- apply(array(bases),1,function(z,y){length(grep(z,y,fixed=T))},y)
	if (sum(ans)>2) {stop("Error in genotype matrix: More than 2 alleles")}
	if (sum(ans)==2) {ref.alt <- bases[which(ans==1)]}
	if (sum(ans)==1) {ref.alt <- c(bases[which(ans==1)],NA)}

	return(ref.alt)
}
#----------------------------------------------------------
# Create plink PED file from gwaspoly genotype 
#----------------------------------------------------------
gwaspTetraGenoToPlinkPed <- function (gwaspGenoFile, markersIdsMap) 
{
	plinkFile  = paste0 ("out/", strsplit (basename(gwaspGenoFile), split="[.]")[[1]][1], "-plink")

	if (file.exists (paste0(plinkFile,".ped"))) {
		msg ("    >>>> Loading plink file...")
		return (plinkFile)
	}else {
		msg ("    >>>> Creating plink PED file...")
		genotype   = read.csv (file=gwaspGenoFile, header=T)
		alleles    <- as.matrix (genotype [,-c(1,2,3)])
		rownames (alleles) <- genotype [,1]

		msg ("    >>>> Creating transposed genotype...")
		markersIds        <- genotype [,1] 
		samplesIds        <- colnames (alleles)

		msg ("    >>>> Getting Ref/Alt Alleles...")
		refAltAlleles <- apply(alleles,1,get.ref)

		msg ("    >>>> Converting tetra to diplo")
		allelesDiplo  <- tetraToDiplos (alleles, refAltAlleles)
		#allelesDiplo  <- old_tetraToDiplos (alleles, refAltAlleles)
		rownames (allelesDiplo) = markersIds
		colnames (allelesDiplo) = samplesIds

		# Adjust for plink PED file
		allelesPlink <- t(allelesDiplo[markersIdsMap,])
		genoPED    <- cbind (0, samplesIds, 0,0,0,-9, allelesPlink)

		msg ("    >>>> Writing plink diplo PED file to ", plinkFile)
		plinkFilePed  = paste0 (plinkFile, ".ped")
		write.table (file=plinkFilePed, genoPED, col.names=F, row.names=F, quote=F, sep="\t")

		# Other files
		# Write tetra to diplo for gwasp2
		gwasp2genotype = cbind (genotype [,1:3], allelesDiplo)
		outFile   = paste0 ("out/",strsplit (basename(gwaspGenoFile), split="[.]")[[1]][1], "-diplo.tbl")
		write.table (file=outFile, gwasp2genotype, col.names=T, row.names=F, quote=F, sep=",")
		runCommand (paste0 ("ln -s ", getwd(),"/",outFile, " out/filtered-gwasp2-genotype.tbl"), formatsLogFile)

		# Transposed gwasp file
		transposedAlleles <- t(alleles)
		rownames (transposedAlleles) = samplesIds
		colnames (transposedAlleles) = markersIds
		write.table (file ="out/tmp-plink-transposed.tbl", transposedAlleles, col.names=T, row.names=T, quote=F, sep="\t")


		# Write diplo matrix for original matrix
		write.table (file="tmp-diplosMatrix.tbl",allelesDiplo,  quote=F, sep="\t")
	}
	
	return (plinkFile)
}

#----------------------------------------------------------
#----------------------------------------------------------
gwasp2shesisGenoPheno <- function (genotypeFile, phenotypeFile) 
{
	msg ("    >>>> Writting gwasp to shesis genopheno...")
	sep <- function (allele) {
		s="";
		for (i in 1:4) s=paste0(s, substr (allele,start=i,stop=i)," ");
		return (s)
	}

	geno    <<- read.csv (file=genotypeFile)
	pheno   <<- read.csv (file=phenotypeFile)
	rownames (pheno) <- pheno [,1]
	rownames (geno)  <- geno [,1]

	alleles    <- geno[,-c(1,2,3)]
	allelesMat <- t(sapply (alleles, sep))

	samples         = rownames (allelesMat)
	markers         = colnames (allelesMat)
	pheno           = pheno [samples,]
	genoPhenoShesis = data.frame (Sample=pheno[,1], Trait=pheno[,2],  allelesMat)

	msg ("    >>>> Writing shesis genopheno...")
	outFile = "out/filtered-shesis-genopheno.tbl"
	write.table (file=outFile, genoPhenoShesis, quote=F,row.names=F,col.names=F, sep="\t")

	msg ("    >>>> Writing shesis marker names...")
	outFile = "out/filtered-shesis-markernames.tbl"
	write.table (file=outFile, rownames(geno), quote=F,row.names=F,col.names=F, sep="\t")
}

#----------------------------------------------------------
# Add tabs to alleels changign AG --> A	G
#----------------------------------------------------------
tetraToDiplos <- function (alleles, refAltAlleles) 
{
	if (file.exists ("tmp-diplosMatrix.tbl")) {
		msg ("    >>>> Loading diplos matrix...")
		allelesMat = as.matrix (read.table ("tmp-diplosMatrix.tbl"))
	}
	else {
		msg ("    >>>> Calculating diplos matrix...")
		msg ("    >>>> Converting tetra to diplos...")
		t2d <- function (allele, ref, alt, snp) {
			if (is.na (allele)  | allele=="0000") return ("00")
			else if (strrep (ref,4) == allele) return (paste0(ref,ref))
			else if (strrep (alt,4) == allele) return (paste0(alt,alt))
			else return (paste0(ref,alt))
		}

		ref <- refAltAlleles [1,]
		alt <- refAltAlleles [2,]

		allelesLst  <- mcmapply (t2d, alleles, ref, alt, rownames (alleles), mc.cores=4)

		allelesMat <- matrix (unlist(allelesLst), ncol=ncol(alleles), 
							   nrow=nrow(alleles), byrow=F)
	}
	return (allelesMat)
}

#-------------------------------------------------------------
# Print a log message with the parameter
#-------------------------------------------------------------
msg <- function (...) 
{
  messages = unlist (list (...))
  cat (">>>>", messages, "\n")
}

#-------------------------------------------------------------
#-------------------------------------------------------------
runCommand <- function (command, logFile="") {
	msg (">>>> ", command)
	if (logFile != "")
		system (paste0 (command, " > ", logFile))
	else
		system (command)
}
#----------------------------------------------------------
# Main
#----------------------------------------------------------
 main <- function () 
 {
	args = c ("agrosavia-genotype-tetra-ACGT.tbl", "agrosavia-phenotype.tbl")

	gwaspGenoFile  = args [1]
	gwaspPhenoFile = args [2]

	message (">>> Converting gwaspoly to plink formats...")
	gwasp2plinkPhenotype (gwaspPhenoFile)
	markersIdsMap = gwasp2plinkGenoMap (gwaspGenoFile)
	gwaspTetraGenoToPlinkPed (gwaspGenoFile, markersIdsMap)
}

