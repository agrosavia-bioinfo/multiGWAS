#!/usr/bin/Rscript

# INFO  : Different functions to convert between different genotype/phenotype formats
# AUTHOR: Luis Garreta (lgarreta@agrosavia.co)
# DATE  : 12/Feb/2020
# LOG   : 
	# r1.0: Used by MultiGWAS tool
	# r1.1: Modified some function names

#options (width=300, stringsAsFactors=F)
#args = commandArgs(trailingOnly=T)
suppressMessages (library (parallel))
suppressMessages (library (dplyr))
formatsLogFile="log-formats.log"

#------------------------------------------------------------------------------
# Convert genotye from plink (.ped, .map) to VCF (Variant Call Format)
#------------------------------------------------------------------------------
plinkToVCFFormat <- function (plinkFile, outFile) {
	#plinkPrefix = strsplit (plinkFile, split="[.]")[[1]][1]
	cmm = sprintf ("plink --file %s --recode vcf-fid --out %s", plinkFile, outFile)
	runCommand (cmm)
}

#------------------------------------------------------------------------------
## Format and write gwasp to tassel phenotype (For rtassel)
#------------------------------------------------------------------------------
gwaspToTasselPhenotype<- function (gwaspPhenotypeFile, outFilename="") 
{
	gwaspPhenotype = read.csv (file=gwaspPhenotypeFile, header=T, check.names=F, )
	taxa           = as.character (gwaspPhenotype [,1])
	trait          = gwaspPhenotype [,2]
	traitName      = colnames (gwaspPhenotype)[2]
	tasselPheno    = cbind (taxa, trait)

	colnames (tasselPheno) = c ("<Trait>", traitName)

	#msgmsg ("    >>>> Writing gwasp to tassel phenotype...")
	if (outFilename=="")
		outFilename = paste0 (strsplit ("out/",gwaspPhenotypeFile, split="[.]")[[1]][1], "-tassel.tbl")

	write.table (file=outFilename, tasselPheno, col.names=T, row.names=F, quote=F, sep="\t")
}

##------------------------------------------------------------------------------
### Format and write gwasp to tassel phenotype (For pipeline)
##------------------------------------------------------------------------------
#gwasp2tasselPhenotype<- function (gwaspPhenotypeFile, outFilename="") 
#{
#	gwaspPhenotype <- read.csv (file=gwaspPhenotypeFile, header=T, check.names=F)
#	Taxa <- as.character (gwaspPhenotype [,1])
#	TRAIT <- gwaspPhenotype [,2]
#	tasselPheno <- cbind (Taxa, TRAIT)
#
#	msgmsg ("    >>>> Writing gwasp to tassel phenotype...")
#	if (outFilename=="")
#		outFilename = paste0 (strsplit ("out/",gwaspPhenotypeFile, split="[.]")[[1]][1], "-tassel.tbl")
#
#	sink (outFilename)
#	cat ("<Phenotype>\n")
#	cat ("taxa\tdata\n")
#	write.table (file="", tasselPheno, col.names=T, row.names=F, quote=F, sep="\t")
#	sink()
#}



#------------------------------------------------------------------------------
## Format gwaspoly phenotype to plink format
#------------------------------------------------------------------------------
gwasp2plinkPhenotype <- function (gwaspPhenoFile, outFile="") 
{
	#msgmsg ("    >>>> Creating plink phenotype...")
	phenotype = read.csv (file=gwaspPhenoFile, header=T, check.names=F)
	idNames = as.character (phenotype [,1])

	samples = phenotype [,1]
	traits  = phenotype [,2]

	#msgmsg ("    >>>> Writing plink phenotype...")
	plinkPheno = data.frame (FID=samples,IID=samples, TRAIT= traits)
	if (outFile=="")
		outFile = paste0 ("out/",strsplit (gwaspPhenoFile, split="[.]")[[1]][1], "-plink.tbl")
	write.table (file=outFile, plinkPheno, col.names=T, row.names=F, quote=F, sep="\t")
}

#------------------------------------------------------------------------------
# Gwasp to plink format (.ped .map)
#------------------------------------------------------------------------------
gwaspToPlinkFormat <- function (genotypeFile, plinkFile) {
	markersIdsMap = createPlinkMapFromGwaspolyGenotype (genotypeFile, plinkFile)
	plinkFile     = createPlinkPedFromGwaspolyGenotype (genotypeFile, plinkFile)
	#plinkFile     = createPlinkPedFromGwaspolyGenotype (genotypeFile, markersIdsMap)
	return (plinkFile)

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

#------------------------------------------------------------------------------
## Create plink MAP file from gwaspoly genotype 
#------------------------------------------------------------------------------
createPlinkMapFromGwaspolyGenotype <- function (gwaspGenoFile, plinkFile) 
{
	#msgmsg ("    >>>> Creating plink MAP file from ", gwaspGenoFile)
	genotype    <- read.table (file=gwaspGenoFile, header=T,stringsAsFactors=T,sep=",", check.names=F)
	map <- genotype [,-(1:3)]
	markers     <- as.character(genotype [,1])
	chromosomes <- genotype [,2]
	positions   <- genotype [,3]

	plinkMap     <- data.frame (chr=chromosomes, iid=markers, dist=0, positions=positions, check.names=F)

	#plinkMapSorted <- plinkMap %>% arrange (chr, positions)
	#write.table (file=outFile, plinkMapSorted, col.names=F, row.names=F, quote=F, sep="\t")
	write.table (file=paste0(plinkFile,".map"), plinkMap, col.names=F, row.names=F, quote=F, sep="\t")
	return (plinkMap$iid)
}


#----------------------------------------------------------
# Create plink PED file from gwaspoly genotype 
#----------------------------------------------------------
createPlinkPedFromGwaspolyGenotype <- function (gwaspGenoFile, plinkFile) 
{
	if (file.exists (paste0(plinkFile,".ped"))) {
		msgmsg ("    >>>> Loading plink file...")
		return (plinkFile)
	}else {
		# ("    >>>> Creating plink PED file...")
		genotype   <- read.csv (file=gwaspGenoFile, header=T, check.names=F)
		alleles    <- as.matrix (genotype [,-c(1,2,3)])
		rownames (alleles) <- genotype [,1]

		# ("    >>>> Creating transposed genotype...")
		samplesIds        <- colnames (alleles)

		# ("    >>>> Getting Ref/Alt Alleles...")
		refAltAlleles <- apply(alleles,1,get.ref)

		# ("    >>>> Converting tetra to diplo")
		allelesDiplo  <- tetraToDiplos (genotype[,-c(2,3)], refAltAlleles)

		# Adjust for plink PED file
		allelesPlink <- t(allelesDiplo)
		genoPED    <- cbind (samplesIds, samplesIds, 0,0,0,-9, allelesPlink)

		# ("    >>>> Writing plink diplo PED file to ", plinkFile)
		write.table (file= paste0 (plinkFile, ".ped"), genoPED, col.names=F, row.names=F, quote=F, sep="\t")
	}
	
	return (plinkFile)
}

#----------------------------------------------------------
#----------------------------------------------------------
gwaspTetraToDiploGenotype <- function (genotypeFile) 
{
	genotype = read.csv (genotypeFile, header=T,check.names=F)
	map      = genotype [,1:3] 
	alleles  = as.matrix (genotype [,-c(1,2,3)])
	rownames (alleles) <- genotype [,1]

	markersIds        <- genotype [,1] 
	samplesIds        <- colnames (alleles)

	#msgmsg ("    >>>> Getting Ref/Alt Alleles...")
	refAltAlleles <- apply(alleles,1,get.ref)

	##msgmsg ("    >>>> Converting tetra to diplo")
	diplosMat  <- tetraToDiplos (genotype[,-c(2,3)], refAltAlleles)
	rownames (diplosMat) = markersIds
	colnames (diplosMat) = samplesIds

	#msgmsg ("    >>>> Writing diplo genotype...")
	genotypeDiplo = data.frame (map, diplosMat, check.names=F)
	outName = addLabel (genotypeFile, "diplo")
	write.csv (file=outName, genotypeDiplo, row.names=F, quote=F)
}

#-------------------------------------------------------------
# Impute NA alleles
#-------------------------------------------------------------
impute.mode <- function(x) {
	ix <- which(is.na(x))
	if (length(ix)>0) {
		x[ix] <- as.integer(names(which.max(table(x))))
	}
	return(x)
}
#----------------------------------------------------------
# Transform table genotype to SHEsis genotype format
#----------------------------------------------------------
gwaspToShesisGenoPheno <- function (genotypeFile, phenotypeFile) 
{
	msgmsg ("    >>>> Writting gwasp to shesis genopheno...")
	sep <- function (allele) {
		s="";
		for (i in 1:4) s=paste0(s, substr (allele,start=i,stop=i)," ");
		return (s)
	}

	geno    <<- read.csv (file=genotypeFile, stringsAsFactors=F, check.names=F)
	pheno   <<- read.csv (file=phenotypeFile, stringsAsFactors=F, check.names=F)
	#stop("LG1")
	rownames (pheno) <- pheno [,1]
	map        <- geno  [,c(1,2,3)]    # Get SNP, Cromosome, Position
	rownames (geno)  <- map [,1] 

	alleles    <- geno[,-c(1,2,3)]
	alleles [is.na(alleles)] = "0000"
	allelesMat <- t(sapply (alleles, sep))

	samples         = rownames (allelesMat)
	pheno           = pheno [samples,]
	pheno [,2]      = impute.mode (pheno [,2])
	genoPhenoShesis = data.frame (Sample=pheno[,1], Trait=pheno[,2],  allelesMat)

	msgmsg ("    >>>> Writing shesis genopheno...")
	outFile = "out/filtered-shesis-genopheno.tbl"
	write.table (file=outFile, genoPhenoShesis, quote=F,row.names=F,col.names=F, sep="\t")

	msgmsg ("    >>>> Writing shesis marker names...")
	outFile = "out/filtered-shesis-markernames.tbl"
	write.table (file=outFile, map[,1], quote=F,row.names=F,col.names=F, sep="\t")
	outFile = "out/filtered-shesis-markernamespos.tbl"
	write.table (file=outFile, map, quote=F,row.names=F,col.names=F, sep="\t")
}

#----------------------------------------------------------
# Add tabs to alleels changign AG --> A	G
#----------------------------------------------------------
tetraToDiplos <- function (allelesMat, refAltAlleles) 
{
	alls <- allelesMat
	if (file.exists ("tmp-diplosMatrix.tbl")) {
		msgmsg ("    >>>> Loading diplos matrix...")
		diplosMat = as.matrix (read.table ("tmp-diplosMatrix.tbl", check.names=F, check.names=F))
	}
	else {
		#msgmsg ("    >>>> Calculating diplos matrix...")
		refs <- refAltAlleles [1,]
		alts <- refAltAlleles [2,]
		
		setB <- function (geno, refs, alts) {
			id      = geno [1]
			alleles = geno [-1]
			ref     = refs [id]
			alt     = alts [id]
			alleles [alleles==strrep (ref,4)] = paste0(ref,ref)
			alleles [alleles==strrep (alt,4)] = paste0(alt,alt)
			alleles [grepl(ref, alleles) & grepl (alt, alleles)] = paste0(alt,ref)
			return (alleles)
		}

		diplosMat  <- t(apply (allelesMat, 1, setB, refs,alts))
		rownames (diplosMat) = allelesMat [,1]
	}
	return (diplosMat)
}
#----------------------------------------------------------
# Convert gwaspoly genotype from ACGT to numeric format
#----------------------------------------------------------
ACGTToNumericGenotypeFormat <- function (genotypeFile) 
{
	geno = read.csv (file=genotypeFile, header=T, check.names=F)
	map <- data.frame(Marker=geno[,1],Chrom=factor(geno[,2],ordered=T),Position=geno[,3],stringsAsFactors=F)
	ploidy=4

	markers <<- as.matrix(geno[,-(1:3)])
	sampleNames <<- colnames (geno[,-(1:3)])
	rownames(markers) <- geno[,1]
	tmp <- apply(markers,1,get.ref)
	map$Ref <- tmp[1,]
	map$Alt <- tmp[2,]

	M <- apply(cbind(map$Ref,markers),1,function(x){
		y <- gregexpr(pattern=x[1],text=x[-1],fixed=T)  
		ans <- as.integer(lapply(y,function(z){ifelse(z[1]<0,ploidy,ploidy-length(z))}))	
		return(ans)
		})

	tM =  (t(M))
	colnames (tM) = sampleNames

	newGeno = data.frame (map[,1:3], tM)
	newName = paste0 (strsplit (genotypeFile, split="[.]")[[1]][1], "-NUM.tbl")
	write.csv (file=newName, newGeno, quote=F, row.names=F)

	return (newName)
}

#----------------------------------------------------------
# Convert numeric (gwaspoly) genotype to ACGT using solcap ref/alt alleles
# Basic genotype: [Markers+Alleles]
# Warning!!! It takes too long
#----------------------------------------------------------
numericToACGTFormatGenotype <- function (genotypeFile, SNPsFile) 
{
	genotype     <- read.csv (genotypeFile, header=T, check.names=F)
	SNPs         <- read.table (SNPsFile, header=T, check.names=F)
	rownames (SNPs) <- SNPs [,1]
	alleles      <- genotype [,-c(2,3)]
	allelesACGT  <- numericToACGTFormatAlleles (alleles, SNPs)
	genoACGT     = data.frame (genotype [,1:3], allelesACGT [,-1])

	outFile = paste0 (strsplit (genotypeFile, split="[.]")[[1]][1],"-ACGT.tbl")
	#msgmsg ("Writing ACGT genotype to ", outFile, "...")
	write.table (file=outFile, genoACGT, quote=F,row.names=F, sep=",")
}

numericToACGTFormatAlleles <- function (alleles, SNPs) 
{
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

##-------------------------------------------------------------
# Convert gwaspoly genotye from numeric tetra to numeric diplo
#-------------------------------------------------------------
numericTetraToNumericDiploGenotype <- function (genotypeFile) {
	toDiplo <- function (markers) {
		id  = markers [1]
		alleles <- as.numeric (markers [-1])
		alleles [alleles==0] = 0
		alleles [alleles==1] = 1
		alleles [alleles==2] = 1
		alleles [alleles==3] = 1
		alleles [alleles==4] = 2
		return (alleles)
	}

	genotype   <- read.csv (genotypeFile, header=T, check.names=F)
	map      = genotype [,1:3]
	alleles  = genotype [,-c(2,3)]

	allelesNum <- t (apply (alleles, 1, toDiplo))
	colnames (allelesNum) = colnames (alleles [-1])
	rownames (allelesNum) = rownames (alleles)

	newGeno = cbind (map, allelesNum)
	newName = addLabel (genotypeFile, "diploNUM")
	write.csv (file=newName, newGeno, row.names=F, quote=F)
}


numericToABGenotypeFormat <- function (genotypeFile) 
{
	geno = read.csv (file=genotypeFile, header=T, check.names=F)
	markers <<- geno [,-(2:3)]
	rownames (markers) = geno [,1]

	toAB <- function (markers) {
		id  = markers [1]
		alleles <- as.numeric (markers [-1])
		ref = "A"
		alt = "B"
		alleles [alleles==0] = strrep(ref,4)
		alleles [alleles==1] = paste0 (strrep(alt,1),strrep(ref,3))
		alleles [alleles==2] = paste0 (strrep(alt,2),strrep(ref,2))
		alleles [alleles==3] = paste0 (strrep(alt,3),strrep(ref,1))
		alleles [alleles==4] = strrep(alt,4)
		return (alleles)
	}

	allelesAB <- t(apply (markers, 1, toAB))
	colnames (allelesAB) = colnames (markers [-1])
	rownames (allelesAB) = rownames (markers)

	newGeno = data.frame (geno [,c(1:3)], allelesAB)
	newName = paste0 (strsplit (genotypeFile, split="[.]")[[1]][1], "-AB.tbl")
	write.csv (file=newName, newGeno, quote=F, row.names=F)
}

#-------------------------------------------------------------
# Add label to filename
#-------------------------------------------------------------
addLabel <- function (filename, label)  {
	nameext = strsplit (filename, split="[.]")
	newName = paste0 (nameext [[1]][1], "-", label, ".", nameext [[1]][2])
	return (newName)
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
# Run a command string using system function and writes output to log file
#-------------------------------------------------------------
runCommand <- function (command, logFile="gwas.log", DEBUG=F) 
{
	if (DEBUG==T) {
		msgmsg (">>>> ", command)
		system (command)
	}else {
		#msgmsg (">>>> ", command)
		errorsLog = paste0 (strsplit(logFile, split="[.]")[[1]], ".errors")
		system (paste0 (command, " > ", logFile," 2> ",errorsLog))
	}
}

#-------------------------------------------------------------
#-------------------------------------------------------------
runCommand <- function (command, logFile="") {
	msgmsg (">>>> ", command)
	if (logFile != "")
		system (paste0 (command, " > ", logFile))
	else
		system (command)
}

hd <- function (data, m=10,n=10) {
	msgmsg (deparse (substitute (data)),":")
	if (is.null (dim (data)))
		print (data [1:10])
	else if (ncol (data) < 5) 
		print (data[1:m,])
	else if (nrow (data) < 10)
		print (data[,1:n])
	else 
		print (data [1:m, 1:n])
}

#----------------------------------------------------------
# Main
#----------------------------------------------------------
main <- function () 
{
	args = c ("filtered-gwasp4-phenotype.tbl", "filtered-plink-genotype.ped")

	pheno  = args [1]
	geno   = args [2]

	gwaspToTasselPhenotype (pheno, "filt-tassel-pheno.tbl")
	plinkToVCFFormat (geno, "filt-tassel-geno")

}
#----------------------------------------------------------
#----------------------------------------------------------
#main ()


