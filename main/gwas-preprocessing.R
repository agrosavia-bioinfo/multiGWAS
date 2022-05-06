#!/usr/bin/Rscript

# INFO  : Different functions to convert between different genotype/phenotype formats
# AUTHOR: Luis Garreta (lgarreta@agrosavia.co)
# DATE  : 12/Feb/2020
# LOG   : 
#         r1.31: Added GWASpoly To GAPIT numeric

#options (width=300, stringsAsFactors=F)
#args = commandArgs(trailingOnly=T)

HOME = Sys.getenv ("MULTIGWAS_HOME")
.libPaths (paste0(HOME, "/opt/Rlibs"))

suppressMessages (library (parallel))
suppressMessages (library (dplyr))
suppressMessages (library (stringi))
formatsLogFile="log-formats.log"
#----------------------------------------------------------
# Main
#----------------------------------------------------------
#msgmsg = message
main <- function () 
{
	DEBUG <<- F
	options (width=300)
	args = commandArgs (trailingOnly=T)

	genotypeFile = args [1]
	ACGTToNumericGenotypeFormat (genotypeFile, 4)

	#convertVCFToACGTByNGSEP (genotypeFile) 

	#convertVCFToGenodiveFormat (genotypeFile)
	#numericTetraToNumericDiploGSMatrix (genotypeFile)
	#convertAABBGWASpolyToNumericFormat (inputFilename)
	#gwaspTetraToDiploGenotype (genotypeFile)
	#getCommonGenoPhenoMap (args[1], args[2])
	#createGwaspolyGenotype (args [1], args [2], 4)
	#convertFitpolyToKMatrix (args [1])
	#convertFitpolyToGwaspolyGenotype (args [1], args [2])

	#convertUpdogToKMatrix (args [1])
	#convertUpdogToGwaspolyGenotype (args [1], args [2])

	## Missing filter
	#args = c ("example-genotype-tetra-gwaspoly-ACGT.csv", "example-phenotype-single-trait.csv")
	#genotypeFile = args [1]
	#geno = filterByMissingMarkersAndSamples (genotypeFile, 0.3, 0.2)

}


#-------------------------------------------------------------
# Get common sample names
#-------------------------------------------------------------
filterByCommonMarkersSamples <- function (genotypeFile, phenotypeFile) {
	msgmsg ("Filtering by missing markers and samples...")
	geno     = read.csv (file=genotypeFile, check.names=F)
	pheno    = read.csv (file=phenotypeFile, check.names=F)

	# From phenotype, remove duplicated and NA samples
	pheno  = pheno [!is.na(pheno[,2]),]
	pheno  = pheno [!duplicated (pheno[,1]),]

	# From genotype, remove duplicated markers
	geno  = geno [!duplicated (geno[,1]),]

	# Get common samples
	genoSamples   = colnames (geno[,-(1:3)])
	phenoSamples  = pheno [,1] 
	commonSamples = intersect (genoSamples, phenoSamples) 
	commonMarkers = as.character(geno [,1])

	# Create new geno, pheno
	genoCommon     = geno  [,c (colnames(geno)[1:3], commonSamples)]
	phenoCommon    = pheno [pheno[,1] %in% commonSamples,]
	mapCommon      = genoCommon [, 1:3]
	trait          = colnames(phenoCommon)[2]

	genoCommonFile   = addLabel (genotypeFile, "COMMON")
	phenoCommonFile  = addLabel (phenotypeFile, "COMMON")
	mapCommonFile    = sprintf ("%s/%s", dirname (genoCommonFile), "filtered-map-COMMON.csv")

	write.csv (file=genoCommonFile, genoCommon, row.names=F)
	write.csv (file=phenoCommonFile, phenoCommon, row.names=F)
	write.csv (mapCommon, mapCommonFile, row.names=F)

	return (list (genotypeFile=genoCommonFile, phenotypeFile=phenoCommonFile, mapFile=mapCommonFile, trait=trait))
}



#-------------------------------------------------------------
# Check call rate of samples and markers
#-------------------------------------------------------------
filterByMissingMarkersAndSamples <- function (genotypeFile, callRateSNPs, callRateSamples) 
{	
	geno          = read.csv (file=genotypeFile, check.names=F)
	allelesMatrix = geno [,-c(1:3)]
	SNPs          = colnames (allelesMatrix)
	samples       = colnames (allelesMatrix)
	markers       = geno [,1]
	rownames (geno) = markers
	rownames (allelesMatrix) = markers

	#------- local function -------
	calculateProportionNAs <- function (alleles) {
		prop = sum (is.na (alleles)) / length (alleles)
		return (prop)
	}
	#------------------------------
	propsSNPs    = apply (allelesMatrix, 1, calculateProportionNAs) # By row (Marker)
	propsSamples = apply (allelesMatrix, 2, calculateProportionNAs) # By column (Samples)
	badSNPs      = which (propsSNPs > callRateSNPs)
	badSamples   = which (propsSamples > callRateSamples)

	if (length (badSNPs) > 0) {
		markers = setdiff (markers, names (badSNPs))
		outFile = addLabel (genotypeFile, "MISSING-SNPs")
		write.csv (names (badSNPs), outFile, quote=F, row.names=F)
	}

	if (length (badSamples) > 0) {
		samples = setdiff (samples, names (badSamples))
		outFile = addLabel (genotypeFile, "MISSING-SAMPLES")
		write.csv (names (badSamples), outFile, quote=F, row.names=F)
	}

	noMissingGeno      = cbind (geno [markers, 1:3], geno [markers, samples])
	noMissingGenoFile  = addLabel (genotypeFile, "noMISSING")
	write.csv (noMissingGeno, noMissingGenoFile, quote=F, row.names=F)

	return (noMissingGenoFile)
}



old_filterByMissingMarkersAndSamples <- function (genotypeFile, callRateSNPs, callRateSamples) 
{	
	geno = read.csv (file=genotypeFile, check.names=F)
	allelesMatrix = geno [,-c(1:3)]
	rownames (allelesMatrix) = geno [,1]
	map  = geno [,1:3]

	#------- local function -------
	calculateProportionNAs <- function (alleles) {
		prop = sum (is.na (alleles)) / length (alleles)
		return (prop)
	}
	#------------------------------
	propsSNPs    = apply (allelesMatrix, 1, calculateProportionNAs) # By row (Marker)
	propsSamples = apply (allelesMatrix, 2, calculateProportionNAs) # By column (Samples)
	badSNPs      = which (propsSNPs > callRateSNPs)
	badSamples   = which (propsSamples > callRateSamples)

	if (length (badSNPs) > 0 && length (badSamples) > 0) {
		noMissingGeno = geno [-badSNPs, -badSamples]
		missingGeno   = geno [badSNPs, badSamples]

		noMissingGenoFile  = addLabel (genotypeFile, "noMISSING")
		write.csv (noMissingGeno, noMissingGenoFile, quote=F, row.names=F)
		genotypeFile = noMissingGenoFile

		outFile = addLabel (genotypeFile, "MISSING-SNPs")
		write.csv (names (badSNPs), outFile, quote=F, row.names=F)

		outFile = addLabel (genotypeFile, "MISSING-SAMPLES")
		write.csv (names (badSamples), outFile, quote=F, row.names=F)
	}else  if (length (badSNPs) > 0 && length (badSamples) == 0) {
		noMissingGeno = geno [-badSNPs,]
		missingGeno   = geno [badSNPs,]

		noMissingGenoFile  = addLabel (genotypeFile, "noMISSING")
		write.csv (noMissingGeno, noMissingGenoFile, quote=F, row.names=F)
		genotypeFile = noMissingGenoFile

		outFile = addLabel (genotypeFile, "MISSING-SNPs")
		write.csv (names (badSNPs), outFile, quote=F, row.names=F)
	}else  if (length (badSNPs) == 0 && length (badSamples) > 0) {
		noMissingGeno = geno [,-badSamples]
		noMissingGenoFile  = addLabel (genotypeFile, "noMISSING")
		write.csv (noMissingGeno, noMissingGenoFile, quote=F, row.names=F)
		genotypeFile = noMissingGenoFile

		outFile = addLabel (genotypeFile, "MISSING-SAMPLES")
		write.csv (names (badSamples), outFile, quote=F, row.names=F)
	}

	return (genotypeFile)
}


#-------------------------------------------------------------
# Impute values by row (SNP) from matrix with rownames at first column
#-------------------------------------------------------------
imputeNumericMatrix <- function (numericMatrix) {
	# Impute NA alleles
	impute.mode <- function(x) {
		ix <- which(is.na(x))
		if (length(ix)>0) {
			x[ix] <- as.integer(names(which.max(table(x))))
		}
		return(x)
	}
	matImputed <- t(apply(numericMatrix[,-1], 1, impute.mode))
	matImputed = data.frame (numericMatrix[,1,drop=F], matImputed, check.names=F)
	return (matImputed)
}

#-------------------------------------------------------------
# Convert GWASpoly phenotype fto GAPIT format
#-------------------------------------------------------------
gwaspolyToGapitPhenotype <- function (phenotypeFile) {
	pheno      = read.csv (phenotypeFile, check.names=F)
	phenoGapit = data.frame (taxa=pheno[,1], pheno [,2, drop=F], check.names=F)

	gapitPhenoFile = addLabel (phenotypeFile, "GAPIT")
	write.csv (phenoGapit, gapitPhenoFile, row.names=F)

	return (gapitPhenoFile)
}
#-------------------------------------------------------------
# Convert GWASpoly genotype to GAPIT format: numeric and map matrix
#-------------------------------------------------------------
gwaspolyToGapitGenotype <- function (genotypeFile, geneAction, FILES=F) {
	# Convert genotype to numeric
	num          = ACGTToNumericGenotypeFormat (genotypeFile, 4, MAP=T) 
	genoNumTetra = num$genoNum [,-c(2,3)]
	genoMap      = num$map [,-c(4,5)]

	genoNumTetraImputed  = imputeNumericMatrix (genoNumTetra)
	colnames (genoMap) = c("SNP", "Chromosome", "Position")

	# Convert to diplo
	genoNumDiplo = numericTetraToNumericDiploMatrix (genoNumTetraImputed)
	genoNumGapit = t (genoNumDiplo[,-1])
	genoNumGapit = data.frame (taxa=colnames (genoNumDiplo[,-1]), genoNumGapit, check.names=F)

	# Convert genotype for additive or dominant geneAction
	# Under the dominant model, both types of homozygous genotypes 
	# are coded 0 and heterozygous genotypes are coded 1.
	if (geneAction=="additive") {
		gapitNumFile = addLabel (genotypeFile, "GAPIT-NUM-ADD")
		write.csv (genoNumGapit, gapitNumFile, quote=F, row.names=F)
	}else if (geneAction=="dominant") {
		# 0,2 -> 0, 1->1
		genoNumGapitValues = genoNumGapit [,-1] 
		genoNumGapitValues [genoNumGapitValues==2] = 0
		#genoNumGapitValues [genoNumGapitValues==1] = 1
		#genoNumGapitValues [genoNumGapitValues==2] = 0
		genoNumGapit = data.frame (genoNumGapit [,1,drop=F], genoNumGapitValues, check.names=F)
		gapitNumFile = addLabel (genotypeFile, "GAPIT-NUM-DOM")
		write.csv (genoNumGapit, gapitNumFile, quote=F, row.names=F)
	}else if (geneAction=="all") {
		genoNumGapitValues = genoNumGapit [,-1] 

		# 0,2 -> 0, 1->1
		genoNumGapitValues [genoNumGapitValues==2] = 0
		genoNumGapit = data.frame (genoNumGapit [,1,drop=F], genoNumGapitValues, check.names=F)
		gapitNumFile = addLabel (genotypeFile, "GAPIT-NUM-HET")
		write.csv (genoNumGapit, gapitNumFile, quote=F, row.names=F)

		# 0,2 -> 1, 1->0
		genoNumGapitValues [genoNumGapitValues==0] = 2
		genoNumGapitValues [genoNumGapitValues==1] = 0
		genoNumGapitValues [genoNumGapitValues==2] = 1
		genoNumGapit = data.frame (genoNumGapit [,1,drop=F], genoNumGapitValues, check.names=F)
		gapitNumFile = addLabel (genotypeFile, "GAPIT-NUM-HOM")
		write.csv (genoNumGapit, gapitNumFile, quote=F, row.names=F)

		# Dominant: AA (0) and Aa(1) to 1, aa (2) to 0
		# 0,1 -> 0, 2 -> 1
		genoNumGapitValues [genoNumGapitValues==2] = 1
		genoNumGapitValues [genoNumGapitValues==1] = 0
		genoNumGapit = data.frame (genoNumGapit [,1,drop=F], genoNumGapitValues, check.names=F)
		gapitNumFile = addLabel (genotypeFile, "GAPIT-NUM-DOM")
		write.csv (genoNumGapit, gapitNumFile, quote=F, row.names=F)

		# Recessive: AA (0) and Aa(1) to 0, aa (2) to 1
		# 0 -> 0, 1,2 -> 1
		genoNumGapitValues [genoNumGapitValues==2] = 1
		genoNumGapit = data.frame (genoNumGapit [,1,drop=F], genoNumGapitValues, check.names=F)
		gapitNumFile = addLabel (genotypeFile, "GAPIT-NUM-REC")
		write.csv (genoNumGapit, gapitNumFile, quote=F, row.names=F)
	}

	gapitMapFile = addLabel (genotypeFile, "GAPIT-MAP")
	write.csv (genoMap, gapitMapFile, quote=F, row.names=F)
	return (list (geno=gapitNumFile, map=gapitMapFile))
}

#-------------------------------------------------------------
# Get common sample names
#-------------------------------------------------------------
getCommonGenoPhenoMapFile <- function (genotypeFile, phenotypeFile, mapFile=NULL) 
{
	geno  = read.csv (file=genotypeFile, check.names=F)
	pheno = read.csv (file=phenotypeFile, check.names=F)
	
	common = commonGenoPheno (geno, pheno)

	genoCommonFile  = addLabel (genotypeFile, "COMMON")
	phenoCommonFile = addLabel (phenotypeFile, "COMMON")

	write.csv (file=genoCommonFile, common$geno, quote=F, row.names=F)
	write.csv (file=phenoCommonFile, common$pheno, quote=F, row.names=F)

	return (list (geno=genoCommon, pheno=phenoCommon))
}

commonGenoPheno <- function (geno, pheno) {
	# From phenotype, remove duplicated and NA samples
	pheno  = pheno [!is.na(pheno[,2]),]
	pheno  = pheno [!duplicated (pheno[,1]),]

	# From genotype, remove duplicated markers
	geno  = geno [!duplicated (geno[,1]),]

	# Get common samples
	genoSamples   = colnames (geno[,-(1:3)])
	phenoSamples  = pheno [,1] 
	commonSamples = intersect (genoSamples, phenoSamples) 

	# Create new geno, pheno
	genoCommon  = geno  [,c (colnames(geno)[1:3], commonSamples)]
	phenoCommon = pheno [pheno[,1] %in% commonSamples,]

	return (list(geno=genoCommon, pheno=phenoCommon))
}
#----------------------------------------------------------
# Get reference/alternate alleles from GWASpoly genotype
#----------------------------------------------------------
getReferenceAlternateAlleles <- function (genotypeFile, ploidy) 
{
	geno = read.csv (file=genotypeFile, check.names=F)
	map  = geno [,2:3]

	markers           = as.matrix(geno[,-(1:3)])
	sampleNames       = colnames (geno[,-(1:3)])
	rownames(markers) = geno[,1]
			
	tmp <- apply(markers,1,getReferenceAllele)
	ref <- tmp[1,]
	alt <- tmp[2,]

	refAltTable = data.frame (Markers=geno[,1], Ref=ref, Alt=alt, map)
	write.csv (refAltTable, addLabel (genotypeFile, "REF-ALT"), quote=F, row.names=F)
}

#-------------------------------------------------------------
# Convert GWASpoly genotype to fitPoly scores
# Write two table: fitPoly scores and map info
#-------------------------------------------------------------
convertGwaspolyToFitpolyScores <- function (gwpFile) 
{
	#---- Local function ------
	convert <- function (row) {
		sampleNames = names (row)[-c(1:2)]
		values      = as.numeric (row [-c(1:2)])
		n           = length (values)
		maxgeno     = geno = values
		geno [is.na(values)] = ""
		rowDF = data.frame (
			marker=as.character (row[1]), MarkerName=as.character (row[2]),
			SampleName=sampleNames, ratio=runif (n,0,1), P0=runif (n,0,1), 
			P1=runif (n,0,1), P2=runif (n,0,1), P3=runif (n,0,1), P4=runif (n,0,1),
			maxgeno=values, maxP=runif (n,0,1), geno=geno)
		return (rowDF)
	}

	gwpTable = read.csv (gwpFile) 
	gwpTable = gwpTable [order(gwpTable[,1]),]
	mapTable = gwpTable [,1:3]
	write.csv (mapTable, addLabel (gwpFile, "MAP"), quote=F, row.names=F)

	gwpTable = gwpTable [,-2:-3]
	gwpTable = data.frame (marker=seq(nrow(gwpTable)), gwpTable)

	fitpolyTable = do.call (rbind.data.frame, apply (gwpTable, 1, convert))
	write.table (fitpolyTable, addLabel (gwpFile, "FITPOLY"), quote=F, row.names=F, sep="\t")
}
#------------------------------------------------------------------------------
# Convert VCF to ACGT files using "NGSEP"
#------------------------------------------------------------------------------
convertVCFToACGTByNGSEP <- function (filename, outFilename="") 
{
	stemName = strsplit (filename, "[.]")[[1]][1]
	cmm=sprintf ("java -jar %s/opt/tools/MultiGWAS_NGSEP.jar VCFConverter -GWASPoly -i %s -o %s", HOME, filename, stemName)
	runCommand (cmm, "log-NGSEP.log")
	outFilename = paste0 (stemName, "_GWASPoly.csv") # Added by NGSEP tool
	geno = read.csv (outFilename, check.names=F)
	geno = arrange (geno, Chrom, Position)
	write.csv (geno, outFilename, quote=F, row.names=F)
	return (outFilename)
}

#-------------------------------------------------------------
# Different routines to check valid data for the different tools
#-------------------------------------------------------------
checkForValidPlinkPhenotype <- function (phenotypeFile) 
{
	phenotypeValues = read.csv (phenotypeFile)[,2]
	levelValues = levels (factor(phenotypeValues))
	if (any (levelValues > 2)) 
		return (list(value=TRUE,info="")) 
	else
		return (list(value=FALSE, info=
"No valid phenotype for PLINK. If phenotype is quantitative, 
it contains values other than 1, 2, 0 or missing")) 
	
}

#-------------------------------------------------------------
# Crete GWASpoly ACGT genotype from numeric kmatrix genotype 
# and map file (Marker, Ref, Alt, Chrom, Pos)
#-------------------------------------------------------------
createGwaspolyGenotype <- function (genotypeFile, mapFile) 
{
	geno  = read.csv (genotypeFile, check.names=F, row.names=1)
	map   = read.csv (mapFile, check.names=F, row.names=1)

	commonMarkers = intersect (rownames(geno),rownames(map))
	genoCommon    = geno [commonMarkers,]
	mapCommon     = map [commonMarkers,]

	gwaspolyGeno  = data.frame (Marker=rownames (genoCommon), 
								Chromosome=map [,3],
								Position=map [,4],
								genoCommon)

	outFile = addLabel (genotypeFile, "GWASPOLY")
	write.csv (gwaspolyGeno, outFile, quote=F, row.names=F)
	outFile = numericToACGTFormatGenotype (outFile, mapFile)
	return (outFile)
}

#----------------------------------------------------------
# Convert updog file to genotype matrix
#----------------------------------------------------------
convertUpdogToGwaspolyGenotype <- function (updogScoresFile, mapFile) 
{
	kmatrixFile  = convertUpdogToKMatrix (updogScoresFile)
	gwaspolyFile = createGwaspolyGenotype (kmatrixFile, mapFile)
	return (gwaspolyFile)
}
#----------------------------------------------------------
# Convert updog file to genotype matrix
#----------------------------------------------------------
convertUpdogToKMatrix <- function (updogScoresFile) 
{
	#--- Function to convert fitPoly table to GWASpoly matrix
	createMatrix <- function (snp, fitGenos) {
		dataSnp = fitGenos [fitGenos$snp==snp,]
		genos   = data.frame (t(dataSnp$geno))
		if (all (is.na (genos))) 
			return (NULL)
		names (genos)    = dataSnp$ind
		rownames (genos) = snp
		return (genos)
	 }
	#------------------------------------------------------

	fitGenos      = read.csv (updogScoresFile)		
	snpList       = levels (fitGenos$snp)
	nCores        = ifelse (detectCores()>1, detectCores()-1, 1)
	nCores = 1

	outs          = mclapply (snpList, createMatrix, fitGenos, mc.cores=nCores)
	updogGenosDF  = do.call (rbind.data.frame, outs)

	outFilename   = addLabel (updogScoresFile, "KMATRIX")

	updogMatrix = data.frame (Markers=rownames (updogGenosDF), updogGenosDF)
	write.csv (updogMatrix, outFilename, quote=F, row.names=F)
	return (outFilename)
}

#----------------------------------------------------------
# Convert fitpoly file to genotype matrix
#----------------------------------------------------------
convertFitpolyToGwaspolyGenotype <- function (fitpolyScoresFile, mapFile) 
{
	kmatrixFile  = convertFitpolyToKMatrix (fitpolyScoresFile)
	gwaspolyFile = createGwaspolyGenotype (kmatrixFile, mapFile)
	return (gwaspolyFile)
}
#----------------------------------------------------------
# Convert fitpoly file to genotype matrix
#----------------------------------------------------------
convertFitpolyToKMatrix <- function (fitpolyScoresFile) 
{
	#--- Function to convert fitPoly table to GWASpoly matrix
	createMatrix <- function (snp, fitGenos) {
		dataSnp = fitGenos [fitGenos$MarkerName==snp,]
		genos   = data.frame (t(dataSnp$geno))
		if (all (is.na (genos))) 
			return (NULL)
		names (genos)    = dataSnp$SampleName
		rownames (genos) = snp
		return (genos)
	 }
	#------------------------------------------------------

	fitGenos      = read.table (fitpolyScoresFile, sep="\t", header=T)		
	snpList       = levels (fitGenos$MarkerName)
	nCores        = ifelse (detectCores()>1, detectCores()-1, 1)

	outs          = mclapply (snpList, createMatrix, fitGenos, mc.cores=nCores)
	fitGenosDF    = do.call (rbind.data.frame, outs)

	outFilename   = addLabel (fitpolyScoresFile, "KMATRIX")

	fitPolyMatrix = data.frame (Markers=rownames (fitGenosDF), fitGenosDF)
	write.csv (fitPolyMatrix, outFilename, quote=F, row.names=F)
	return (outFilename)
}

#-------------------------------------------------------------
# Return the format type of genotype
# Checks if VCF, GWASpoly, or k-matrix.
#-------------------------------------------------------------
checkGenotypeFormat <- function (genotypeFile, ploidy) {
	msg ("Checking genotype file format...")
	con = file (genotypeFile, "r")
	firstLine = readLines (con, n=1)
	close (con)

	# Check if VCF
	if (grepl ("VCF", firstLine)) {
		genotypeFile = convertVCFToACGTByNGSEP (genotypeFile) #output: filename.csv
		return (genotypeFile)
	}

	# Compare if all cells have the same lengths (nchars)
	# if True --> matrix type, False --> GWASpoly type
	equalLength <- function (cell, ploidy) {
		cellChar = as.character (cell)
		return (nchar (cellChar)==ploidy)
	}	

	# Check if k-matrix
	data   = read.csv (genotypeFile, row.names=1, check.names=F)
	sample = data [1:10,1:10]

	if (all(sapply(sample, equalLength, ploidy), na.rm=T)==TRUE) {
		N = nrow (data)
		newFilename = paste0 (strsplit (genotypeFile, "[.]")[[1]][1], "_GWASPoly.csv")
		genoGwaspoly = data.frame (Markers=rownames (data), Chrom=1, Position=1:N, data)
		write.csv (genoGwaspoly, newFilename, row.names=F)
		return (newFilename)
	}
		
	# Check if GWASpoly
	sample = data [1:10,3:10]
	if (all(sapply(sample, equalLength, ploidy), na.rm=T)==TRUE)
		return (genotypeFile)

	return ("Unknown-genotype")
}

#------------------------------------------------------------------------------
# Convert VCF to ACGT files
#------------------------------------------------------------------------------
convertVCFtoACGTByVCFR <- function (filename, outFilename="") {
	suppressMessages (library (vcfR))
	vcf = read.vcfR (filename, verbose=F)

	# Extract genotipe strings" 
	gt = extract.gt (vcf, return.alleles=T)

	# Extract map info (id, chrom, pos) Merge gt + info
	map = vcf@fix [,c(1:2,4:5)]
	gtmap = cbind (MARKERS=rownames (gt), map, gt)

	# >>>>>>>> Convert VCF Genotypes to ACGT. First by row, then by cell
	vcfToACGTForCell <- function (allelesCell) {
		if (is.na (allelesCell))
			return (NA)

		numsall    = strsplit (allelesCell, split="[|/]")
		str = sapply (numsall[[1]], substring, 1, 1)
		acgt = stri_reverse (paste (str, collapse=""))
		return (acgt)
	}
	vcfToACGTForRow <- function (allelesRow) {
		rows = sapply (allelesRow, vcfToACGTForCell)
		return (rows)
	} # ">>>>>>>>>>>"

	allelesMat = gtmap [,-1:-5]
	allelesACGT <- t(apply (allelesMat, 1, vcfToACGTForRow))
	colnames (allelesACGT) = colnames (allelesMat)
	rownames (allelesACGT) = rownames (allelesMat)

	newAlleles <- data.frame (gtmap[,1:3], allelesACGT)
	if (outFilename=="")
		outFilename = paste0 (strsplit (filename, "[.]")[[1]][1], ".csv")

	write.csv (newAlleles, outFilename, row.names=F, quote=F)
	return (outFilename)
}
#------------------------------------------------------------------------------
# Convert genotye from plink (.ped, .map) to VCF (Variant Call Format)
#------------------------------------------------------------------------------

plinkToVCFFormat <- function (plinkFile) {
	#plinkPrefix = strsplit (plinkFile, split="[.]")[[1]][1]
	tasselFile = gsub ("PLINK", "TASSEL", plinkFile)
	cmm = sprintf ("plink --file %s --allow-extra-chr --recode vcf-fid --out %s", plinkFile, tasselFile)
	runCommand (cmm)
	outFile = paste0 (tasselFile, ".vcf")

	return (outFile)
}

#------------------------------------------------------------------------------
## Format and write gwasp to tassel phenotype (For rtassel)
#------------------------------------------------------------------------------
gwaspolyToTasselPhenotype<- function (gwaspPhenotypeFile) 
{
	gwaspPhenotype = read.csv (file=gwaspPhenotypeFile, header=T, check.names=F, )
	taxa           = as.character (gwaspPhenotype [,1])
	trait          = gwaspPhenotype [,2]
	traitName      = colnames (gwaspPhenotype)[2]
	tasselPheno    = cbind (taxa, trait)

	colnames (tasselPheno) = c ("<Trait>", traitName)

	outFilename = addLabel (gwaspPhenotypeFile, "TASSEL")
	write.table (file=outFilename, tasselPheno, col.names=T, row.names=F, quote=F, sep="\t")
	
	return (outFilename)
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
gwaspolyToPlinkPhenotype <- function (gwaspPhenoFile) {
	#msgmsg ("    >>>> Creating plink phenotype...")
	phenotype = read.csv (file=gwaspPhenoFile, check.names=F)
	idNames = as.character (phenotype [,1])

	samples = phenotype [,1]
	traits  = phenotype [,2]

	#msgmsg ("    >>>> Writing plink phenotype...")
	plinkPheno = data.frame (FID=samples,IID=samples, TRAIT= traits)
	#if (outFile=="")
	#	outFile = paste0 ("out/",strsplit (gwaspPhenoFile, split="[.]")[[1]][1], "-plink.tbl")

	plinkPhenoFile = addLabel (gwaspPhenoFile, "PLINK")
	write.table (plinkPheno, plinkPhenoFile, col.names=T, row.names=F, quote=F, sep="\t")
	return (plinkPhenoFile)
}

#------------------------------------------------------------------------------
# Gwasp to plink format (.ped .map)
#------------------------------------------------------------------------------
gwaspolyToPlinkGenotype <- function (genotypeFile) {
	plinkFile     = addLabel (genotypeFile, "PLINK")
	plinkFile     = strsplit (plinkFile, "[.]")[[1]][1]
	markersIdsMap = createPlinkMapFromGwaspolyGenotype (genotypeFile, plinkFile)
	plinkFile     = createPlinkPedFromGwaspolyGenotype (genotypeFile, plinkFile)
	cmm = paste ("plink --file", plinkFile, "--allow-extra-chr --make-bed", "--out", plinkFile)
	runCommand (cmm, "log-filtering.log")
	#plinkFile     = createPlinkPedFromGwaspolyGenotype (genotypeFile, markersIdsMap)
	return (plinkFile)
}
#----------------------------------------------------------
# Get ref/alt alleles for SNPs
# Replaces the GWASpoly ref.alt that get random ref/alt alleles
#----------------------------------------------------------
bases <- c("A","C","G","T")
getReferenceAllele <- function (x) {
	y   <- paste(na.omit(x),collapse="")
	ans <- apply(array(bases),1,function(z,y){length(grep(z,y,fixed=T))},y)

	if (sum(ans)>2) {stop("Error in genotype matrix: More than 2 alleles")}

	if (sum(ans)==2) {
		ref.alt <- bases [which (ans==1)]
		countAllele1 = stri_count (y, fixed=ref.alt[1])
		countAllele2 = stri_count (y, fixed=ref.alt[2])
		if (countAllele2 > countAllele1)
			ref.alt = c(ref.alt[2],ref.alt[1])
	}
	if (sum(ans)==1) {ref.alt <- c(bases[which(ans==1)],NA)}

	if (!"ref.alt" %in% ls()) {
		message ("------------------------------------------------")
		message ("ERROR: The genotype has errors. Program stopped!")
		message ("------------------------------------------------")
		quit ()
	}
	return (ref.alt)
}


#----------------------------------------------------------
# Wrong: get random get/ref alleles
#----------------------------------------------------------
get.ref <- function(x) {
	y <- paste(na.omit(x),collapse="")
	ans <- apply(array(bases),1,function(z,y){length(grep(z,y,fixed=T))},y)

	if (sum(ans)>2) {stop("Error in genotype matrix: More than 2 alleles")}

	if (sum(ans)==2) {
		ref.alt <- bases[which(ans==1)]
	}
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
	map         <- genotype [,-(1:3)]
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
	# Creating plink PED file
	genotype   <- read.csv (file=gwaspGenoFile, header=T, check.names=F)
	alleles    <- as.matrix (genotype [,-c(1,2,3)])
	rownames (alleles) <- genotype [,1]

	# Creating transposed genotype
	samplesIds        <- colnames (alleles)

	# Getting Ref/Alt Alleles
	#refAltAlleles <- apply(alleles,1,get.ref)
	refAltAlleles <- apply(alleles,1,getReferenceAllele)

	# Converting tetra to diplo
	allelesDiplo  <- tetraToDiplos (genotype[,-c(2,3)], refAltAlleles)

	# Adjust for plink PED file
	allelesPlink <- t(allelesDiplo)
	genoPED    <- cbind (samplesIds, samplesIds, 0,0,0,-9, allelesPlink)

	# Writing plink diplo PED file to plinkFile
	write.table (file= paste0 (plinkFile, ".ped"), genoPED, col.names=F, row.names=F, quote=F, sep="\t")
	
	return (plinkFile)
}

#----------------------------------------------------------
#----------------------------------------------------------
gwaspTetraToDiploGenotype <- function (genotypeFile) {
	#msgmsg ("Converting genotype from tetraploid to diploid...")
	genotype = read.csv (genotypeFile, header=T,check.names=F)
	map      = genotype [,1:3] 
	alleles  = as.matrix (genotype [,-c(1,2,3)])
	rownames (alleles) <- genotype [,1]

	markersIds        <- genotype [,1] 
	samplesIds        <- colnames (alleles)

	#msgmsg ("    >>>> Getting Ref/Alt Alleles...")
	#refAltAlleles <- apply(alleles,1,get.ref)
	refAltAlleles <- apply(alleles,1,getReferenceAllele)

	#msgmsg ("    >>>> Converting tetra to diplo")
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
# Add tabs to alleels changign AG --> A	G
#----------------------------------------------------------
tetraToDiplos <- function (allelesMat, refAltAlleles) {
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
# Convert AABB to 0,1,2,3,4
#----------------------------------------------------------
convertAABBGWASpolyToNumericFormat <- function (genotypeFile) 
{
	geno   = read.csv (file=genotypeFile, check.names=F)[,1:13]

	allelesAB           = as.matrix(geno[,-(1:3)])
	map                 = geno [,1:3]
	sampleNames         = colnames (geno[,-(1:3)])
	markerNames         = geno [1,]
	rownames(allelesAB) = geno[,1]

	setNum <- function (allelesABMarker) {
		allelesNumMarker = sapply (allelesABMarker, stri_count, fixed="B")
	}

	allelesNumList = unlist (mclapply (allelesAB,  setNum))
	allelesMat     = matrix (allelesNumList, nrow=nrow(allelesAB), ncol=ncol(allelesAB))
	colnames (allelesMat) = colnames (allelesAB)

	newGeno = cbind (map, allelesMat)
	newName = gsub (".csv", "-NUM.csv", genotypeFile)
	write.csv (newGeno, newName, quote=F, row.names=F)
}

#-------------------------------------------------------------
# Convert genotype in VCF format to kmatrix numeric format (A:1, C:2, G:3, T:4)
#----------------------------------------------------------
convertVCFToGenodiveFormat <- function (genotypeFile) {
	message ("Converting to kmatrix....")
	kmatrixFile = convertVCFToACGTByNGSEP (genotypeFile)
	message ("Converting to Genodive")

	NCORES = detectCores()
	geno   = read.csv (file=kmatrixFile, header=T, check.names=F)

	markers           = as.matrix(geno[,-(1:3)])
	sampleNames       = colnames (geno[,-(1:3)])
	rownames(markers) = geno[,1]

	#>>>> Convert one row from ACGT to Num. x=RefAllele|...Alleles...
	acgtToNum <- function(ACGTList){
		changeGenotype <- function (ACGT) {
			newGeno = gsub ("A","1",gsub ("C","2",gsub("G","3",gsub("T","4",ACGT))))
		}
		numbersList = sapply (ACGTList, changeGenotype)
	}
	
	# Convert All ACGT matrix to Num
	ACGTList       = mclapply(seq(ncol(markers)), function(i) markers[,i], mc.cores=NCORES)
	numList        = mclapply(ACGTList, acgtToNum, mc.cores=NCORES)

	M = as.data.frame (numList,col.names=colnames (markers), check.names=F)

	newGeno = data.frame (Markers=geno[,1], M, check.names=F) # Check=FALSE names but not row names
	newName = paste0 (strsplit (kmatrixFile, split="[.]")[[1]][1], "-GENODIVE.csv")
	write.csv (file=newName, newGeno, quote=F, row.names=F)

	return (newName)
}

#----------------------------------------------------------
# Convert gwaspoly genotype from ACGT to numeric format
#----------------------------------------------------------
ACGTToNumericGenotypeFormat <- function (genotypeFile, ploidy, MAP=F) 
{
	NCORES = ifelse (detectCores()>1, detectCores()-1, 1)
	geno = read.csv (file=genotypeFile, header=T, check.names=F)

	markers           = as.matrix(geno[,-(1:3)])
	sampleNames       = colnames (geno[,-(1:3)])
	rownames(markers) = geno[,1]

	tmp     <- apply(markers,1,getReferenceAllele)
	map     <- data.frame(Marker=geno[,1],Chrom=factor(geno[,2],ordered=T),Position=geno[,3], 
						  Ref=tmp[1,], Alt=tmp[2,], stringsAsFactors=F)
	# cache change: Remove fixed out dir
	#write.csv (map, "out/map.csv", quote=F, row.names=F)

	map$Ref <- tmp[1,]
	map$Alt <- tmp[2,]

	#>>>> Convert all genotypes of a marker row from ACGT to Num
	acgtToNum <- function(x){
		x   <- unlist (x)
		y   <- gregexpr(pattern=x[1],text=x[-1],fixed=T)
		ans <- as.integer(lapply(y,function(z){ifelse(z[1]<0,ploidy,ploidy-length(z))}))
		return(ans)
	}
	
	# Convert All ACGT matrix to Num
	matRefMarkers   = cbind(map$Ref,markers)
	matTransposed   = t(matRefMarkers)
	#ACGTList        = mclapply(seq_len(ncol(matTransposed)), function(i) matTransposed[,i],mc.cores=NCORES)
	ACGTList        = lapply(seq_len(ncol(matTransposed)), function(i) matTransposed[,i])
	numList         = mclapply(ACGTList, acgtToNum, mc.cores=NCORES)
	
	numericMatrixM  = as.data.frame (numList,col.names=rownames (matRefMarkers), optional=F, make.names=F)
	colnames (numericMatrixM) = rownames (matRefMarkers)

	tM              =  (t(numericMatrixM))
	colnames (tM)   = sampleNames

	newGeno     = data.frame (map[,1:3], tM, check.names=F) # Check=FALSE names but not row names

	if (MAP == TRUE)
		return (list(geno=geno, genoNum=newGeno, map=map))
	else {
		newGenoFile = addLabel (genotypeFile, "NUM")
		write.csv (newGeno, newGenoFile, quote=F, row.names=F)
		return (newGenoFile)
	}
}

#----------------------------------------------------------
# Convert numeric (gwaspoly) genotype to ACGT using solcap ref/alt alleles
# Basic genotype: [Markers+Alleles]
# Warning!!! It takes too long
#----------------------------------------------------------
numericToACGTFormatGenotype <- function (genotypeFile, mapFile) 
{
	genotype     <- read.csv (genotypeFile, header=T, check.names=F)
	map         <- read.csv (mapFile, header=T, check.names=F)
	rownames (map) <- map [,1]
	alleles      <- genotype [,-c(2,3)]
	allelesACGT  <- numericToACGTFormatAlleles (alleles, map)
	genoACGT     = data.frame (genotype [,1:3], allelesACGT [,-1])

	#outFile = paste0 (strsplit (genotypeFile, split="[.]")[[1]][1],"-ACGT.csv")
	outFile = addLabel (genotypeFile, "ACGT")
	write.csv (genoACGT, outFile, quote=F,row.names=F)
	return (outFile)
}

numericToACGTFormatAlleles <- function (alleles, map) 
{
	setA <- function (allelesVec, refs, alts) {
		id  = allelesVec [1]
		gnt <- as.numeric (allelesVec [-1])
		ref = refs [id,2]
		alt = alts [id,2]
		gnt [gnt==4] = strrep(alt,4)
		gnt [gnt==3] = paste0 (strrep(ref,1),strrep(alt,3))
		gnt [gnt==2] = paste0 (strrep(ref,2),strrep(alt,2))
		gnt [gnt==1] = paste0 (strrep(ref,3),strrep(alt,1))
		gnt [gnt==0] = strrep(ref,4)
		return (gnt)
	}
	refs <- data.frame (map [,c(1,2)])
	rownames (refs) <- map [,1]
	alts <- data.frame (map [,c(1,3)])
	rownames (alts) <- map [,1]
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
	genotype   <- read.csv (genotypeFile, header=T, check.names=F)
	map      = genotype [,1:3]
	alleles  = genotype [,-c(2,3)]

	allelesNum = numericTetraToNumericDiploMatrix (alleles)	

	newGeno = cbind (map, allelesNum[,-1])
	newName = addLabel (genotypeFile, "diploNUM")
	write.csv (file=newName, newGeno, row.names=F, quote=F)
}

# Convert tetra numeric matrix to diplo numeric matrix
numericTetraToNumericDiploMatrix <- function (tetraMatrix) {
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

	diploMatrix <- t (apply (tetraMatrix, 1, toDiplo))
	colnames (diploMatrix) = colnames (tetraMatrix [-1])
	rownames (diploMatrix) = rownames (tetraMatrix)
	return (data.frame (tetraMatrix[,1, drop=F], diploMatrix, check.names=F))
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

##-------------------------------------------------------------
# Convert gwaspoly genotye from numeric tetra to numeric diplo for GS (-1,0,1)
#-------------------------------------------------------------
numericTetraToNumericDiploGSMatrix <- function (genotypeFileMatrix) 
{
	toDiplo <- function (alleles) {
		alleles [alleles==0] = -1
		alleles [alleles==1] = 0
		alleles [alleles==2] = 0
		alleles [alleles==3] = 0
		alleles [alleles==4] = 1
		return (alleles)
	}

	genotype   <- read.csv (genotypeFileMatrix, header=T, check.names=F)
	markers  = as.character (genotype [,1])
	alleles  = genotype [,-1]

	allelesNum <- t (apply (alleles, 1, toDiplo))
	#colnames (allelesNum) = colnames (alleles )
	#rownames (allelesNum) = rownames (alleles)

	newGeno = cbind (markers, allelesNum)
	newName = addLabel (genotypeFileMatrix, "diploNUM-GS")
	write.csv (file=newName, newGeno, row.names=F, quote=F)
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
# Run a command string using system function and writes output to log file
#-------------------------------------------------------------
runCommand <- function (command, logFile="gwas.log") 
{
	if (DEBUG==T) {
		msgmsg (">>>> ", command)
		system (command)
	}else {
		errorsLog = paste0 (strsplit(logFile, split="[.]")[[1]], ".errors")
		if (logFile=="gwas.log")
			system (paste0 (command, " >> ", logFile," 2>> ",errorsLog))
		else
			system (paste0 (command, " > ", logFile," 2> ",errorsLog))
	}
}

#----------------------------------------------------------
# Convert fitpoly file to genotype matrix
#----------------------------------------------------------
commandArgsFitpoly <- function (args) {
	msg ("Converting fitPoly data table to genotype matrix format...")
	#--- Function to convert fitPoly table to GWASpoly matrix
	createMatrix <- function (snp, fitGenos) {
		dataSnp = fitGenos [fitGenos$MarkerName==snp,]
		genos = data.frame (t(dataSnp$geno))
		if (all (is.na (genos))) 
			return (NULL)
		names (genos) = dataSnp$SampleName
		rownames (genos) = snp
		return (genos)
	 }
	#------------------------------------------------------

	fitGenosFile = args[1]
	fitGenos     = read.table (fitGenosFile, sep="\t", header=T)		
	nCores       = ifelse (detectCores()>1, detectCores()-1, 1)

	outs         = mclapply (snpList, createMatrix, fitGenos, mc.cores=detectCores())
	fitGenosDF   = do.call (rbind.data.frame, outs)

	outFilename  = addLabel (fitGenosFile, "MATRIX")

	fitPolyMatrix = data.frame (Makers=rownames (fitGenosDF), fitGenosDF)
	write.csv (fitPolyMatrix, outFilename, quote=F, row.names=F)
}

#----------------------------------------------------------
# Util to print head of data
# Fast head, for debug only
#----------------------------------------------------------
view <- function (data, n=5,m=6, outPrefix="") {
	filename = deparse (substitute (data))
	name = paste (deparse (substitute (data)),":  ")
	if (is.null (dim (data))) {
		dimensions = paste (length (data))
		message (name, class(data), " : (", paste0 (dimensions),")")
		if (length (data) < 6) n = length(data)
		print (data[1:n])
	}else {
		dimensions = paste0 (unlist (dim (data)),sep=c(" x ",""))
		message (name, class(data), " : (", paste0 (dimensions),")")
		if (n==0 || nrow (data) < 5) n = nrow(data)
		if (m==0 || ncol (data) < 6) m = ncol(data)
		print (data[1:n,1:m])
	}
	write.csv (data, paste0("x-", outPrefix,"-", filename, ".csv"), quote=F, row.names=T)
}
#----------------------------------------------------------
#----------------------------------------------------------
#main ()
