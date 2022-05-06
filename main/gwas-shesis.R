
#-------------------------------------------------------------
# SHEsis tool
# Shesis uses a low significance level of 0.01 as their pvalues are inflates
# LOG: 
#      v1.01: Relatedness using kinship matrix from GWASpoly (>0.9) 
#-------------------------------------------------------------
runToolShesis <- function (params) 
{
	geneAction = params$geneAction
	#if (geneAction != "automatic" & geneAction != "additive")
	#	return ()

	geneAction = "additive"

	if (params$gwasModel == "naive") {
		dataShesis = gwaspolyToShesisGenoPheno (params$genotypeFile, params$phenotypeFile, params$ploidy, params$traitType)
	} else if (params$gwasModel == "full") {
		closestIndividuals   = getClosestIndividuals ("out/GWASpoly-kinship-matrix.csv", params)

		genotype             = read.csv (file=params$genotypeFile, header=T, check.names=F)
		phenotype            = read.csv (file=params$phenotypeFile, header=T, check.names=F)
		genotypeKinship      = genotype [, !(colnames (genotype) %in% closestIndividuals)]
		phenotypeKinship     = phenotype [!(phenotype [,1] %in% closestIndividuals),]

		genotypeFileKinship  = addLabel (params$genotypeFile,  "kinship-shesis")
		phenotypeFileKinship = addLabel (params$phenotypeFile, "kinship-shesis")
		write.table (file=genotypeFileKinship,  genotypeKinship,  row.names=F, quote=F, sep=",")
		write.table (file=phenotypeFileKinship, phenotypeKinship, row.names=F, quote=F, sep=",")

		dataShesis = gwaspolyToShesisGenoPheno (genotypeFileKinship, phenotypeFileKinship, params$ploidy, params$traitType)
	}else {
		message ("ERROR: GWAS model: ", params$gwasModel, " not supported")
		quit ()
	}

	shesis  = runShesisCommand (params$traitType, dataShesis$genoPhenoFile, dataShesis$markersFile, geneAction, params)

	return (list (tool="SHEsis", scoresFile=shesis$scoresFile, scores=shesis$scoresMgwas))
}

#-------------------------------------------------------------
# Get closest relatives using kinship matrix from GWASpoly
#-------------------------------------------------------------
getClosestIndividuals <- function (kinshipMatrixFile, params) {
	#----- Flatten kinship matrix to pairs of individuals ----
	flattenKMat <- function(kmat) {
	  ut <- upper.tri(kmat)
	  data.frame( IND1 = rownames(kmat)[row(kmat)[ut]],
				  IND2 = rownames(kmat)[col(kmat)[ut]],
				  K  =(kmat)[ut], stringsAsFactors=F)
	}
	#------ Create kinship matrix using GWASpoly -------------
	createKinshipMatrix <- function (kinshipMatrixFile, params) {
		data1 = read.GWASpoly (ploidy = params$ploidy, pheno.file = params$phenotypeFile, 
								geno.file = params$genotypeFile, format = "ACGT", n.traits = 1, delim=",")
		data2 = set.K (data1)
		kmat  = data2@K  
		write.csv (kmat, kinshipMatrixFile, quote=F)
		return (kmat)
	}
	#---------------------------------------------------------
	if (!file.exists (kinshipMatrixFile)) 
		kmat = createKinshipMatrix (kinshipMatrixFile, params)
	else 
		kmat = read.csv (kinshipMatrixFile, row.names=1)

	kTable = flattenKMat (kmat)
	closest = kTable$IND2 [abs (kTable$K) > 0.9 ]
	return (closest)
}

#-------------------------------------------------------------
#-------------------------------------------------------------
runShesisCommand <- function (traitType, genoPhenoFile, markersFile, geneAction, params) 
{
	#genoPheno  = "out/filtered-shesis-genopheno.tbl"
	#inMarkers    = "out/filtered-shesis-markernames.tbl"
	scoresFile   = paste0 ("out/tool-SHEsis-scores-", params$gwasModel)
	outFile      = paste0 (scoresFile, "-SOURCES")
	flagQTL      = ifelse (traitType=="quantitative", "--qtl", "")

	cmm=sprintf ("%s/main/scripts/script-shesis-associations-qtl.sh %s %s %s %s %s", HOME, genoPhenoFile, params$ploidy, markersFile, outFile, flagQTL)
	runCommand (cmm, "log-SHEsis.log")	

	if (traitType=="quantitative")
		scoresMgwas = createTableFromQuantitativeResults (outFile, params, geneAction)
	else # case-control
		scoresMgwas = createTableFromBinaryResults (outFile, params, geneAction)

	scoresMgwasFile = paste0 (scoresFile, ".csv")
	write.table (scoresMgwas, scoresMgwasFile, row.names=F, quote=F, sep="\t")

	return (list (tool="SHEsis", scoresFile=scoresMgwasFile, scoresMgwas=scoresMgwas))
}
#-------------------------------------------------------------
# Create table from quantitative .txt file 
#-------------------------------------------------------------
createTableFromQuantitativeResults <- function (outFile, params, geneAction) {
	resultsFile = paste0 (outFile, ".txt")
	results  = read.table (file=resultsFile, header=T, sep="\t", check.names=T) # TRUE as SHEsis colnames have spaces

	# LG: Added 1e-10 to avoid "inf" values in scores
	#pValues  = results[,"P.value"] + 1e-10
	pValues  = results[,"P.value"] 
	m        = length (pValues)
	adj       = adjustPValues (params$significanceLevel, pValues, params$correctionMethod)
	pValues   = adj$pValues
	threshold = adj$threshold
	scores    = -log10 (pValues)

	# Set Columns
	#map <- read.csv ("out/map.csv", check.names=F)
	map <- read.csv (params$mapFile, check.names=F)
	rownames (map) = map [,1]
	Marker <- as.character (results$SNP)
	CHR <- map [Marker, 2]
	POS <- map [Marker, 3]
	GC  = calculateInflationFactor (scores)

	model = geneAction
	resultsAll <- data.frame (MODEL=model, GC=GC$delta, Marker, CHR, POS, P=pValues, SCORE=round (scores,6), THRESHOLD=round (threshold,6), DIFF=round (scores-threshold, 6))
	resultsAll <- resultsAll [order (resultsAll$DIFF, decreasing=T),]

	return (resultsAll)
}

#-------------------------------------------------------------
# Create table from binart .txt results file
#-------------------------------------------------------------
createTableFromBinaryResults <- function (outFile, params, geneAction) {
	resultsFile = paste0 (outFile, ".txt")

	conn = file (resultsFile, open="r")
	lines = readLines (conn)
	close (conn)
	n  = length (lines) 

	inAllele=F; genotypeList = list (); allelesList  = list (); currentMarker = NULL
	for (i in 1:n)  {
		li = lines [i]
		if (grepl ("Allele", li)) {
			inAllele=T
			marker        = strsplit (li, "(", fixed=T)[[1]][1]
			currentMarker = c(SNP=marker)
		}else if (grepl ("Genotype", li)) {
			inAllele=F
			marker        = strsplit (li, "(", fixed=T)[[1]][1]
			currentMarker = c(SNP=marker)
		}else if (grepl ( "Pearson", li)) {
			p = strsplit (li, " ")
			currentMarker = c (currentMarker, pPearson=p[[1]][length(p[[1]])])
			if (inAllele)
				allelesList = append (allelesList, currentMarker)
			else
				genotypeList = append (genotypeList, currentMarker)
		}
	}

	columns = c("SNP", "pPearson")
	allelesMat   = matrix (allelesList, ncol=2, byrow=T)
	colnames (allelesMat) = columns
	write.csv (allelesMat, "shesis-alleles-pvalues.csv", quote=F, row.names=F)

	results   = matrix (genotypeList, ncol=2, byrow=T)
	colnames (results) = columns
	write.csv (results, "shesis-genotype-pvalues.csv", quote=F, row.names=F)
	results = read.csv ("shesis-genotype-pvalues.csv")

# LG: Added 1e-10 to avoid "inf" values in scores
	#pValues  = results[,"P.value"] + 1e-10

	pValues  = results[,"pPearson"] 
	m        = length (pValues)
	adj       = adjustPValues (0.01, pValues, params$correctionMethod)
	message ("Formating 0 ...")
	pValues   = adj$pValues
	message ("Formating 1 ...")
	threshold = adj$threshold
	message ("Formating 2 ...")
	scores    = -log10 (pValues)

	# Set Columns
	#map <- read.table (file="out/map.tbl", sep="\t", check.names=F)
	map <- read.csv (params$mapFile, check.names=F)
	rownames (map) = map [,1]
	SNP <- as.character (results$SNP)
	CHR <- map [SNP, 2]
	POS <- map [SNP, 3]
	GC  = calculateInflationFactor (scores)

	model = geneAction
	resultsAll <- data.frame (MODEL=model, GC=GC$delta, SNP, CHR, POS, P=pValues, SCORE=round (scores,6), THRESHOLD=round (threshold,6), DIFF=round (scores-threshold, 6), results)
	resultsAll <- resultsAll [order (resultsAll$DIFF, decreasing=T),]

	return (resultsAll)
}



#----------------------------------------------------------
# Transform table genotype to SHEsis genotype format
#----------------------------------------------------------
gwaspolyToShesisGenoPheno <- function (genotypeFile, phenotypeFile, ploidy, traitType) 
{
	msgmsg ("    >>>> Writting gwasp to shesis genopheno...")
	sep <- function (allele) {
		s="";
		for (i in 1:ploidy) s=paste0(s, substr (allele,start=i,stop=i)," ");
		#s = paste (strsplit (allele, "")[[1]], collapse=" ")
		return (s)
	}
	geno    = read.csv (file=genotypeFile, stringsAsFactors=F, check.names=F)
	pheno   = read.csv (file=phenotypeFile, stringsAsFactors=F, check.names=F)
	if (traitType=="case-control") {
		pheno [pheno[,2]==1,2]=2
		pheno [pheno[,2]==0,2]=1
	}

	rownames (pheno) <- pheno [,1]
	map        <- geno  [,c(1,2,3)]    # Get SNP, Cromosome, Position
	rownames (geno)  <- map [,1] 

	alleles    <- geno[,-c(1,2,3)]
	alleles [is.na(alleles)] = paste (rep ("0", ploidy), collapse="") # NAs as "00" or "0000"

	allelesMat <- t(sapply (alleles, sep))

	samples         = rownames (allelesMat)
	pheno           = pheno [samples,]
	pheno [,2]      = impute.mode (pheno [,2])
	genoPhenoShesis = data.frame (Sample=pheno[,1], Trait=pheno[,2],  allelesMat)

	msgmsg ("    >>>> Writing SHEsis genopheno...")
	#outFile = "out/filtered-shesis-genopheno.tbl"
	outFileGenoPheno = addLabel (genotypeFile, "SHESIS-GENOPHENO")
	write.table (genoPhenoShesis, outFileGenoPheno, quote=F,row.names=F,col.names=F, sep="\t")

	msgmsg ("    >>>> Writing SHEsis marker names...")
	#outFile = "out/filtered-shesis-markernames.tbl"
	outFileMarkerNames = addLabel (genotypeFile, "SHESIS-MARKERNAMES")
	write.table (map[,1], outFileMarkerNames, quote=F,row.names=F,col.names=F, sep="\t")

	return (list(genoPhenoFile=outFileGenoPheno, markersFile=outFileMarkerNames))
}

