#!/usr/bin/Rscript

#-------------------------------------------------------------
#-------------------------------------------------------------
main <- function () {
	HOME="/home/lg/agrosavia/GWAS-TOOL/multiGWAS-dev"

	source (paste0(HOME, "/main/gwas-lib.R"))
	source (paste0(HOME, "/main/gwas-preprocessing.R"))
	msgmsg <<- message


	args = commandArgs(trailingOnly = TRUE)
	args = c("example-genotype-tetra-gwaspoly-ACGT.csv", "example-phenotype-single-trait.csv")

	params = list()
	params$genotypeFile  = args [1]
	params$phenotypeFile = args [2]
	params$gwasModel     = "Naive"
	params$correctionMethod  = "Bonferroni"
	params$significanceLevel = 0.05
	params$geneAction = "additive"

	res = runToolGapit (params)
	write.table (res$scores, "scores-gapit.csv", quote=F, row.names=F)
}

#-------------------------------------------------------------
#-------------------------------------------------------------
runToolGapit <- function (params) {

	if (DEBUG==F) sink ("log-GAPIT-outputs.log")

	scoresFile = paste0 ("out/tool-GAPIT-scores-", params$gwasModel, ".csv")
	if (params$geneAction=="additive") {
		scoresMgwas = runGapit ("additive", params$genotypeFile, params$phenotype, scoresFile, params)
	}else if (params$geneAction=="dominant") {
		scoresMgwas = runGapit ("dominant", params$genotypeFile, params$phenotype, scoresFile, params)
	}else if (params$geneAction %in% c("all", "automatic")){
		res = mclapply (c("additive","dominant"), runGapit, params$genotypeFile, params$phenotype, scoresFile, params, mc.cores=NCORES)
		scoresMgwas  = do.call (rbind.data.frame, res)
	}

	write.table (scoresMgwas, scoresFile, sep="\t", row.names=F, quote=F)
	if (DEBUG==F) sink ()

	return (list (tool="GAPIT", scoresFile=scoresFile, scores=scoresMgwas))
}

#-------------------------------------------------------------
#-------------------------------------------------------------
runGapit <- function (geneAction, genotypeFile, phenotypeFile, scoresFile, params) {
	#source (paste0 (HOME, "/main/gapit_functions.R"), echo=F)      # Module with functions to convert between different genotype formats 
	library (compiler)
	source (paste0 (HOME, "/main/gwas-gapit_functions.R"), echo=F)      # Module with functions to convert between different genotype formats 

	#gapit = gwaspolyToGapitFormat (genotypeFile, phenotypeFile, geneAction, FILES=F)
	genotype  = read.csv (params$gapitGenotypeFile)
	phenotype = read.csv (params$gapitPhenotypeFile)
	map       = read.csv (params$gapitMapFile)

	if (tolower (params$gwasModel)=="full") 
		gapitModel = "MLM"
	else 
		gapitModel = "GLM"

	#out <- GAPIT(Y=phenotype, GM=map, GD=genotype, model="MLM", file.output=F)#, kinship.algorithm="None")
	out <- GAPIT(Y=phenotype, GM=map, GD=genotype, model=gapitModel, file.output=F, QC=T)#, kinship.algorithm="None")
	scoresGapit = out$GWAS

	# Write source scores
	scoresGapitFile = addLabel (scoresFile, paste0("SOURCE-",geneAction))
	write.csv (scoresGapit, scoresGapitFile, quote=F, row.names=F)

	# Create multiGWAs scores
	scoresMgwas  = createMultigwasTableFromGapit (scoresGapit, geneAction, params)
	scoresMgwas  = scoresMgwas [order (scoresMgwas$DIFF, decreasing=T),]

	return (scoresMgwas)
}

#-------------------------------------------------------------
#-------------------------------------------------------------
createMultigwasTableFromGapit <- function (scoresGapit, geneAction, params) {
	scoresColumns = c("MODEL", "GC", "Marker", "CHR", "POS", "P", "SCORE", "THRESHOLD", "DIFF")

	MODEL     = geneAction

	Marker    = as.character (scoresGapit$SNP)
	pValues   = scoresGapit$P.value
	adj       = adjustPValues (params$significanceLevel, pValues, params$correctionMethod)
	P         = adj$pValues

	THRESHOLD = round (adj$threshold, 6)
	SCORE     = round (-log10 (P), 6)
	GC        = calculateInflationFactor (SCORE)$delta
	CHR       = scoresGapit$Chromosome
	POS       = scoresGapit$Position
	DIFF      = SCORE - THRESHOLD

	multiGwasTable = data.frame (MODEL=MODEL, GC=GC, Marker=Marker, CHR=CHR, POS=POS, P=P, SCORE=SCORE, THRESHOLD=THRESHOLD, DIFF=DIFF)

	return (multiGwasTable)
}

#-------------------------------------------------------------
#-------------------------------------------------------------
#main ()

