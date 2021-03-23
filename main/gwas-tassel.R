
#-------------------------------------------------------------
# Run Tassel pipeline (GDM and MLM)
#-------------------------------------------------------------
runToolTassel <- function (params) {
	model = params$gwasModel

	# Parameters for the scripts
	#inGenoVCF     = "out/filtered-tassel-genotype.vcf"
	#inPhenoTBL    = "out/filtered-tassel-phenotype.tbl"
	inGenoVCF     = params$tasselGenotypeFile
	inPhenoTBL    = params$tasselPhenotypeFile
	outFile       = paste0 ("out/tool-TASSEL-scores-", model)
	scoresFile    = paste0 (outFile, ".csv")

	if (model=="naive") {
		cmm=sprintf ("%s/main/scripts/script-tassel-NaiveModel.sh %s %s %s", HOME,inGenoVCF, inPhenoTBL, outFile)
		runCommand (cmm, "log-tassel.log")
		tasselFile   = list.files("out/", pattern=sprintf("^(.*(%s).*(1).*(txt)[^$]*)$",model), full.names=T)
	}else if (model=="full") {
		cmm=sprintf ("%s/main/scripts/script-tassel-FullModel.sh %s %s %s", HOME, inGenoVCF, inPhenoTBL, outFile)
		runCommand (cmm, "log-tassel.log")
		tasselFile   = list.files("out/", pattern=sprintf("^(.*(%s).*(stats).*(txt)[^$]*)$",model), full.names=T)
	}
	
	# Rename output file
	msgmsg ("Tassel output file: ", tasselFile)

	# Create table by group of p-values: p, pAdd, and pDom
	results      <- read.table (file=tasselFile, header=T, sep="\t", check.names=F)
	# Returns a new table according to "varP" values
	createTableTassel  <- function (results, model, varP) {
		allP    = results [,varP]
		results = results [complete.cases (allP),]
		pValues = results [,varP]

		#scores = -log10(results [,varP])
		#m      = length (scores)


		# Adjust scores, threshold using Bonferroni or FDR
		adj       = adjustPValues (params$significanceLevel, pValues, params$correctionMethod)
		pValues   = adj$pValues
		threshold = round (adj$threshold, 6)
		scores    = round (-log10 (pValues), 6)

		# Obsolete
		#if (params$correctionMethod=="FDR") 
		#	threshold    <- calculateThreshold (level=SIGNIFICANCE_LEVEL, scores=scores, method="FDR")
		#else if  (params$correctionMethod=="Bonferroni") 
		#	threshold    <- -log10 (SIGNIFICANCE_LEVEL / m)
			#threshold    <- calculateThreshold (level=SIGNIFICANCE_LEVEL, scores=scores, method="Bonferroni")

		# Compose the results table
		GC = calculateInflationFactor (scores)
		pTable  <-  cbind (MODEL=model, GC=GC$delta, results [, c("Marker", "Chr", "Pos")], P=pValues, 
						   SCORE=scores, THRESHOLD=threshold, DIFF=(scores-threshold))
		pTable  <- pTable [order (scores, decreasing=T),]    
		return (pTable)
	}
	# Create tables for each gene action

	# Join the tables, order by DIFF, and write it
	if (params$geneAction == "additive") 
		scoresTableAll <- createTableTassel (results, "additive", "add_p")
	else if (params$geneAction == "general") 
		scoresTableAll <- createTableTassel (results, "general", "p")
	else if (params$geneAction == "dominant") 
		scoresTableAll <- createTableTassel (results, "dominant", "p")
	else if (params$geneAction == "all") { 
		addTable <- createTableTassel (results, "additive", "add_p")
		domTable <- createTableTassel (results, "general", "p")
		gnrTable <- createTableTassel (results, "dominant", "p")
		scoresTableAll <- rbind (addTable, domTable, gnrTable)
	}
	colnames (scoresTableAll)[colnames(scoresTableAll) %in% c("Chr","Pos")] = c ("CHR","POS")
	write.table (scoresTableAll, scoresFile, quote=F, sep="\t", row.names=F)
	msg ("... Ending TASSEL")

	return (list (tool="TASSEL", scoresFile=scoresFile, scores=scoresTableAll))
}

