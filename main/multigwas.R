#!/usr/bin/Rscript
DEBUG = F; SILENT=T

options (width=300, error=traceback)
if (DEBUG) {
	SILENT <- FALSE;options (warn=2)
	message ("Running in debug mode...")
}

# Get environment GWAS variable 
HOME   <<- Sys.getenv ("MULTIGWAS_HOME")
OUTDIR <- getwd()
params <- list ()
LOGS   <- list ()

# INFO  : Tool for running GWAS integratind four GWAS tools: 
#         GWASpoly and SHEsis for polyploids species, and Plink and Tassel for diploids.
# AUTHOR: Luis Garreta (lgarreta@agrosavia.co) 
# DATE  : 12/feb/2020
# LOGS  :   
	# r1.5: Added cacheDir. Modified ACGTToNumber. Modified report LD. 
	# r1.3: Improved check input files. Added try catch for running gwas tools.
	# r1.1: Added create summary table with scored SNPs. Deleted config and pheno files
	# r1.0: Added genotype formats, LD, and HWE
	# r0.8: Add linkage disequilibrium analysis
	# r0.5: Using VCF files
	# r0.3: Added column gene action model to tables of results by tool

#-------------------------------------------------------------
# Return string with usage instructions
#-------------------------------------------------------------
usageInstructions <- function () {
	USAGE="USAGE: multigwas <config file>\n"
	return (USAGE)
}

#-------------------------------------------------------------
# Main for multi traits
#-------------------------------------------------------------
main <- function () {
	message ("MultiGWAS 1.5 is running!\n")
	args = commandArgs(trailingOnly = TRUE)
	if (length (args) < 1) {
		message (usageInstructions())
		quit ()
	}

	# Check dependencies and init environments
	checkDependencies (params)
	initGlobalsLoadSources ()

	# Read and check config file arguments
	configFile = args [1]
	params <<- readCheckConfigParameters (configFile)
	params <<- initWorkingCacheDir (configFile)

	printParameters (params)

	# Preprocess input files, convert to tools format, and create config files for traits
	params <<- convertGenoPhenoToMultigwas (params$cacheDir) 
  	data = filterGenoPhenoMapFiles (params$genotypeFile, params$phenotypeFile, params$cacheDir)
	params <<- convertGenotypeToToolFormats (data$genotypeFile, params$cacheDir)

	# Analyze each trait independiently
	mainDir        = getwd ()
	bestSNPsList   = list()
	allSNPsList    = list()
	allGScoreList  = list()
	allPhenotypes  = read.csv (data$phenotypeFile, row.names=1)
	traitNames     = colnames (allPhenotypes)
	sampleNames    = row.names (allPhenotypes)
	for (traitName in traitNames) {
		createDir (traitName)
		setwd (traitName)

		phenotype            = cbind (sampleNames, allPhenotypes [,traitName])
		colnames (phenotype) = c ("NAMES", traitName)
		traitConfigFile      = createConfigFileForTrait (phenotype, traitName)

		resultsTables       = mainSingleTrait (traitConfigFile)

		bestSNPsList  = append (bestSNPsList, list (resultsTables$best))
		allSNPsList   = append (allSNPsList, list (resultsTables$all))
		allGScoreList = append (allGScoreList, list (resultsTables$gscore))
		write.csv (resultsTables$gscore, sprintf ("%s/%s-Scores.csv", mainDir, traitName))
		setwd (mainDir)
	}
	bestSNPsTable   = do.call ("rbind", bestSNPsList)
	allSNPsTable    = do.call ("rbind", allSNPsList)
	allGScoreTable  = do.call ("rbind", allGScoreList)

	createGlobalReport (configFile, bestSNPsTable, allSNPsTable, allGScoreTable) 
	printLOGS ("")
}

#-------------------------------------------------------------
# Create output working dir and copy input files into it
#-------------------------------------------------------------
initWorkingCacheDir <- function (configFile) {
	msg ("Initializing working dir...")
	runningDir   = paste0 ("out-", strsplit (configFile, split="[.]") [[1]][1])
	createDir (runningDir)

	copyFilesToRunningDir (params, runningDir)

	# Change to the output dir and set global runningDir 
	baseDir = getwd()
	setwd (runningDir)

	msgmsg ("Initializing cache dir...")
	if (DEBUG) {
		cacheDir = sprintf ("%s/%s", baseDir, "cacheDir")
		if (dir.exists (cacheDir)==FALSE)
			createDir (cacheDir)
	}else {
		cacheDir = sprintf ("%s/%s", getwd(), "cacheDir")
		createDir (cacheDir)
	}

	params$cacheDir    = cacheDir
	params$runningPath = getwd()
	return (params)
}

#----------------------------------------------------------
# Copy input files to running dir
# Remove path and copy files to trait directory (runningDir)
#----------------------------------------------------------
copyFilesToRunningDir <- function (params, runningDir) {
	runCommand(sprintf ("cp %s %s", params$genotypeFile, runningDir))
	params$genotypeFile  = basename (params$genotypeFile)

	runCommand(sprintf ("cp %s %s", params$phenotypeFile, runningDir))
	params$phenotypeFile = basename (params$phenotypeFile)

	if (tolower (params$genotypeFormat) %in% c("matrix", "fitpoly", "updog")) {
		file.copy (params$mapFile, runningDir)
		params$mapFile = basename (params$mapFile)
	}
}

#-------------------------------------------------------------
# Write global tables for all and best markers for all traits
#-------------------------------------------------------------
createGlobalReport <- function (configFile, bestSNPsTable, allSNPsTable, allGScoreTable){
	msg ("Writing global report: table with all and best markers for all traits...")
	outFile = "SCORES-BEST.csv"
	write.csv (bestSNPsTable, outFile, row.names=F)
	
	#outFile = "ALL-SCORES-ALL.csv"
	#write.csv (allSNPsTable, outFile, row.names=F)

	outFile = "SCORES-ALL.csv"
	write.csv (allGScoreTable, outFile, row.names=F)
}

#-------------------------------------------------------------
# Main for a single trait
#-------------------------------------------------------------
mainSingleTrait <- function (traitConfigFile) {
	msg ("Processing parameters file: ", traitConfigFile)

	# Copy files to working dirs and create output dirs
	paramsTrait   		 = getTraitConfig (traitConfigFile)
	params.phenotypeFile = paramsTrait$phenotypeFile

	# Run the four tools in parallel
	listOfResults = runGWASTools (paramsTrait)

	# Create reports
	msg ("Creating reports (Table, Venn diagrams, Manhattan&QQ plots, SNP profiles)...")
	resultsTables  = createReports (listOfResults, paramsTrait)
	snpTables.best = data.frame (TRAIT=paramsTrait$trait, resultsTables$best)
	allScores      = data.frame (TRAIT=paramsTrait$trait, resultsTables$all)
	allGScores     = data.frame (TRAIT=paramsTrait$trait, resultsTables$gscore)

	# Move out files to output dir
	msg ("Sending trait files to output folders...")
	moveOutFiles (paramsTrait$outputDir, paramsTrait$reportDir, traitConfigFile, params.phenotypeFile, paramsTrait)

	return (list (best=snpTables.best, all=allScores, gscore=allGScores))
}

#-------------------------------------------------------------
# Print logs manually writen during the process
#-------------------------------------------------------------
printLOGS <- function (e) {
	message ("#------------------------------------------------------------")
	message ("#------------------- WARNINGS SUMARY ------------------------")
	message ("#------------------------------------------------------------")
	for (l in LOGS) 
		message (l)
	message ("#----------------- END OF EXECUTION -------------------------")
	message ("#--------------- Moving files to outpud dirs -------------------------")
	sink (sprintf("%s/%s", OUTDIR, "gwas.errors"))
	print (e)
	sink()
}
#-------------------------------------------------------------
#-------------------------------------------------------------
# Used to run in parallel the other functions
#-------------------------------------------------------------
runGWASTools <- function (params) {
	results = NULL
	runOneTool <- function (tool, params) {
		tryCatch ({
			if (tool=="gwaspoly") 
				results = runToolGwaspoly (params)
			else if (tool=="shesis")   
				results = runToolShesis (params)
			else if (tool=="gapit") 
				results = runToolGapit (params) 
			else if (tool=="tassel")   
				results = runToolTassel (params)
			else if (tool=="plink")    
				results = runToolPlink (params)
			else                       
				stop ("Tool not supported")
			return (results)
		},error = function (e){
			text =  "WARNING: One of the four tools without results."
			message (text)
			LOGS <<- append (LOGS, text)
			return (NULL)
		})
	}

	# A string containing the names of the tools to run (e.g. "GWASpoly SHEsis PLINK TASSEL")
	msg ("Preparing to execute in parallel the GWAS tools:") 
	tools = strsplit(params$tools ,split=" ")[[1]]
	for (tool in tools) 
		msgmsg ("Running ", tool)

	listOfResults     = mclapply (tools, runOneTool, params, mc.cores=NCORES, mc.silent=F)
	# Remove NULLs if they were running tool
	listOfResults     = listOfResults [lengths (listOfResults) != 0]
	listOfResultsBest = selectBestGeneActionModelAllTools (listOfResults, params$geneAction, params$nBest)
	listOfResultsLD   = runLinkageDisequilibriumAnalysis (listOfResultsBest, params$nBest, params$genotypeNumFile, params$R2)

	#return (listOfResultsBest)
	return (listOfResultsLD)
}

#-------------------------------------------------------------
# Create configuration files for individual traits in multitrait file
#-------------------------------------------------------------
createConfigFileForTrait <- function (phenotype, traitName) {
	msg ("Creating configuration files for trait: ", traitName)

	phenoTraitFile = paste0 (traitName,".csv")
	write.csv (phenotype, phenoTraitFile, row.names=F)

	params = convertPhenotypeToToolFormats (phenoTraitFile)

	# Create config file for trait
	params$phenotypeFile = phenoTraitFile
	params$trait         = traitName
	params$outputDir     = traitName
	configLines = c()
	for (n in names (params))
		configLines = c (configLines, paste0 (n, " : ", params[n]))

	traitConfigFile  = paste0 (traitName, ".config")
	writeLines (configLines, traitConfigFile)
	return (traitConfigFile)
}

#-------------------------------------------------------------
# Check dependencies before running GWAS tools
#-------------------------------------------------------------
checkDependencies <- function (params) {
		# Check if JAVA is installed
		if (Sys.which ("java")=="") {
				message ("ERROR: Java not found. Java is needed to run TASSEL and GUI interface.")
				quit ()
		}
}
#-------------------------------------------------------------
# Define global variables, load packages, and load main
#-------------------------------------------------------------
initGlobalsLoadSources <- function () {
	source (paste0 (HOME, "/main/gwas-lib.R")) 
	msg ("Loading libraries, setting globals, sourcing files...")

	.libPaths (paste0(HOME, "/opt/Rlibs"))

	# Load packages
	suppressMessages (library (GWASpoly)) #
	suppressMessages (library (parallel)) #
	suppressMessages (library (yaml))  # For read config file

	NCORES <<- ifelse (DEBUG==T, 1, detectCores ())

	# New class for gwaspoly
	setClass ("GWASpolyStruct", slots=c(params="list"), contains="GWASpoly.K")

	# Load main
	source (paste0 (HOME, "/main/gwas-parameters.R"))      # Module for checking configuration parameters
	source (paste0 (HOME, "/main/gwas-preprocessing.R"))      # Module with functions to convert between different genotype formats 
	source (paste0 (HOME, "/main/gwas-reports.R"))            # Module with functions to create summaries: tables and venn diagrams
	source (paste0 (HOME, "/main/gwas-heatmap.R"))            # Module with functions to create heatmaps for shared SNPs
	source (paste0 (HOME, "/main/gwas-gwaspoly.R"))           # Module with gwaspoly functions
	#source (paste0 (HOME, "/main/gwas-plink.R"))              # Module with plink functions
	source (paste0 (HOME, "/main/gwas-tassel.R"))             # Module with tassel functions
	source (paste0 (HOME, "/main/gwas-shesis.R"))             # Module with shesis functions
	source (paste0 (HOME, "/main/gwas-gapit.R"))             # Module with shesis functions
	#source (paste0 (HOME, "/main/gwas-parameters.R"))             # Module with shesis functions
}

#-------------------------------------------------------------
# Read configuration parameters and set working directories
#   - Create dir, if it exists, it is renamed as old-XXXX
#   - Copy and make links of geno/pheno to output dirs
#-------------------------------------------------------------
getTraitConfig <- function (traitParamsFile) {
	paramsTrait = yaml.load_file (traitParamsFile, merge.precedence="order") 
	traitDir    = paramsTrait$outputDir
	#outDir      = paste0 (traitDir, "/out")
	outDir      = "out"

	# Copy files two times to both trait and out dir
	outDirs    = c(outDir)
	for (dir in outDirs) {
		createDir (dir)
		file.copy (traitParamsFile, dir)
		#file.copy (paramsTrait$genotypeFile, dir )
		file.symlink (sprintf ("../%s", paramsTrait$genotypeFile, dir), dir)
		file.copy (paramsTrait$phenotypeFile, dir)
		if (paramsTrait$genotypeFormat %in% c("matrix", "fitpoly", "updog"))
			runCommand (sprintf ("cp %s %s", paramsTrait$mapFile, dir))
	}
	# Change to the working dir and set dirs in paramsTrait
	paramsTrait$outputDir = "out/"
	paramsTrait$reportDir = "report/"
	paramsTrait$traitDir  = traitDir

	return (paramsTrait)
}

#-----------------------------------------------------------
# Select best N SNPs from multiple action models (for GWASpoly and TASSEL)
# Uses three criteria: best GC, best replicability, and best significants
# PLINK also can produce info of more action models using options
#-----------------------------------------------------------
selectBestGeneActionModelAllTools <- function (listOfResults, geneAction, nBest) {
	msg ("Selecting best gene action model for all tools...")
	i = 1
	for (res in listOfResults) {
		msgmsg ("Best gene action for: ", res$tool)
		bestScoresTool   = selectBestGeneActionModelTool (res$scores, nBest, res$tool, geneAction)
		scoresFileBest   = addLabel (res$scoresFile, "BEST")
		write.table (bestScoresTool, scoresFileBest, sep="\t", quote=F, row.names=F)
		listOfResults [[i]]$scoresFile = scoresFileBest
		listOfResults [[i]]$scores     = bestScoresTool
		i = i + 1
	}
	return (listOfResults)
}


selectBestGeneActionModelTool <- function (scoresTool, nBest, tool, geneAction) {
	# Select main columns

	dataSNPs = scoresTool [,c("MODEL", "GC", "Marker", "CHR", "POS", "P", "SCORE", "THRESHOLD", "DIFF")]; 

	if (geneAction=="all") {
		bestScoresTool = distinct (arrange (dataSNPs, -DIFF), Marker, .keep_all=T)
		return (bestScoresTool)
	}

	if (geneAction %in% c("all", "additive", "dominant", "general") || tool %in% c("SHEsis")){
		bestScoresTool = distinct (arrange (dataSNPs, -DIFF), Marker, .keep_all=T)
		return (bestScoresTool)
	}
		#return (dataSNPs)

	# Order by nBest, DIFF, GC
	orderedSNPs = dataSNPs [order (dataSNPs$MODEL,-dataSNPs$DIFF),]; 

	for (N in c(200, 100, 50, nBest)) {
		# Reduce to groups of nBest
		groupedSNPs = Reduce (rbind, by(orderedSNPs, orderedSNPs["MODEL"], head, n=N)); 

		# Add Count of SNPs between groups
		countedSNPs   = data.frame (add_count (groupedSNPs, Marker, sort=T, name="nSharedSNPs")); 
		# Add count of significatives
		countedSignificantSNPs   = data.frame (add_count (countedSNPs [countedSNPs$DIFF >0, ], MODEL, name="nSign", .drop=F))
		# Add count of shared SNPs
		countedSharedSNPs = aggregate (x=countedSNPs$nSharedSNPs, by=list(MODEL=countedSNPs$MODEL, GC=countedSNPs$GC), 	FUN=sum)
		colnames (countedSharedSNPs) = c("MODEL", "GC", "SHAREDNSPS")

		# Add fraction of shared SNPs between all models
		summMdlSign = cbind (countedSharedSNPs, nSIGN=0)
		rownames (summMdlSign) = summMdlSign [,1]

		if (length (countedSignificantSNPs$nSharedSNPs)==0)
			summMdlSign [as.character (countedSignificantSNPs$MODEL),"nSIGN"] = 0
		else
			summMdlSign [as.character (countedSignificantSNPs$MODEL),"nSIGN"] = countedSignificantSNPs$nSharedSNPs / sum (countedSignificantSNPs$nSharedSNPs)

			
		# Calculate best model score
		totalNs     = length (summMdlSign$MODEL) * N
		scoreGC     = 1 - abs (1-summMdlSign$GC)
		scoreShared = summMdlSign$SHAREDNSPS/totalNs  
		scoreSign   = summMdlSign$nSIGN 

		modelScore  = scoreGC + scoreShared + scoreSign
		summScores  = cbind (countedSharedSNPs, scoreGC, scoreShared, scoreSign, score=modelScore)
		summScores  = summScores [order (summScores$score, summScores$MODEL, decreasing=T),]

		outFilename = addLabel ("out/tmp-bestModel.csv", sprintf ("%s-%0.3d", tool, N))
		write.csv (summScores, file=outFilename, quote=F, row.names=F)
	}

	bestModel = summScores [1, "MODEL"]

	# Select SNPs for model and sort by DIFF
	bestScoresTool = dataSNPs [dataSNPs[,"MODEL"] %in% bestModel,]
	bestScoresTool = bestScoresTool [order (-bestScoresTool$DIFF),]

	return (bestScoresTool)
}
#-----------------------------------------------------------------------
# Calculates LD for each Tool's SNP and rename pairs of SNPs
# added the new file to the listOfResults
#-----------------------------------------------------------------------
runLinkageDisequilibriumAnalysis <- function (listOfResults, nBest, genotypeNumFile, R2) {
	msg ("Linkage disequilibrium analysis...")
	ldTable = createLinkageDisequilibriumTable (listOfResults, nBest, genotypeNumFile, R2) 
	if (is.null (ldTable))
		return (listOfResults)

	i = 1
	for (res in listOfResults) {
		msgmsg ("LD in ", res$tool)
		scoresLD = res$scores
		rownames (scoresLD) = scoresLD$Marker
		scoresLD$Marker     = as.character (scoresLD$Marker)
		for (k in 1:nrow (ldTable)) {
			snp1    = ldTable [k, "SNP1"]
			snp2    = ldTable [k, "SNP2"]
			newName = as.character (ldTable [k, "LD_SNP"])
			scoresLD [scoresLD$Marker %in% snp1 ,"Marker"] = newName
			scoresLD [scoresLD$Marker %in% snp2 ,"Marker"] = newName
		}
		scoresLD     = arrange (scoresLD, -DIFF)
		scoresLD     = scoresLD [!duplicated (scoresLD$Marker),]
		scoresLD     = scoresLD [1:nBest,]
		scoresLD     = scoresLD [grepl("LD", scoresLD$Marker),]
		if (nrow (scoresLD) > 0) {
				scoresLDFile = addLabel (res$scoresFile, "LD")
				write.table (scoresLD, scoresLDFile, sep="\t", quote=F, row.names=F)
				listOfResults [[i]]["scoresLDFile"] = scoresLDFile
		}
		i = i+1
	}
	return (listOfResults)
}

#-------------------------------------------------------
# Create a table for SNPs in all tools in high LD
#-------------------------------------------------------
createLinkageDisequilibriumTable <- function (listOfResults, nBest, genotypeNumFile, R2) {
	# Join all SNPs in one vector
	allSNPs = NULL
	for (res in listOfResults) {
		msgmsg ("Joining results from tool: ", res$tool)
		markers = res$scores$Marker
		N = length (markers)
		if (N >= nBest) N=nBest
		toolSNPs = markers [1:N] 
		allSNPs  = union (allSNPs, toolSNPs)
	}

	# Run linkage disequilibrium
	msgmsg ("Running linkage disequilibrium...")
	genotypeNum = read.csv (genotypeNumFile, row.names=1);
	genomat     = as.matrix (genotypeNum [,-1:-2]); 
	genomatSNPs = genomat [allSNPs,]; 

	ldMatrixAll = mldest(genomatSNPs, K = 4, nc = 1, type = "comp", se=F);
	ldMatrix    = ldMatrixAll [,c(3,4,7)]

	# Create a table with LD SNPs
	msgmsg ("Create a table with LD SNPs...")
	i=1;k=1
	snpsLDTable = NULL
	while (i <= nrow (ldMatrix)) {
		if (ldMatrix[i,"r2"] > R2) {
			snpi = ldMatrix[i,"snpi"]
			snpj = ldMatrix[i,"snpj"]
			r2   = ldMatrix[i,"r2"]
			msgmsg (sprintf ("SNPs in high LD: %s and %s with R2=%s", snpi, snpj, r2))
			df          = data.frame (LD_SNP=paste0 ("LD_SNP_",k),SNP1=snpi, SNP2=snpj, R2=r2, stringsAsFactors=F)
			snpsLDTable = rbind (snpsLDTable, df)
			k = k+1
		}
		i = i+1
	}

	if (is.null (snpsLDTable))
		params$SNPsHighLDFile <<- NULL
	else {
		SNPsHighLDFile = "out/out-SNPsHighLD-Table.csv"
		write.csv (snpsLDTable, SNPsHighLDFile, quote=F, row.names=F)
		params$SNPsHighLDFile <<- SNPsHighLDFile
	}
	return (snpsLDTable)
}
	
#-------------------------------------------------------------
# Move output files to specific directories
#-------------------------------------------------------------
moveOutFiles <- function (outputDir, reportDir, traitConfigFile, phenotypeFile, params) {
	system (sprintf ("mv %s/%s/multiGWAS-report.html %s/%s-report.html", params$runningPath, params$trait, params$runningPath, params$trait))
	system (sprintf ("cp %s/tool*csv %s &> /dev/null", outputDir, reportDir))
	system (sprintf ("mv %s/out*pdf %s > /dev/null 2>&1", outputDir, reportDir))
	system ("mkdir logs")
	system ("mv *.log* logs > /dev/null 2>&1", ignore.stdout=T, ignore.stderr=T)
	system ("mv *.errors* logs > /dev/null 2>&1", ignore.stdout=T, ignore.stderr=T)
	system ("mv *PCs* logs > /dev/null 2>&1", ignore.stdout=T, ignore.stderr=T)
	system ("mv ../*log* logs > /dev/null 2>&1", ignore.stdout=T, ignore.stderr=T)
	system ("mv ../*errors* logs > /dev/null 2>&1", ignore.stdout=T, ignore.stderr=T)
	system ("mv out-mapping* report/ > /dev/null 2>&1", ignore.stdout=T, ignore.stderr=T)

	system (sprintf ("rm ../%s > /dev/null 2>&1", phenotypeFile))
	system (sprintf ("rm ../%s > /dev/null 2>&1", traitConfigFile))
	#system ("rm -rf out & > /dev/null 2>&1")

}

#-------------------------------------------------------------
# First convert genotype to GwasPoly format
# Filters the genotype by different quality control filters
# Read and check files, sample samples
#-------------------------------------------------------------
filterGenoPhenoMapFiles <- function (genotypeFile, phenotypeFile, cacheDir) {
	# If non-model organisms sort and select the number of chromosomes to show
	genoFileNonModel = filterByNonModelOrganisms (genotypeFile, params)

	# Filter by valid markers (alleles lenght, call rate SNPs and Samples)
	genoFileValid = filterByValidMarkers (genoFileNonModel, params$ploidy, cacheDir)

	if (params$filtering==FALSE) {
		msg ("Without filters")
		genoFileFiltered = genoFileValid
	}else {
		msg ("Using filters...")
		genoFileFiltered  = filterByMissingMarkersAndSamples (genoFileValid, params$GENO, params$MIND) 
	}

	# Filter by common markers, samples and remove duplicated and NA phenos
	common = filterByCommonMarkersSamples (genoFileFiltered, phenotypeFile)
	trait  = common$trait

	# Remove no polymorohic markers (MAF > 0.0) and Write chromosome info (map.tbl)
	msgmsg ("Filtering no polymorohic markers (MAF > THRESHOLD) ...") 
	maf    = filterByMAF (common$genotypeFile, params, params$thresholdMAF, cacheDir)

	# Print info trait, N samples
	msgmsg("Evaluating following trait: ", trait) 

	# Add new parameters to params
	params$genotypeFile    <<- maf$genotypeFile
	params$genotypeNumFile <<- maf$genotypeNumFile
	params$phenotypeFile   <<- common$phenotypeFile
	params$mapFile         <<- common$mapFile
	params$trait           <<- trait

	return (list (genotypeFile=params$genotypeFile, phenotypeFile=params$phenotypeFile, 
				  mapFile=params$mapFile, trait=params$trait, genotypeNumFile=params$genotypeNumFile))
}

#-------------------------------------------------------------
# For non-model organisms:
# Sort by chromosome lenght and select N first chromosomes 
#-------------------------------------------------------------
filterByNonModelOrganisms <- function (genotypeFile, params) {
	if (params$nonModelOrganism == FALSE) 
		return (genotypeFile)

	# Genotype is for non-model organism
	# Function for get chromosome limites

	limits <- function (x) {
		start = head (x,n=1)
		end   = tail (x,n=1)
		return (end-start+1)
	}

	data         = read.csv (genotypeFile, check.names=F)
	data         = arrange (data, Chrom, Position)
	grouped      = group_by (data, Chrom)
	sortedChroms = grouped %>% summarize(SIZE=limits(Position)) %>% arrange (-SIZE)
	nLarger      = head (sortedChroms, n=params$numberOfChromosomes)
	selected     = data %>% filter (Chrom %in% nLarger$Chrom)
	dataSorted   = selected %>% arrange (match(Chrom, nLarger$Chrom)) 

	## Change chromosome names
	nLarger = as.data.frame (nLarger)
	chrs    = as.character (nLarger$Chrom)
	dfContigs = NULL
	for (i in 1:length (chrs)) {
		contigName = paste0 ("contig",i)
		dataSorted$Chrom [dataSorted$Chrom %in% chrs[i]] = contigName
		dfContigs = rbind (dfContigs, data.frame (CHROMOSOME= chrs[i], CONTIG=contigName))
	}

	write.csv (dfContigs, "out-mapping-chromosome-names-table.csv", quote=F, row.names=F)

	nonModelFile = addLabel (genotypeFile,"NonMODEL")
	write.csv (dataSorted, nonModelFile, quote=F, row.names=F)
	return (nonModelFile)
}

#-------------------------------------------------------------
# Convert genotype to tool formats, and add to params
#-------------------------------------------------------------
convertPhenotypeToToolFormats <- function (phenotypeFile) {
	# Create PLINK geno/pheno (For SHEsis and TASSEL)
	if (grepl("plink", params$tool) || grepl("tassel",params$tools)) {
		msgmsg ("Converting phenotype to PLINK format...")
		plinkPhenotype = gwaspolyToPlinkPhenotype  (phenotypeFile) 
		params$plinkPhenotypeFile =  plinkPhenotype
	}

	if (grepl("tassel", params$tool)) {
		msgmsg ("Converting phenotype to TASSEL format (.vcf)...")
		tasselPhenotype = gwaspolyToTasselPhenotype (phenotypeFile) 
		params$tasselPhenotypeFile =  tasselPhenotype
	}

	if (grepl("gapit", params$tool)) {
		msgmsg ("Converting geno/pheno to GAPIT format (.vcf)...")
		gapitPhenotype = gwaspolyToGapitPhenotype (phenotypeFile)
		params$gapitPhenotypeFile =  gapitPhenotype
	}
	return (params)
}

#-------------------------------------------------------------
# Returns filename from cache if it exists, otherwise returns NULL
#-------------------------------------------------------------
getFileFromCache <- function (filename, cacheDir) {
	cacheFile = sprintf ("%s/%s", cacheDir, basename (filename))
	if (file.exists (cacheFile)) {
		msgmsg ("Loading file from cache: ", cacheFile)
		return (cacheFile)
	}
	return ("")
}
#-------------------------------------------------------------
# Convert genotype to tool formats, and add to params
#-------------------------------------------------------------
convertGenotypeToToolFormats <- function (genotypeFile, cacheDir) {
	msg ("Converting genotype to tools format...")
	# Create PLINK geno/pheno (For SHEsis and TASSEL)
	if (grepl("plink", params$tool) || grepl("tassel",params$tools)) {
		plinkGenotype = gsub (".csv", "-PLINK", genotypeFile)
		pedFile       = paste0 (plinkGenotype, ".ped")
		if (getFileFromCache (pedFile, cacheDir) == "") {
			msgmsg ("Converting genotype to PLINK format...")
			plinkGenotype = gwaspolyToPlinkGenotype (genotypeFile)
		}

		#if (file.exists (paste0(plinkGenotype, ".ped"))==FALSE) {
		#	msgmsg ("Converting genotype to PLINK format...")
		#	plinkGenotype = gwaspolyToPlinkGenotype (genotypeFile)
		#}
		params$plinkGenotypeFile <<- plinkGenotype
	}

	if (grepl("tassel", params$tool)) {
		tasselGenotype = gsub (".csv", "-TASSEL.vcf", genotypeFile)
		if (file.exists (tasselGenotype)==FALSE){
			msgmsg ("Converting genotype to TASSEL format (.vcf)...")
			tasselGenotype  = plinkToVCFFormat (plinkGenotype)
		}
		params$tasselGenotypeFile <<- tasselGenotype
	}

	if (grepl("gapit", params$tool)) {
		msgmsg ("Converting genotype to GAPIT format...")
		gapit = gwaspolyToGapitGenotype (genotypeFile, params$geneAction, FILES=T, cacheDir)
		params$gapitGenotypeFile  <<- gapit$geno
		params$gapitMapFile       <<- gapit$map
	}
	return (params)
}

#-------------------------------------------------------------
# Check valid markers: 
#   alleles length==ploidy, proportion of NAs in marker and samples (call rate),
#   chromosomes names as number, and sorted by chromosome and position
#-------------------------------------------------------------
filterByValidMarkers <- function (genotypeFile, ploidy, cacheDir) {
	msgmsg ("Checking valid markers...")
	validGenoFile = addLabel (genotypeFile, "VALIDSNPs")
	if (DEBUG) {
		msgmsg ("Loading valid geno file from cache...")
		genoCache = sprintf ("%s/%s", cacheDir, basename (validGenoFile))
		if (file.exists (genoCache))
			return (genoCache)
		else
			validGenoFile = genoCache
	}
	msgmsg ("Calculating valid markers...")

	geno                  = read.csv (file=genotypeFile, check.names=F)
	allelesGeno           = as.matrix(geno[,-(1:3)])
	allelesNames          = colnames (allelesGeno)
	map                   = geno [,1:3]

	sampleNames           = colnames (geno[,-(1:3)])
	markerNames           = geno [1,]
	rownames(allelesGeno) = geno [,1]

	#--------------------------------------------
	# Check and convert chromosome text names to numbers using factors
	#--------------------------------------------
	##anyNonNumericChrom <- function (chrs) {
	##	suppressWarnings (any (is.na (as.numeric (chrs))))
	##}
 
	##chrs = as.character (map [,2])
	##if  (anyNonNumericChrom (chrs)==TRUE) {
	##	msgmsg ("!!!Mapping chromosome names to numbers (see 'out-mapped-chromosome-names.csv') file...")
	##	chrs            = as.factor (chrs)
	##	levels (chrs)   = 1:length (levels (chrs))
	##	write.csv (data.frame (ORIGINAL_CHROM=map[,2], NEW_CHROM=chrs), "out-mapped-chromosome-names.csv", quote=F, row.names=F)
	##	map [,2] = chrs
	##}
	#--------------------------------------------

	#--------------------------------------------
	# Set NAs to alleles with length < ploidy
	#--------------------------------------------
	setNAs <- function (alleles, ploidy) {
		if (is.na (alleles) || nchar (alleles) < ploidy)
			return (NA)
		return (alleles)
	}
	allelesList   = unlist (mclapply (allelesGeno,  setNAs, ploidy, mc.cores=NCORES))
	if (sum (is.na(allelesList)) > length (allelesList)/2) {
		message ("-------------------------------------------------------------------------------------------")
		message ("ERROR: Most alleles are not ploidy ", ploidy, ". Check the ploidy number in the config file!") 
		message ("------------------------------------------------------------------------------------------")
		quit ()
	}
		
	allelesMat    = matrix (allelesList, nrow=nrow(allelesGeno), ncol=ncol(allelesGeno))

	# Remove NA columns 
	NACols        = which (colSums(!is.na (allelesMat))==0) 
	if (length (NACols)>0) {
		allelesMat   = allelesMat [,-NACols] 
		allelesNames = allelesNames [-NACols]
	}
	colnames (allelesMat) = allelesNames

	# Remove NA rows
	NARows        = which (rowSums(!is.na (allelesMat))==0) 
	if (length (NARows)>0) {
		allelesMat = allelesMat [-NARows,] 
		map        = map [-NARows,]
	}

	validGeno     = cbind (map, allelesMat)
	write.csv (validGeno, validGenoFile, quote=F, row.names=F)
	#--------------------------------------------

	return (validGenoFile)
}
#-------------------------------------------------------------
# Convert input genotype to multigwas format (gwaspoly format)
# Checks if VCF, GWASpoly(k-matrix-chrom-pos), k-matrix, and fitPoly
# Phenotype is assume in table format
#-------------------------------------------------------------
convertGenoPhenoToMultigwas <- function (cacheDir) {
	genoFile = sprintf ("%s/%s", cacheDir, "genotype.csv")
	if (file.exists (genoFile)){
		msg ("Loading genotype from cache...")
		newGenotypeFile = genoFile
	}else {
		type = params$genotypeFormat
		msg ("Converting genotype format in ", type , " to ", "GWASpoly format...")

		if (type=="gwaspoly"){#Don't do anything, same genotype
			newGenotypeFile = params$genotypeFile
		}else if (type=="matrix") {# Only for tetraploids
			newGenotypeFile = createGwaspolyGenotype (params$genotypeFile, params$mapFile)
		}else if (type=="vcf") {
			newGenotypeFile = convertVCFToACGTByNGSEP (params$genotypeFile) #output: filename.csv
		}else if (type=="fitpoly"){
			newGenotypeFile = convertFitpolyToGwaspolyGenotype (params$genotypeFile, params$mapFile) #output: filename.csv
		}else if (type=="updog"){
			newGenotypeFile = convertUpdogToGwaspolyGenotype (params$genotypeFile, params$mapFile) #output: filename.csv
		}else {
			message ("--------------------------------------------")
			message ("ERROR: Unknown genotype file format: ", type)
			message ("--------------------------------------------")
			quit ()
		}
	}
	# Copy converted geno/pheno to base names into cacheDir
	newGenoFile  = sprintf ("%s/%s", cacheDir, "genotype.csv")
	newPhenoFile = sprintf ("%s/%s", cacheDir, "phenotype.csv")
	file.copy (newGenotypeFile,  newGenoFile)
	file.copy (params$phenotypeFile, newPhenoFile)
	params$genotypeFile  = newGenoFile
	params$phenotypeFile = newPhenoFile

	return (params)
}

#-------------------------------------------------------------
# Get alternate allele from a row of alleles
#-------------------------------------------------------------
get.ref <- function(x) {
	y <- paste(na.omit(x),collapse="")
	ans <- apply(array(bases),1,function(z,y){length(grep(z,y,fixed=T))},y)
	if (sum(ans)>2) {stop("Error in genotype matrix: More than 2 alleles")}
	if (sum(ans)==2) {ref.alt <- bases[which(ans==1)]}
	if (sum(ans)==1) {ref.alt <- c(bases[which(ans==1)],NA)}
	
	return(ref.alt)
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

.qvalue <- function(p) {
        smooth.df = 3
        if (min(p) < 0 || max(p) > 1) {
            print("ERROR: p-values not in valid range.")
            return(0)
        }
        lambda = seq(0, 0.9, 0.05)
        m <- length(p)
        pi0 <- rep(0, length(lambda))
        for (i in 1:length(lambda)) {
            pi0[i] <- mean(p >= lambda[i])/(1 - lambda[i])
        }

        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0 <- predict(spi0, x = max(lambda))$y
        pi0 <- min(pi0, 1)
        if (pi0 <= 0) {
            print("ERROR: The estimated pi0 <= 0. Check that you have valid p-values.")
            return(0)
        }
        u <- order(p)
        qvalue.rank <- function(x) {
            idx <- sort.list(x)
            fc <- factor(x)
            nl <- length(levels(fc))
            bin <- as.integer(fc)
            tbl <- tabulate(bin)
            cs <- cumsum(tbl)
            tbl <- rep(cs, tbl)
            tbl[idx] <- tbl
            return(tbl)
        }
        v <- qvalue.rank(p)
        qvalue <- pi0 * m * p/v
        qvalue[u[m]] <- min(qvalue[u[m]], 1)
        for (i in (m - 1):1) {
            qvalue[u[i]] <- min(qvalue[u[i]], qvalue[u[i + 1]], 
                1)
        }
        return(qvalue)
    }

#-------------------------------------------------------------
# Impute, filter by MAF, unify geno and pheno names
# Only for "ACGT" format (For other formats see GWASpoly main)
#-------------------------------------------------------------
filterByMAF <- function(genotypeFile, params, thresholdMAF, cacheDir) {
	ploidy = params$ploidy

	if (is.null (thresholdMAF)) {
		if (params$filtering == TRUE)
			thresholdMAF = params$MAF
		else
			thresholdMAF = 0.0
	}

	outNum   = ACGTToNumericGenotypeFormat (genotypeFile, ploidy, cacheDir)
	geno     = outNum$geno
	genoNum  = outNum$genoNum
	genoMap  = outNum$map
	nMarkers = nrow (geno)
	nSamples = ncol (geno[,-1:-3])

	# Get only genotypes without mapping
	numericMatrixM = t(genoNum[,-1])
	colnames (numericMatrixM) = genoNum [,1]

	# Check LG Global MAF (AF?)
	msgmsg ("Checking minor allele frecuency, MAF=", thresholdMAF)
	#-------------------------------------------------------------
	mafFunction <- function(x) {
		AF <- mean(x,na.rm=T)/ploidy;
		MAF <- ifelse(AF > 0.5,1-AF,AF)
	}
	MAF         = apply (numericMatrixM,2,mafFunction)
	polymorphic = which (MAF>thresholdMAF)
	genoNum     = genoNum [polymorphic, ]
	genoMap     = genoMap [polymorphic,]

	genoMap = genoMap[order(genoMap$Chrom,genoMap$Position),]
	msgmsg ("Number of polymorphic markers:", nrow (genoMap),"\n")

	numericMatrixM = numericMatrixM[,polymorphic]
	#numericMatrixM <- numericMatrixM[,genoMap$Marker]
	
	missing <- which(is.na(numericMatrixM))
	if (length(missing)>0) {
		msgmsg("Missing marker data imputed with population mode...")
		numericMatrixM <- apply(numericMatrixM,2,impute.mode)
	}
	# Write geno MAF
	rownames (geno) = geno [,1]
	genoMAF         = geno [colnames(numericMatrixM),]
	genoMAFFile     = addLabel (genotypeFile, "MAF")
	write.csv (genoMAF, genoMAFFile, row.names=F)

	# Write numeric geno
	rownames (genoNum) = genoNum [,1]
	genotypeNumFile    = addLabel (genotypeFile, "MAF-NUM")  
	genotypeMapFile    = addLabel (genotypeFile, "MAF-MAP")  
	write.csv (genoNum, genotypeNumFile, row.names=F)
	write.csv (genoMap, genotypeMapFile, row.names=F)
							  
	return (list (genotypeFile=genoMAFFile, geno=genoMAF, genotypeNumFile=genotypeNumFile, 
				  nMarkers=nMarkers, nSamples=nSamples))
}

#-------------------------------------------------------------
# Add label to filename and new extension (optional)
#-------------------------------------------------------------
addLabelExt <- function (filename, label, newExt=NULL)  {
	nameext = strsplit (filename, split="[.]")
	name    = nameext [[1]][1] 
	if (is.null (newExt))
		ext     = nameext [[1]][2] 
	else
		ext     = newExt
	newName = paste0 (nameext [[1]][1], "-", label, ".", ext )
	return (newName)
}
#----------------------------------------------------------
# Create dir, if it exists the it is renamed old-XXX
#----------------------------------------------------------
createDir <- function (newDir) {
	checkOldDir <- function (newDir) {
		name  = basename (newDir)
		path  = dirname  (newDir)
		if (dir.exists (newDir) == T) {
			oldDir = sprintf ("%s/old-%s", path, name)
			if (dir.exists (oldDir) == T) {
				checkOldDir (oldDir)
			}
			file.rename (newDir, oldDir)
		}
	}

	checkOldDir (newDir)
	system (sprintf ("mkdir %s", newDir))
}

#--------------------------------------------------------
# Print parameters configuration file
#--------------------------------------------------------
printParameters <- function (params) {
	# Print params file
	msgmsg ("-----------------------------------------------------------")
	msgmsg ("Summary of configuration parameters:")
	msgmsg ("-----------------------------------------------------------")
	for (i in 1:length (params)) 
		msgmsg (sprintf ("%-18s : %s", names (params[i]), 
			if (is.null (params [i][[1]])) "NULL" else params [i][[1]]    ))
	msgmsg ("-----------------------------------------------------------")

}

#-------------------------------------------------------------
# Call main 
#-------------------------------------------------------------
main ()
quit()

tryCatch ({main ()}, 
	error=function (e){
		printLOGS(e)
	})
quit()
withCallingHandlers (
	main (), 
	warning = function (w) { 
		message (geterrmessage ())
	},
	error = function (e) { 
		print (sys.calls());quit()
			msg = geterrmessage ()
			msgBox = paste0 ("\n-----------------------------------------------------------------------------",
							 "\n", msg, "\n",
							 "Check errors.log file\n",
							 "-----------------------------------------------------------------------------\n")
			message (msgBox)
			errorsFile = file (paste0 (OUTDIR,"/errors.log"), open="w")
			sink (errorsFile, type="message")
			if (!grepl ("MG Error:", msg))
				message (msgBox, head (paste0(sys.calls()[-1], "\n\n"),-2))
			quit ()
	}
)
