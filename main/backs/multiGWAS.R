#!/usr/bin/Rscript
DEBUG = F; SILENT=T
message ("DEBUG = ", DEBUG, "\n")

options (width=300, error=traceback)
if (DEBUG) {SILENT <- FALSE;options (warn=2)}

# Get enviornment GWAS variable 
HOME   <<- Sys.getenv ("MULTIGWAS_HOME")
OUTDIR <- getwd()
params <- list ()

# INFO  : Tool for running GWAS integratind four GWAS tools: 
#         GWASpoly and SHEsis for polyploids species, and Plink and Tassel for diploids.
# AUTHOR: Luis Garreta (lgarreta@agrosavia.co) 
# DATE  : 12/feb/2020
# LOGS  :   
	# r1.8: Add linkage disequilibrium analysis
	# r1.5: Using VCF files
	# r1.3: Added column gene action model to tables of results by tool

#-------------------------------------------------------------
# Return string with usage instructions
#-------------------------------------------------------------
usageInstructions <- function () {
	USAGE="USAGE: multiGWAS <config file>"
	return (USAGE)
}

#-------------------------------------------------------------
# Main for multi traits
#-------------------------------------------------------------
main <- function () {
	source (paste0 (HOME, "/sources/gwas-lib.R"))             # Module with shesis functions

	msg ("MultiGWAS 1.0")
	msg ("Working dir: ", getwd())
	args = commandArgs(trailingOnly = TRUE)

	if (length (args) < 1) 
		stop (usageInstructions())

	# Process config file and run multiGWAS
	msg ("Processing config file...")
	initGlobalEnvironment ()

	configFile = args [1]
	if (file.exists (configFile)==F) 
		stop ("Configuration file not found")

	# Read and check config file arguments
	msg ("Reading configuration file...")
	params     = readCheckConfigParameters (configFile)

	mainDir = getwd ()
	for (traitParamsFile in params$traitConfigList) {
		setwd (mainDir)
		mainSingleTrait (traitParamsFile)
	}
}

#-------------------------------------------------------------
# Define global variables, load packages, and load sources
#-------------------------------------------------------------
initGlobalEnvironment <- function () 
{
	message (">>>> Loading libraries and setting globals...")

	.libPaths (paste0(HOME, "/opt/Rlibs"))

	# Load packages
	suppressMessages (library (GWASpoly)) #
	suppressMessages (library (parallel)) #
	suppressMessages (library (yaml))  # For read config file

	NCORES <<- ifelse (DEBUG==T, 1, detectCores ())

	# New class for gwaspoly
	setClass ("GWASpolyStruct", slots=c(params="list"), contains="GWASpoly.K")

	# Load sources
	source (paste0 (HOME, "/sources/gwas-preprocessing.R"))      # Module with functions to convert between different genotype formats 
	source (paste0 (HOME, "/sources/gwas-reports.R"))            # Module with functions to create summaries: tables and venn diagrams
	source (paste0 (HOME, "/sources/gwas-heatmap.R"))            # Module with functions to create heatmaps for shared SNPs
	source (paste0 (HOME, "/sources/gwas-gwaspoly.R"))           # Module with gwaspoly functions
	source (paste0 (HOME, "/sources/gwas-plink.R"))              # Module with plink functions
	source (paste0 (HOME, "/sources/gwas-tassel.R"))             # Module with tassel functions
	source (paste0 (HOME, "/sources/gwas-shesis.R"))             # Module with shesis functions
	source (paste0 (HOME, "/sources/gwas-gapit.R"))             # Module with shesis functions
}

#-------------------------------------------------------------
# Create individual phenotype files from a multitrait phenotype file
#-------------------------------------------------------------
createTraitConfigFiles <- function (phenotype, traitName, configFile, params) 
{
	# Create phenotype file for trait
	phenotypeTrait     = phenotype [, c(colnames(phenotype)[1], traitName)]
	phenotypeTraitPath = paste0 (traitName,".csv")
	write.csv (phenotypeTrait, phenotypeTraitPath, quote=F, row.names=F)

	# Create config file for trait
	params$phenotypeFile = basename (phenotypeTraitPath)
	configLines = c()
	for (n in names (params))
		configLines = c (configLines, paste0 (n, " : ", params[n]))

	configLines      = c (configLines, paste0 ("outputDir : ", traitName))
	traitParamsFile  = paste0 (traitName, ".config")
	writeLines (configLines, traitParamsFile)

	return (traitParamsFile)
}

#-------------------------------------------------------------
# Main for a single trait
#-------------------------------------------------------------
mainSingleTrait <- function (traitParamsFile) {
	# Copy files to working dirs and create output dirs
	params   <<- getTraitConfig (traitParamsFile)

	# Read, filter, and check phenotype and genotype
	msg ("Preprocessing genomic data (Filtering and Formating data)...")

	data <- genoPhenoMapProcessing (params$genotypeFile, params$genotypeFormat,
									params$phenotypeFile, params$mapFile)

	params$genotypeFile    = data$genotypeFile
	params$phenotypeFile   = data$phenotypeFile
	params$trait           = data$trait
	params$genotypeNumFile = data$genotypeNumFile
	params$nBest           = as.integer (params$nBest)

	# Run the four tools in parallel
	results           = runGWASTools (params)
	listOfResultsFile = results$BestList
	SNPsLDTable       = results$SNPsLDTable

	# Create reports
	msg ("Creating reports (Table, Venn diagrams, Manhattan&QQ plots, SNP profiles)...")
	createReports (listOfResultsFile, SNPsLDTable, params)

	# Move out files to output dir
	msg ("Moving files to output folders...")
	moveOutFiles (params$outputDir, params$reportDir)
}
#-------------------------------------------------------------
# Read configuration parameters and set working directories
#   - Create dir, if it exists, it is renamed as old-XXXX
#   - Copy and make links of geno/pheno to output dirs
#-------------------------------------------------------------
getTraitConfig <- function (traitParamsFile) {
	msg ("Processing params file: ", traitParamsFile)
	#params     = config::get (file=traitParamsFile, config="advanced") 
	params     = yaml.load_file (traitParamsFile, merge.precedence="order") 
	traitDir   = params$outputDir
	params$traitDir = traitDir
	outDir     = paste0 (traitDir, "/out")


	# Copy files two times to both trait and out dir
	outDirs    = c(traitDir, outDir)
	for (dir in outDirs) {
		createDir (dir)
		file.copy (traitParamsFile, dir)
		file.copy (params$genotypeFile, dir )
		file.copy (params$phenotypeFile, dir)
		if (params$genotypeFormat %in% c("kmatrix", "fitpoly", "updog"))
			runCommand (sprintf ("cp %s %s", params$mapFile, dir))
	}
	# Change to the working dir and set dirs in params
	setwd (traitDir)
	params$outputDir <- "out/"
	params$reportDir <- "report/"

	return (params)
}

#-------------------------------------------------------------
# Get params from config file and define models according to ploidy
#-------------------------------------------------------------
readCheckConfigParameters <- function (paramsFile) {
	params = tryCatch (yaml.load_file (paramsFile), 
					   error=function (cond){
						   stop ("Some errors in configuration file. Check for valid or repeated names!!!")
					   })

	params$paramsFilename = paramsFile

	# Set default values if not set
	if (is.null (params$geneAction)) params$geneAction = "additive"
	if (is.null (params$traitType)) params$traitType = "quantitative"

	# Change to lower case text parameters
	params$genotypeFormat   = tolower (params$genotypeFormat) 
	params$gwasModel   = tolower (params$gwasModel) 
	params$filtering   = ifelse (tolower (params$filtering)=="true", T, F) 
	params$tools       = tolower (params$tools) 
	params$geneAction  = tolower (params$geneAction) 
	params$correctionMethod   = tolower (params$correctionMethod) 
	if (params$correctionMethod == "bonferroni") params$correctionMethod = "Bonferroni"
	else if (params$correctionMethod == "fdr")   params$correctionMethod = "FDR"
	else stop (paste0 ("MG Error: Unknown correction method: ", params$correctionMethod), call.=T)

	`%notin%` <- Negate(`%in%`)

	# Check possible errors in ploidy 
	if (params$ploidy %notin% c("2", "4"))   stop ("MG Error: Ploidy not supported")

	# Create output dir, check input files, and copy files to output dir
	outDir   = paste0 ("out-", strsplit (paramsFile, split="[.]") [[1]][1])
	createDir (outDir)
	if (!file.exists (params$genotypeFile)) 
		stop (sprintf ("MG Error: Genotype file not found: '%s'", params$genotypeFile), call.=T)
	runCommand(sprintf ("cp %s %s", params$genotypeFile, outDir))
	params$genotypeFile  = basename (params$genotypeFile)

	if (!file.exists (params$phenotypeFile)) 
		stop (sprintf ("MG Error: Phenotype file not found: '%s'", params$phenotypeFile), call.=T)
	runCommand(sprintf ("cp %s %s", params$phenotypeFile, outDir))
	params$phenotypeFile = basename (params$phenotypeFile)

	if (tolower (params$genotypeFormat) %in% c("kmatrix", "fitpoly", "updog")) {
		if (is.null (params$mapFile) | !file.exists (params$mapFile))      
			stop ("MG Error: Map file not found or not specified in the config file", call.=T)
		file.copy (params$mapFile, outDir)
		params$mapFile = basename (params$mapFile)
	}
	# Change to the output dir and set global OUTDIR 
	setwd (outDir)
	#OUTDIR <<- paste0 (getwd(), "/", outDir)
	
	# Create params files for each trait
	phenotype = read.csv (params$phenotypeFile, check.names=F)
	traitList = colnames (phenotype)[-1]
	traitConfigList = c()
	for (traitName in traitList) {
		params$trait = traitName
		traitConfig     = createTraitConfigFiles (phenotype, traitName, paramsFile, params)
		traitConfigList = c (traitConfigList, traitConfig)
	}

	params$traitConfigList = traitConfigList 

	# Print params file
	msgmsg ("-----------------------------------------------------------")
	msgmsg ("Summary of configuration parameters:")
	msgmsg ("-----------------------------------------------------------")
	for (i in 1:length (params)) 
		msgmsg (sprintf ("%-18s : %s", names (params[i]), 
			if (is.null (params [i][[1]])) "NULL" else params [i][[1]]    ))
	msgmsg ("-----------------------------------------------------------")



	return (params)
}

#-------------------------------------------------------------
# Used to run in parallel the other functions
#-------------------------------------------------------------
runGWASTools <- function (params) {
	runOneTool <- function (tool, params) {
		if      (tool=="gwaspoly") runToolGwaspoly (params)
		else if (tool=="plink")    runToolPlink (params)
		else if (tool=="shesis")   runToolShesis (params)
		else if (tool=="tassel")   runToolTassel (params)
		else if (tool=="gapit")    runToolGapit (params)
		else                       stop ("Tool not supported")
	}

	# A string containing the names of the tools to run (e.g. "GWASpoly SHEsis PLINK TASSEL")
	params$tools       = strsplit(tolower (params$tools) ,split=" ")[[1]]
	msg ("Preparing to execute in parallel the GWAS tools:") 
	for (i in 1:length(params$tools)) 
		msgmsg ("Running ", params$tools [i])

	listOfResultsFile     = mclapply (params$tools, runOneTool, params, mc.cores=NCORES, mc.silent=SILENT)

	listOfResultsFileBest = selectBestGeneActionModelAllTools (listOfResultsFile, params$geneAction, params$nBest)

	SNPsLDTable           = getSNPsHighLDAllTools (listOfResultsFileBest, params)

	#return (listOfResultsFileLD)
	return (list (BestList=listOfResultsFileBest, SNPsLDTable=SNPsLDTable))
}

#-----------------------------------------------------------
# Select best N SNPs from multiple action models (for GWASpoly and TASSEL)
# Uses three criteria: best GC, best replicability, and best significants
# PLINK also can produce info of more action models using options
#-----------------------------------------------------------
selectBestGeneActionModelAllTools <- function (listOfResultsFile, geneAction, nBest) {
	msg ("Selecting best gene action model for all tools...")
	i = 1
	for (res in listOfResultsFile) {
		msgmsg ("Best gene action for: ", res$tool)
		scoresFileBest   = addLabel (res$scoresFile, "BEST")
		bestScoresTool   = selectBestGeneActionModelTool (res$scores, nBest, res$tool, geneAction)
		write.table (bestScoresTool, scoresFileBest, sep="\t", quote=F, row.names=F)
		listOfResultsFile [[i]]$scoresFile = scoresFileBest
		listOfResultsFile [[i]]$scores     = bestScoresTool
		i = i + 1
	}
	return (listOfResultsFile)
}


selectBestGeneActionModelTool <- function (scoresTool, nBest, tool, geneAction) {
	# Select main columns

	dataSNPs = scoresTool [,c("MODEL", "GC", "Marker", "CHR", "POS", "P", "SCORE", "THRESHOLD", "DIFF")]; 

	if (geneAction=="all") {
		bestScoresTool = distinct (arrange (dataSNPs, -DIFF), Marker, .keep_all=T)
		return (bestScoresTool)
	}

	if (geneAction %in% c("additive", "dominant", "general") | tool %in% c("SHEsis"))
		return (dataSNPs)


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

#-------------------------------------------------------
# Create a table for SNPs in all tools in high LD
#-------------------------------------------------------
getSNPsHighLDAllTools <- function (listOfResultsFile, params) {
	msg ("Analyzing linkage disequilibrium for SNPs in each tool...")
	# Create table with all scores
	snpsLDTable = NULL
	for (res in listOfResultsFile) {
		scoresTool = data.frame (TOOL=res$tool, res$scores [,1:9]) # 1:9 are the first common columns for scores
		snpsLD     = getSNPsHighLDTool (params$genotypeNumFile, scoresTool, 0.9, params$nBest, res$tool)
		snpsLDTable = rbind (snpsLDTable, snpsLD)
	}

	return (snpsLDTable)
}

#-----------------------------------------------------------------------
# Return a vector of SNPs in high LD
#-----------------------------------------------------------------------
getSNPsHighLDTool <- function (genoNumFile, scores, maxLD, maxBest, tool) {
	if (!exists ("genotypeNum")) 
		genotypeNum    <<- read.csv (genoNumFile, row.names=1);

	genomat = as.matrix (genotypeNum [,-1:-2]); 

	# Get Top SNPs from score file
	N = if (length (scores$Marker) > 2*maxBest) 2*maxBest else length (scores$Marker) 
	snpList = as.character (scores$Marker [1:N]) 

	# Get genotypes for SNPs and calculate LD matrix (r2)
	genomatSNPs = genomat [snpList,]; 
	ldMatrixAll = mldest(genomatSNPs, K = 4, nc = NCORES, type = "comp", se=F);
	ldMatrix    = ldMatrixAll [,c(3,4,7)]

	# Create a table with LD SNPs
	i=1
	snpsLD = NULL
	msgmsg (tool)
	while (i <= nrow (ldMatrix)) {
		if (ldMatrix[i,"r2"] > maxLD) {
			msgmsg (sprintf ("SNPs of %s in high LD: %s and %s with R2=%s", tool,  ldMatrix[i,"snpi"], ldMatrix[i,"snpj"], ldMatrix[i,"r2"]))
			df     = data.frame (TOOL=tool,SNP1=ldMatrix [i, "snpi"], SNP2=ldMatrix [i, "snpj"], R2=ldMatrix[i,"r2"])
			snpsLD = rbind (snpsLD, df)
			#snpj     = ldMatrix [i, "snpj"]
			#ldMatrix = ldMatrix [ldMatrix$snpi != snpj,]
			#snpsLD   = c(snpsLD, snpj)
		}
		i = i+1
	}
	return (snpsLD)
}

#-------------------------------------------------------
# Analyze and remove for each tool SNPs in high LD
#-------------------------------------------------------
removeLinkageDisequilibriumSNPsTools <- function (listOfResultsFile, params) {
	msg ("Analyzing linkage disequilibrium for SNPs in each tool...")
	# Create table with all scores
	scoresAll = NULL
	for (res in listOfResultsFile) {
		df = data.frame (TOOL=res$tool, res$scores [,1:9]) # 1:9 are the first common columns for scores
		scoresAll = rbind (scoresAll, df)
	}

	# Remove SNPs in LD from table with all scores
	'%ni%' <- Negate('%in%')
	for (res in listOfResultsFile) { 
		scoresTool = filter (scoresAll, TOOL==res$tool)
		snpsLD = matchSNPsByLDSingleTool (params$genotypeNumFile, scoresTool, 0.9, params$nBest, res$tool)
		scoresAll = scoresAll [scoresAll$Marker %in% setdiff (scoresAll$Marker, snpsLD),]
	}

	# Create new versions of score files
	i = 1; 
	for (res in listOfResultsFile) {
		scoresTool   = scoresAll [scoresAll$TOOL %in% res$tool,] 
		scoresFileLD = addLabel (res$scoresFile, "LD")
		write.table (scoresTool, scoresFileLD, sep="\t", quote=F, row.names=F)
		listOfResultsFile [[i]]$scoresFile = scoresFileLD
		i = i + 1
	}

	return (listOfResultsFile)
}

#-------------------------------------------------------------
# Move output files to specific directories
#-------------------------------------------------------------
moveOutFiles <- function (outputDir, reportDir) 
{
	system (sprintf ("cp %s/tool*csv %s &> /dev/null", outputDir, reportDir))
	system (sprintf ("mv %s/out*pdf %s > /dev/null 2>&1", outputDir, reportDir))
	system ("mkdir logs")
	system ("mv *.log* logs > /dev/null 2>&1", ignore.stdout=T, ignore.stderr=T)
	system ("mv *.errors* logs > /dev/null 2>&1", ignore.stdout=T, ignore.stderr=T)
	system ("mv *PCs* logs > /dev/null 2>&1", ignore.stdout=T, ignore.stderr=T)
	system ("mv ../*log* logs > /dev/null 2>&1", ignore.stdout=T, ignore.stderr=T)
	system ("mv ../*errors* logs > /dev/null 2>&1", ignore.stdout=T, ignore.stderr=T)
}

#-------------------------------------------------------------
# Filters the genotype by different quality control filters
# Read and check files, sample samples
# Convert geno an pheno to other tool formats 
#-------------------------------------------------------------
genoPhenoMapProcessing <- function (genotypeFile, genotypeFormat, phenotypeFile, mapFile) {
	# Convert to gwaspoly format from other formats: VCF, GWASpoly, k-matrix, or fitPoly.
	newData       = convertGenotypeToGWASpolyFormat (genotypeFile, phenotypeFile, mapFile, genotypeFormat, params$ploidy, outDir="out/")
	genotypeFile  = newData$geno
	phenotypeFile = newData$pheno
	mapFile       = newData$map

	# Filter by valid markers (alleles lenght, call rate SNPs and Samples)
	msgmsg ("Checking valid markers...")
	validGenotypeFile = filterByValidMarkers (genotypeFile, params$ploidy)

	# Remove no polymorohic markers (MAF > 0.0) and Write chromosome info (map.tbl)
	msgmsg ("Filtering no polymorohic markers (MAF > THRESHOLD) ...") 
	maf             = filterByMAF (validGenotypeFile, params)
	genotypeFile    = maf$genotypeFile
	genotypeNumFile = maf$genotypeNumFile

	if (params$filtering==FALSE) {
		msgmsg ("Without filters")
	}else {
		msgmsg ("Using filters...")
		msgmsg ("Filtering by missing markers and samples...")
		genotypeFile  = filterByMissingMarkersAndSamples (genotypeFile, params$GENO, params$MIND) 
	}

	# Filter by common markers, samples and remove duplicated and NA phenos
	msgmsg ("Selecting common markers and samples names...")
	common          = filterByCommonMarkersSamples (genotypeFile, phenotypeFile, genotypeNumFile)
	genotypeFile    = common$genotypeFile
	genotypeNumFile = ACGTToNumericGenotypeFormat (genotypeFile, params$ploidy)
	phenotypeFile   = common$phenotypeFile
	trait           = common$trait

	# Print info trait, N samples
	msgmsg("Evaluating following trait: ", trait) 
	nSamples = ncol (common$genotype)
	msgmsg ("N =",nSamples,"individuals with phenotypic and genotypic information \n")

	# Create PLINK geno/pheno (For SHEsis and TASSEL)
	msgmsg ("Converting phenotype to PLINK format (.ped, .map, .bim, .fam, .bed)...")
	plinkPhenotype = gwaspolyToPlinkPhenotype  (phenotypeFile) 
	plinkGenotype  = gwaspolyToPlinkGenotype (genotypeFile)
	params$plinkPhenotypeFile <<- plinkPhenotype
	params$plinkGenotypeFile <<- plinkGenotype

	params$tools = strsplit(tolower (params$tools) ,split=" ")[[1]]
	if ("tassel" %in% params$tools) {
		msgmsg ("Converting phenotype to TASSEL format (.vcf)...")
		tasselPhenotype = gwaspolyToTasselPhenotype (phenotypeFile) 
		tasselGenotype  = plinkToVCFFormat (plinkGenotype)
		params$tasselPhenotypeFile <<- tasselPhenotype
		params$tasselGenotypeFile <<- tasselGenotype
	}

	if ("gapit" %in% params$tools) {
		msgmsg ("Converting geno/pheno to GAPIT format (.vcf)...")
		gapit = gwaspolyToGapitFormat (genotypeFile, phenotypeFile, params$geneAction, FILES=T)
		params$gapitGenotypeFile  <<- gapit$geno
		params$gapitPhenotypeFile <<- gapit$pheno
		params$gapitMapFile       <<- gapit$map
	}

	return (list (genotypeFile=genotypeFile, phenotypeFile=phenotypeFile, trait=common$trait, genotypeNumFile=genotypeNumFile))
}

#-------------------------------------------------------------
# Check valid markers: 
#   alleles length==ploidy, proportion of NAs in marker and samples (call rate),
#   chromosomes names as number, and sorted by chromosome and position
#-------------------------------------------------------------
filterByValidMarkers <- function (genotypeFile, ploidy) 
{
	geno                  = read.csv (file=genotypeFile, check.names=F)
	allelesGeno           = as.matrix(geno[,-(1:3)])
	map                   = geno [,1:3]

	sampleNames           = colnames (geno[,-(1:3)])
	markerNames           = geno [1,]
	rownames(allelesGeno) = geno [,1]

	#--------------------------------------------
	# Check and convert chromosome text names to numbers using factos
	#--------------------------------------------
	anyNonNumericChrom <- function (chrs) {
		suppressWarnings (any (is.na (as.numeric (chrs))))
	}
 
	chrs = as.character (map [,2])
	if  (anyNonNumericChrom (chrs)==TRUE) {
		msgmsg ("!!!Mapping chromosome names to numbers (see 'out-mapped-chromosome-names.csv') file...")
		chrs            = as.factor (chrs)
		levels (chrs)   = 1:length (levels (chrs))
		write.csv (data.frame (ORIGINAL_CHROM=map[,2], NEW_CHROM=chrs), "out-mapped-chromosome-names.csv", quote=F, row.names=F)
		map [,2] = chrs
	}
	#--------------------------------------------

	#--------------------------------------------
	# Set NAs to alleles with length < ploidy
	#--------------------------------------------
	setNAs <- function (alleles, ploidy) {
		if (is.na (alleles) | nchar (alleles) < ploidy)
			return (NA)
		return (alleles)
	}
	allelesList   = unlist (mclapply (allelesGeno,  setNAs, ploidy, mc.cores=NCORES))
	allelesMat    = matrix (allelesList, nrow=nrow(allelesGeno), ncol=ncol(allelesGeno))
	colnames (allelesMat) = colnames (allelesGeno)
	validGeno     = cbind (map, allelesMat)
	validGenoFile = addLabel (genotypeFile, "VALIDSNPs")
	write.csv (validGeno, validGenoFile, quote=F, row.names=F)
	#--------------------------------------------

	return (validGenoFile)
}
#-------------------------------------------------------------
# Return the format type of genotype
# Checks if VCF, GWASpoly(k-matrix-chrom-pos), k-matrix, and fitPoly
#-------------------------------------------------------------
convertGenotypeToGWASpolyFormat <- function (genotypeFile, phenotypeFile, mapFile, type, ploidy, outDir) {
	msg ("Converting genotype format in ", type, " to ", "GWASpoly format...")

	if (type=="gwaspoly"){#Don't do anything, same genotype
		genotypeFile = genotypeFile
	}else if (type=="kmatrix") {# Only for tetraploids
		genotypeFile = createGwaspolyGenotype (genotypeFile, mapFile)
	}else if (type=="vcf") {
		genotypeFile = convertVCFToACGTByNGSEP (genotypeFile) #output: filename.csv
	}else if (type=="fitpoly"){
		genotypeFile = convertFitpolyToGwaspolyGenotype (genotypeFile, mapFile) #output: filename.csv
	}else if (type=="updog"){
		genotypeFile = convertUpdogToGwaspolyGenotype (genotypeFile, mapFile) #output: filename.csv
	}else {
		msgmsg ("Error: Unknown genotype file format")
		return (NULL)
	}

	newGenotypeFile  = paste0 (outDir, "genotype.csv")
	newPhenotypeFile = paste0 (outDir, "phenotype.csv")

	file.copy (genotypeFile, newGenotypeFile)
	file.copy (phenotypeFile, newPhenotypeFile)

	return (list (geno=newGenotypeFile, pheno=newPhenotypeFile, map=mapFile))
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
# Only for "ACGT" format (For other formats see GWASpoly sources)
#-------------------------------------------------------------
filterByMAF <- function(genotypeFile, params) {
	ploidy = params$ploidy

	if (params$filtering == TRUE)
		thresholdMAF = params$MAF
	else
		thresholdMAF = 0.0

	outNum       = ACGTToNumericGenotypeFormat (genotypeFile, ploidy, MAP=T)
	geno         = outNum$geno
	numericGeno  = outNum$genoNum
	map          = outNum$map

	# Get only genotypes without mapping
	numericMatrixM = t(numericGeno[,-c(1:3)])

	# Check LG Global MAF (AF?)
	msgmsg ("Checking minor allele frecuency, MAF=", thresholdMAF)
	#-------------------------------------------------------------
	mafFunction <- function(x) {
		AF <- mean(x,na.rm=T)/ploidy;
		MAF <- ifelse(AF > 0.5,1-AF,AF)
	}
	MAF         <- apply (numericMatrixM,2,mafFunction)
	polymorphic <- which (MAF>thresholdMAF)

	map = map[polymorphic,]
	map = map[order(map$Chrom,map$Position),]
	msgmsg ("Number of polymorphic markers:", nrow (map),"\n")

	numericMatrixM = numericMatrixM[,polymorphic]
	#numericMatrixM <- numericMatrixM[,map$Marker]
	
	missing <- which(is.na(numericMatrixM))
	if (length(missing)>0) {
		msgmsg("Missing marker data imputed with population mode...")
		numericMatrixM <- apply(numericMatrixM,2,impute.mode)
	}
	# Write geno MAF
	rownames (geno) = geno [,1]
	genoMAF         = geno [colnames(numericMatrixM),]
	genoMAFFile     = addLabel (genotypeFile, "MAF")
	write.csv (genoMAF, genoMAFFile, quote=F, row.names=F)

	# Write chromosome info 
	#write.table (file="out/map.tbl", map, quote=F, row.names=F, sep="\t")

	# Write numeric geno
	rownames (numericGeno) = numericGeno [,1]
	numericGeno            = numericGeno [polymorphic, ]
	genotypeNumFile        = addLabel (genotypeFile, "NUM")  
	write.csv (numericGeno, genotypeNumFile, quote=F, row.names=F)
							  
	return (list (genotypeFile=genoMAFFile, geno=genoMAF, genotypeNumFile=genotypeNumFile))
}

#----------------------------------------------------------
# Util to print head of data
# Fast head, for debug only
#----------------------------------------------------------
view <- function (data, n=5,m=6) {
	name = paste (deparse (substitute (data)),":  ")
	if (is.null (dim (data))) {
		dimensions = paste (length (data))
		message (name, "(", paste0 (dimensions),")")
		if (length (data) < 6) n = length(data)
		print (data[1:n])
	}else {
		dimensions = paste0 (unlist (dim (data)),sep=c(" x ",""))
		message (name, "(", paste0 (dimensions),")")
		if (nrow (data) < 5) n = nrow(data)
		if (ncol (data) < 6) m = ncol(data)
		print (data[1:n,1:m])
	}
	#write.csv (data, paste0("x-", filename, ".csv"), quote=F, row.names=F)
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

#-------------------------------------------------------------
# Call main 
#-------------------------------------------------------------
main ()
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
