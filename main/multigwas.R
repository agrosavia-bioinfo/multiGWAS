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
# Used to run in parallel the other functions
#-------------------------------------------------------------
runGWASTools <- function () {
	runOneTool <- function (tool, params) {
		if      (tool=="gwaspoly") runToolGwaspoly (params)
		else if (tool=="shesis")   runToolShesis (params)
		else if (tool=="gapit")    runToolGapit (params)
		else if (tool=="tassel")   runToolTassel (params)
		else if (tool=="plink")    runToolPlink (params)
		else                       stop ("Tool not supported")
	}

	# A string containing the names of the tools to run (e.g. "GWASpoly SHEsis PLINK TASSEL")
	msg ("Preparing to execute in parallel the GWAS tools:") 
	tools = strsplit(params$tools ,split=" ")[[1]]
	for (tool in tools) 
		msgmsg ("Running ", tool)

	#listOfResultsFile     = mclapply (tools, runOneTool, params, mc.cores=NCORES, mc.silent=SILENT)
	listOfResultsFile     = mclapply (tools, runOneTool, params, mc.cores=NCORES, mc.silent=FALSE)
	listOfResultsFileBest = selectBestGeneActionModelAllTools (listOfResultsFile, params$geneAction, params$nBest)
	listOfResultsFileLD   = runLinkageDisequilibriumAnalysis (listOfResultsFileBest, params$nBest, params$genotypeNumFile, params$R2)

	return (listOfResultsFileLD)
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
	params$trait           <<- data$trait
	params$genotypeNumFile <<- data$genotypeNumFile
	params$nBest           <<- as.integer (params$nBest)

	msg ("Converting genotype/phenotype to tools formats...")
	data = convertGenotypeToToolFormats (data$genotypeFile, data$phenotypeFile)
	params$genotypeFile    <<- data$genotypeFile
	params$phenotypeFile   <<- data$phenotypeFile

	# Run the four tools in parallel
	listOfResultsFile = runGWASTools ()

	# Create reports
	msg ("Creating reports (Table, Venn diagrams, Manhattan&QQ plots, SNP profiles)...")
	createReports (listOfResultsFile, params)

	# Move out files to output dir
	msg ("Moving files to output folders...")
	moveOutFiles (params$outputDir, params$reportDir, params)
}

#-------------------------------------------------------------
# Main for multi traits
#-------------------------------------------------------------
main <- function () {
	source (paste0 (HOME, "/main/gwas-lib.R"))             # Module with shesis functions

	msg ("MultiGWAS 1.0")
	msg ("Working dir: ", getwd())
	args = commandArgs(trailingOnly = TRUE)

	if (length (args) < 1) 
		stop (usageInstructions())

	# Check dependencies
	checkDependencies (params)

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

	# Load main
	source (paste0 (HOME, "/main/gwas-preprocessing.R"))      # Module with functions to convert between different genotype formats 
	source (paste0 (HOME, "/main/gwas-reports.R"))            # Module with functions to create summaries: tables and venn diagrams
	source (paste0 (HOME, "/main/gwas-heatmap.R"))            # Module with functions to create heatmaps for shared SNPs
	source (paste0 (HOME, "/main/gwas-gwaspoly.R"))           # Module with gwaspoly functions
	source (paste0 (HOME, "/main/gwas-plink.R"))              # Module with plink functions
	source (paste0 (HOME, "/main/gwas-tassel.R"))             # Module with tassel functions
	source (paste0 (HOME, "/main/gwas-shesis.R"))             # Module with shesis functions
	source (paste0 (HOME, "/main/gwas-gapit.R"))             # Module with shesis functions
	#source (paste0 (HOME, "/main/gwas-parameters.R"))             # Module with shesis functions
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
	if (is.null (params$R2)) params$R2 = 1

	# Change to lower case text parameters
	params$genotypeFormat   = tolower (params$genotypeFormat) 
	params$gwasModel        = tolower (params$gwasModel) 
	params$filtering        = ifelse (tolower (params$filtering)=="true", T, F) 
	params$tools            = tolower (params$tools) 
	params$geneAction       = tolower (params$geneAction) 

	params$correctionMethod   = tolower (params$correctionMethod) 
	if (params$correctionMethod == "bonferroni") params$correctionMethod = "Bonferroni"
	else if (params$correctionMethod == "fdr")   params$correctionMethod = "FDR"
	else stop (paste0 ("MG Error: Unknown correction method: ", params$correctionMethod), call.=T)

	`%notin%` <- Negate(`%in%`)

	# Check possible errors in ploidy 
	if (params$ploidy %notin% c("2", "4"))   stop ("MG Error: Ploidy not supported")

	# Create output dir, check input files, and copy files to output dir
	runningDir   = paste0 ("out-", strsplit (paramsFile, split="[.]") [[1]][1])
	createDir (runningDir)
	if (!file.exists (params$genotypeFile)) 
		stop (sprintf ("MG Error: Genotype file not found: '%s'", params$genotypeFile), call.=T)
	runCommand(sprintf ("cp %s %s", params$genotypeFile, runningDir))
	params$genotypeFile  = basename (params$genotypeFile)

	if (!file.exists (params$phenotypeFile)) 
		stop (sprintf ("MG Error: Phenotype file not found: '%s'", params$phenotypeFile), call.=T)
	runCommand(sprintf ("cp %s %s", params$phenotypeFile, runningDir))
	params$phenotypeFile = basename (params$phenotypeFile)

	if (tolower (params$genotypeFormat) %in% c("kmatrix", "fitpoly", "updog")) {
		if (is.null (params$mapFile) | !file.exists (params$mapFile))      
			stop ("MG Error: Map file not found or not specified in the config file", call.=T)
		file.copy (params$mapFile, runningDir)
		params$mapFile = basename (params$mapFile)
	}
	# Change to the output dir and set global OUTDIR 
	setwd (runningDir)
	params$runningPath = getwd()
	
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

	if (geneAction %in% c("all", "additive", "dominant", "general") | tool %in% c("SHEsis")){
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
runLinkageDisequilibriumAnalysis <- function (listOfResultsFile, nBest, genotypeNumFile, R2) {
	msg ("Linkage disequilibrium analysis...")
	ldTable = createLinkageDisequilibriumTable (listOfResultsFile, nBest, genotypeNumFile, R2) 
	if (is.null (ldTable))
		return (listOfResultsFile)

	i = 1
	for (res in listOfResultsFile) {
		msgmsg ("LD in ", res$tool)
		scoresLD = res$scores
		view (scoresLD)
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
				view (listOfResultsFile [[i]]$scores)
				listOfResultsFile [[i]]["scoresLDFile"] = scoresLDFile
		}
		i = i+1
	}
	return (listOfResultsFile)
}

#-------------------------------------------------------
# Create a table for SNPs in all tools in high LD
#-------------------------------------------------------
createLinkageDisequilibriumTable <- function (listOfResultsFile, nBest, genotypeNumFile, R2) {
	# Join all SNPs in one vector
	allSNPs = NULL
	for (res in listOfResultsFile) {
		toolSNPs = res$scores$Marker [1:nBest] 
		allSNPs  = union (allSNPs, toolSNPs)
	}

	# Run linkage disequilibrium
	genotypeNum = read.csv (genotypeNumFile, row.names=1);
	genomat     = as.matrix (genotypeNum [,-1:-2]); 
	genomatSNPs = genomat [allSNPs,]; 

	ldMatrixAll = mldest(genomatSNPs, K = 4, nc = 1, type = "comp", se=F);
	ldMatrix    = ldMatrixAll [,c(3,4,7)]

	# Create a table with LD SNPs
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
moveOutFiles <- function (outputDir, reportDir, params) 
{
	system (sprintf ("mv %s/%s/multiGWAS-report.html %s/%s-report.html", params$runningPath, params$trait, params$runningPath, params$trait))
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
	maf             = filterByMAF (validGenotypeFile, params, 0.0)
	#maf             = filterByMAF (validGenotypeFile, params, params$thresholdMAF)
	genotypeFile    = maf$genotypeFile
	genotypeNumFile = maf$genotypeNumFile

	if (params$filtering==FALSE) {
		msg ("Without filters")
	}else {
		msg ("Using filters...")
		msgmsg ("Filtering by missing markers and samples...")
		genotypeFile  = filterByMissingMarkersAndSamples (genotypeFile, params$GENO, params$MIND) 
	}

	# Filter by common markers, samples and remove duplicated and NA phenos
	msgmsg ("Filtering by common markers and samples...")
	common          = filterByCommonMarkersSamples (genotypeFile, phenotypeFile, genotypeNumFile)
	genotypeFile    = common$genotypeFile
	genotypeNumFile = common$genotypeNumFile
	#genotypeNumFile = ACGTToNumericGenotypeFormat (genotypeFile, params$ploidy)
	phenotypeFile   = common$phenotypeFile
	trait           = common$trait

	# Print info trait, N samples
	msgmsg("Evaluating following trait: ", trait) 
	nSamples = ncol (common$genotype)
	msgmsg ("N =",nSamples,"individuals with phenotypic and genotypic information \n")

	return (list (genotypeFile=genotypeFile, phenotypeFile=phenotypeFile, trait=common$trait, genotypeNumFile=genotypeNumFile))
}

#-------------------------------------------------------------
# Convert genotype to tool formats
#-------------------------------------------------------------
convertGenotypeToToolFormats <- function (genotypeFile, phenotypeFile) {
	# Create PLINK geno/pheno (For SHEsis and TASSEL)
	if (grepl("plink", params$tool) | grepl("tassel",params$tools)) {
		msgmsg ("Converting phenotype to PLINK format (.ped, .map, .bim, .fam, .bed)...")
		plinkPhenotype = gwaspolyToPlinkPhenotype  (phenotypeFile) 
		plinkGenotype  = gwaspolyToPlinkGenotype (genotypeFile)
		params$plinkPhenotypeFile <<- plinkPhenotype
		params$plinkGenotypeFile <<- plinkGenotype
	}

	if (grepl("tassel", params$tool)) {
		msgmsg ("Converting phenotype to TASSEL format (.vcf)...")
		tasselPhenotype = gwaspolyToTasselPhenotype (phenotypeFile) 
		tasselGenotype  = plinkToVCFFormat (plinkGenotype)
		params$tasselPhenotypeFile <<- tasselPhenotype
		params$tasselGenotypeFile <<- tasselGenotype
	}

	if (grepl("gapit", params$tool)) {
		msgmsg ("Converting geno/pheno to GAPIT format (.vcf)...")
		gapit = gwaspolyToGapitFormat (genotypeFile, phenotypeFile, params$geneAction, FILES=T)
		params$gapitGenotypeFile  <<- gapit$geno
		params$gapitPhenotypeFile <<- gapit$pheno
		params$gapitMapFile       <<- gapit$map
	}

	return (list (genotypeFile=genotypeFile, phenotypeFile=phenotypeFile))
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
# Only for "ACGT" format (For other formats see GWASpoly main)
#-------------------------------------------------------------
filterByMAF <- function(genotypeFile, params, thresholdMAF) {
	ploidy = params$ploidy

	if (is.null (thresholdMAF)) {
		if (params$filtering == TRUE)
			thresholdMAF = params$MAF
		else
			thresholdMAF = 0.0
	}

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
	genotypeNumFile        = addLabel (genotypeFile, "MAF-NUM")  
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
