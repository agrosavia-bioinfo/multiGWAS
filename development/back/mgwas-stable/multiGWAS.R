#!/usr/bin/Rscript

# INFO  : Tool for running GWAS using four GWAS tools: GWASpoly and SHEsis for polyploids species, and Plink and Tassel for diploids.
# AUTHOR: Luis Garreta (lgarreta@agrosavia.co) 
# DATA  : 12/feb/2020
# LOGS  :   
	# r0.91: Working with rmarkdown (without QR profile)
	# r0.9: Full functional, but it needs to be improved: own filtes (HWE), SHEsis Full model, corrections and thresholds (FDR, BONF), Manhattan plots
	# r12.0: Fixed preprocessing using read.GWASpoly. Fixed getQTLs for gwaspoly gwas. Full functionality
	# r11.0: Improved naive with unity kinship matrix. Results according to plink and tassel. It needs to show in summary
	# r10.0: Fixed convertion of tetra to diplos. Excelent results with gwaspoly data
	# r8.0: Full summary, reorganized scores tables, checked scores/thresholds, two correction methods: FDR, BONF. 

# Constants
DEBUG              = F
LOAD_DATA          = FALSE
SIGNIFICANCE_LEVEL = 0.05      # Minimun level of significance (alpha) to considerer a significant SNP
MAX_BEST           = 8         # Max number of SNPs of best scored SNPs to show in tables and graphics

args = commandArgs(trailingOnly = TRUE)
#args = c("config-TraitPA-ModelStructure-FiltersFalse-CorrectionBONF.config")
#args = c("in/config-Gota-Naive-filtersNone-impute.config")
args = c("config-TraitTuberShape-ModelNaive-FiltersTrue-CorrectionBONF.config")

USAGE="USAGE: Rscript gwas-polypiline.R <config file>"
if (length (args) != 1) {
	message (USAGE)
	quit()
}
#-------------------------------------------------------------
# Load required packages
#-------------------------------------------------------------
HOME = Sys.getenv ("MULTIGWAS_HOME")
library (GWASpoly) #
library (parallel) #
suppressMessages(library (config))  # For read config file

# New class for gwaspoly
setClass ("GWASpolyStruct", slots=c(params="list"), contains="GWASpoly.K")

options (width=300)
source (paste0 (HOME, "/sources/gwas-preprocessing.R"))     # Module with functions to convert between different genotype formats 
source (paste0 (HOME, "/sources/gwas-summary.R"))           # Module with functions to create summaries: tables and venn diagrams
source (paste0 (HOME, "/sources/scripts/script-gwaspoly.R"))       # Module with gwaspoly functions

#-------------------------------------------------------------
# Global configs
#-------------------------------------------------------------
data = data1 = data2 = data3 = NULL
#-------------------------------------------------------------
# Main
#-------------------------------------------------------------
main <- function (args) 
{
	# Read and check config file arguments
	configFile       <-  args [1]
	config           <- getConfigurationParameters (configFile)
	config$outputDir <- "out/"
	config$reportDir <- "report/"

	# Create output dir "out/" and move to it
	workingDir = initOutputDir (configFile, config$genotypeFile, config$phenotypeFile)

	msg (">>>>>>>>>>>>", config$gwasModel, "<<<<<<<<<<<") 

	# Read, filter, and check phenotype and genotype
	data <- dataPreprocessing (config$genotypeFile, config$phenotypeFile, config)
	config$genotypeFile  =data$genotypeFile
	config$phenotypeFile =data$phenotypeFile

	config$trait = data$trait

	# Run the four tools in parallel
	#runPlinkGwas (config)
	#runGwaspolyGwas (config)
	#runShesisGwas (config)
	mclapply (c("Gwasp", "Plink", "SHEsis", "Tassel"), runGWASTools, config, mc.cores=4)

	# Create reports
	createReports (config$outputDir , config$gwasModel, config$reportDir, nBest=7)

	# Create outputs: tables, figures
	##title = gsub(".*\\config-(.*)\\..*", "\\1", c(configFile))
	##markersSummaryTable ("out/", config$gwasModel, title,  "out/", nBEST=MAX_BEST, significanceLevel=config$signficanceLevel)

	# Move out files to output dir
	moveOutFiles (config$outputDir, config$reportDir)

	# Call to rmarkdown report
	msg ("Creating markdown report...")
	outputFile = paste0 (workingDir, "/multiGWAS-report.html")
	rmarkdown::render (paste0(HOME,"/sources/gwas-markdown.Rmd"), output_file=outputFile,  params=list (workingDir=workingDir))
}
#-------------------------------------------------------------
# Create dir, if it exists, it is renamed as old-XXXX
# Copy and make links of geno/pheno to output dirs
#-------------------------------------------------------------
initOutputDir <- function (configFile, genotypeFile, phenotypeFile) 
{
	outDir   = gsub ("config","test", configFile)
	msg ("Output dir: ", outDir)
	createDir (outDir)
	system (sprintf ("cp %s %s", configFile, outDir))

	runCommand (sprintf ("ln -s %s/%s %s/%s", getwd(), genotypeFile, outDir, genotypeFile))
	runCommand (sprintf ("ln -s %s/%s %s/%s", getwd(), phenotypeFile, outDir, phenotypeFile))

	setwd (outDir)
	system ("mkdir out")
	system (sprintf ("cp %s %s", configFile, "out"))

	# Copy geno/pheno to "out/" dir 
	runCommand (sprintf ("cp -a %s %s", genotypeFile, "out/"))
	runCommand (sprintf ("cp -a %s %s", phenotypeFile, "out/"))

	workingDir = getwd ()
	return (workingDir)

}

#-------------------------------------------------------------
# Move output files to specific directories
#-------------------------------------------------------------
moveOutFiles <- function (outputDir, reportDir) 
{
	msg ("Moving output files to output directories...")
	system (sprintf ("mv %s/out*scores %s", outputDir, reportDir))
	system (sprintf ("mv %s/out*pdf %s", outputDir, reportDir))
	system ("mkdir logs")
	system ("mv *.log* logs")
	system ("mv *.errors* logs")
	system ("mv *PCs* logs")
}

#-------------------------------------------------------------
# Used to run in parallel the other functions
#-------------------------------------------------------------
runGWASTools <- function (tool, config) 
{
	if (tool=="Gwasp")
		runGwaspolyGwas (config)
	else if (tool=="Plink")
		runPlinkGwas (config)
	else if (tool=="SHEsis")
		runShesisGwas (config)
	else if (tool=="Tassel")
		runTasselGwas (config)
	else
		stop ("Tool not supported")
}
#-------------------------------------------------------------
#-------------------------------------------------------------
runGwaspolyGwas <- function (params) 
{
	msg("Running GWASpoly...")

	genotypeFile  = params$genotypeFile
	phenotypeFile = params$phenotypeFile

	# Only for tetra ployds
	ploidy = 4

	#snpModels=testModels = ("general")
	snpModels  = c("general","additive","1-dom", "2-dom")
	testModels = c("general", "additive","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")

	params = append (params, list (snpModels=snpModels, testModels=testModels))

	# Read input genotype and genotype (format: "numeric", "AB", or "ACGT")
	data1 <- initGWAS (phenotypeFile, genotypeFile, ploidy, "ACGT", data1)

	# Control population structure
	data2 = controlPopulationStratification (data1, params$gwasModel, data2)

	# GWAS execution
	data3 <- runGwaspoly (data2, params$gwasModel, params$snpModels, data3)
	showResults (data3, params$testModels, params$trait, params$gwasModel, 
				 params$correctionMethod, params$phenotypeFile, ploidy)
}

#-------------------------------------------------------------
# Plink tool
#-------------------------------------------------------------
runPlinkGwas <- function (params) 
{
	msg ("Running Plink GWAS...")
	model = params$gwasModel

	inGeno   = "out/filtered-plink-genotype"       # Only prefix for Plink
	inPheno  = "out/filtered-plink-phenotype.tbl"  
	outFile  = paste0 ("out/out-Plink-", model)
	outPlink = paste0(outFile,".TRAIT.assoc.linear.adjusted")

	if (model=="Naive") 
		cmm=sprintf ("%s/sources/scripts/script-plink-NaiveModel.sh %s %s %s", HOME,inGeno, inPheno, outFile)
	else if (model=="Structure") {
		# First: kinship filtering
		cmm=sprintf ("%s/sources/scripts/script-kin-kinship.sh %s %s", HOME, inGeno, outFile)
		runCommand (cmm, "log-Kin.log")

		# Second: Structure by PCs 
		inGeno = paste0 (inGeno,"-kinship")
		cmm=sprintf ("%s/sources/scripts/script-plink-FullModel.sh %s %s %s", HOME, inGeno, inPheno, outFile)
	}
	else
		quit (paste ("Type of GWAS:\"", model, "\", not supported"))

	# Data preprocessing
	runCommand (cmm, "log-Plink.log")
	results      <- read.table (file=outPlink,  header=T) 
	pValues      <- results$UNADJ
	if (params$correctionMethod=="FDR") {
		scores       = -log10 (pValues)
		threshold    <- calculateThreshold (level=SIGNIFICANCE_LEVEL, scores=scores, method="FDR")
	}else if (params$correctionMethod=="Bonferroni") {
		# Adjust scores to number of NMISS	(Number of genotype calls considered)
		outPlinkLinear = paste0(outFile,".TRAIT.assoc.linear")
		resultsLinear  = read.table (file=outPlinkLinear,  header=T) 
		resultsLinear  <- resultsLinear [!duplicated (resultsLinear$SNP),]
		rownames (resultsLinear) = resultsLinear$SNP
		nMiss     <- resultsLinear [results$SNP, "NMISS"]
		threshold = -log10 (SIGNIFICANCE_LEVEL/nMiss)
		scores    = -log10 (nMiss*pValues)
		positions =  resultsLinear [results$SNP, "BP"]
	}else
		stop ("Unknown correction method: ", params$correctionMethod)

	resultsAll    <- cbind (results, POS=positions, P=pValues, SCORE=round (scores,6), THRESHOLD=round (threshold,6), DIFF=round (scores-threshold, 6))
	resultsAll    <- resultsAll [order (resultsAll$DIFF, decreasing=T),]
	outScores = paste0 (outFile, ".scores")
	write.table (file=outScores, resultsAll, row.names=F, quote=F, sep="\t")
}
	
#-------------------------------------------------------------
# SHEsis tool
#-------------------------------------------------------------
runShesisGwas <- function (params) 
{
	msg ("Running SHEsis GWAS...")
	model = params$gwasModel

	inGenoPheno  = "out/filtered-shesis-genopheno.tbl"
	inMarkers    = "out/filtered-shesis-markernames.tbl"
	inMarkersPos = "out/filtered-shesis-markernamespos.tbl"
	shesisFile      = paste0 ("out/out-SHEsis-", model)
	scoresFile   = paste0(shesisFile,".scores")

	if (params$gwasModel == "Naive") {
		cmm=sprintf ("%s/sources/scripts/script-shesis-NaiveModel.sh %s %s %s", HOME, inGenoPheno, inMarkers, shesisFile)
		runCommand (cmm, "log-SHEsis.log")
	} else if (params$gwasModel == "Structure") {
	}

	

	# Format data to table with scores and threshold
	results = read.table (file=scoresFile, header=T, sep="\t")
	pValues = results[,"P.value"]
	if (params$correctionMethod=="FDR") {
		scores    <- -log10(p.adjust (pValues, method="fdr"))
		threshold <- calculateThreshold (level=SIGNIFICANCE_LEVEL, scores=scores, method="FDR")
	}else {
		nMiss     = results [results$SNP, "Nonmissing"]
		threshold = -log10 (SIGNIFICANCE_LEVEL/nMiss)
		scores    = -log10 (nMiss*pValues)
	}
	inMarkersPosData <- read.table (file=inMarkersPos, header=F, sep="\t")
	rownames (inMarkersPosData) = inMarkersPosData [,1]
	SNP <- as.character (results$SNP)
	CHR <- inMarkersPosData [SNP, 2]
	POS <- inMarkersPosData [SNP, 3]

	resultsAll <- data.frame (SNP, CHR, POS, P=pValues, SCORE=round (scores,6), THRESHOLD=round (threshold,6), DIFF=round (scores-threshold, 6), results)
	resultsAll <- resultsAll [order (resultsAll$DIFF, decreasing=T),]
	write.table (file=scoresFile, resultsAll, row.names=F, quote=F, sep="\t")
}

#-------------------------------------------------------------
# Run Tassel pipeline (GDM and MLM)
#-------------------------------------------------------------
runTasselGwas <- function (params) 
{
	msg ("Running Tassel GWAS...")
	model = params$gwasModel

	# Parameters for the scripts
	inGenoPED  = "out/filtered-plink-genotype.ped"
	inGenoMAP  = "out/filtered-plink-genotype.map"
	inPhenoTBL = "out/filtered-tassel-phenotype.tbl"
	outPrefix  = paste0("out/out-Tassel-", model)

	if (model=="Naive") {
		cmm=sprintf ("%s/sources/scripts/script-tassel-NaiveModel.sh %s %s %s %s",
					 HOME,inGenoPED, inGenoMAP, inPhenoTBL, outPrefix)
		runCommand (cmm, "log-tassel.log")
		outFile   = list.files("out/", pattern=sprintf("^(.*(%s).*(1).*(txt)[^$]*)$",model), full.names=T)
	}else if (model=="Structure") {
		cmm=sprintf ("%s/sources/scripts/script-tassel-FullModel.sh %s %s %s %s",
					 HOME, inGenoPED, inGenoMAP, inPhenoTBL, outPrefix)
		runCommand (cmm, "log-tassel.log")
		outFile   = list.files("out/", pattern=sprintf("^(.*(%s).*(stats).*(txt)[^$]*)$",model), full.names=T)
	}else 
		quit (paste ("Type of GWAS:\"", model, "\", not implemented"))
	
	# Rename output file
	msg ("Tassel output file: ", outFile)
	outTassel = sprintf ("out/out-Tassel-%s.scores", model)

	# Data preprocessing: Test correction, sort and rename output file
	#ts      = read.table (file=outFile, header=T, sep="\t")
	#tsAdj   = ts %>% mutate (adjP=p.adjust(p,"fdr"),adjPadd=p.adjust(add_p, "fdr"),adjPdom=p.adjust(dom_p, "fdr"))
	#tsMin   = tsAdj %>% rowwise %>% mutate (minP=min(adjP, adjPadd, adjPdom, na.rm=T)) 
	#tsArr   = tsMin %>% arrange (minP)
	#write.table (file=outTassel, tsArr, quote=F, sep="\t", row.names=F)

	# Create tree tables for p, pAdd, and pDom
	results      <- read.table (file=outFile, header=T, sep="\t")
	createTableTassel  <- function (results, model, var) {
		P       <- results [,var]
		scores  <- -log10(results [,var])
		if (params$correctionMethod=="FDR") 
			threshold    <- calculateThreshold (level=SIGNIFICANCE_LEVEL, scores=scores, method="FDR")
		else if  (params$correctionMethod=="Bonferroni") 
			threshold    <- calculateThreshold (level=SIGNIFICANCE_LEVEL, scores=scores, method="Bonferroni")
		else 
			stop (paste0 ("Unknown correction method: ", params$correctionMethod))

		# Compose the results table
		pTable  <-  cbind (Model=model, results [, c("Marker", "Chr", "Pos")], P=P, 
						   SCORE=scores, THRESHOLD=threshold, DIFF=(scores-threshold))
		pTable  <- pTable [order (scores, decreasing=T),]    
		return (pTable)
	}
	# Create tables for each model
	gnrTable <- createTableTassel (results, "Gnr", "p")
	addTable <- createTableTassel (results, "Add", "add_p")
	domTable <- createTableTassel (results, "Dom", "dom_p")

	# Join the tables, order by DIFF, and write it
	scoresTable <- rbind (addTable, domTable, gnrTable)
	scoresTable <- scoresTable [order (scoresTable$DIFF, decreasing=T),]
	scoresTable = scoresTable [!duplicated (scoresTable$Marker),]
	write.table (file=outTassel, scoresTable, quote=F, sep="\t", row.names=F)
}

#-------------------------------------------------------------
# Filters the genotype by different quality control filters
# Read and check files, sample samples
# Convert geno an pheno to other tool formats 
#-------------------------------------------------------------
dataPreprocessing <- function (genotypeFile, phenotypeFile, config) 
{
	msg();msg("Data preprocessing...");msg()

	# Filter by common sample names
	common <- filterByMAFCommonNames (genotypeFile, phenotypeFile)
	#common = filterByCommonNames (genotypeFile, phenotypeFile)
	genotypeFile  = common$genotypeFile
	phenotypeFile = common$phenotypeFile
	gwaspolyData  = common$data

	if (config$filtering == F) {
		# check if paths are asigned to geno/pheno vars
		msg (    "Without filters")
		runCommand (sprintf ("ln -s %s %s", basename (genotypeFile), "out/filtered-gwasp4-genotype.tbl"))
		runCommand (sprintf ("ln -s %s %s", basename (phenotypeFile), "out/filtered-gwasp4-phenotype.tbl"))

		# Create plink files
		markersIdsMap = gwaspToPlinkGenoMap (genotypeFile, "out/")
		plinkFile     = gwaspTetraGenoToPlinkPed (genotypeFile, markersIdsMap, "out/")

		# Recode to plink format adjusted for tassel and plink
		cmm = sprintf ("plink --file %s --recode tab --out %s", plinkFile, plinkFile)
		runCommand (cmm, "log-filtering.log")

		# Copy links of filtered plink files to main dir
		runCommand (sprintf ("ln -s %s.ped out/filtered-plink-genotype.ped", basename (plinkFile)), "log-filtering.log")
		runCommand (sprintf ("ln -s %s.map out/filtered-plink-genotype.map", basename (plinkFile)), "log-filtering.log")
	}
	else {
		msg (    "With filters")
		# Apply filters to genotype (markers and samples) by calling external program
		filtered = filterByQCFilters (common$genotypeFile, common$phenotypeFile, config)
		genotypeFile  = filtered$genotypeFile
		phenotypeFile = filtered$phenotypeFile
	}

	msg  (">>>> Evaluating trait ", common$trait)
	msg ("Converting and writing plink pheno filtered to other tools formats...")
	msg (phenotypeFile)
	gwasp2plinkPhenotype  (phenotypeFile,"out/filtered-plink-phenotype.tbl") 
	gwasp2tasselPhenotype (phenotypeFile,"out/filtered-tassel-phenotype.tbl") 
	#gwaspToShesisGenoPheno ("out/filtered-gwasp4-genotype.tbl", "out/filtered-gwasp4-phenotype.tbl") 
	gwaspToShesisGenoPheno (genotypeFile, phenotypeFile)
	return (list (genotypeFile=genotypeFile, phenotypeFile=phenotypeFile, trait=common$trait))
}
#-------------------------------------------------------------
# Filter by common sample names
#-------------------------------------------------------------
filterByCommonNames <- function (genotypeFile, phenotypeFile) 
{
	geno  = read.csv (file=genotypeFile, header=T)
	pheno = read.csv (file=phenohenotypeFile, header=T)
	rownames (pheno) = pheno [,1]

	genoColumns   = colnames (geno)
	phenoSamples  = pheno [,1] 
	commonSamples <- intersect (genoColumns[-(1:3)], phenoSamples) 

	genoCommon  <- geno  [,c(genoColumns[1:3], commonSamples)]
	map = genoCommon [, (1:3)]
	write.table (file="out/map.tbl", map)

	phenoCommon <- pheno [commonSamples,]
	trait  <- colnames (phenoCommon)[2]
	genoCommonFile  = paste0 ("out/", addLabel (genotypeFile, "COMMON"))
	phenoCommonFile = paste0 ("out/", addLabel (phenotypeFile, "COMMON"))
	write.csv (file=genoCommonFile, genoCommon, quote=F, row.names=F)
	write.csv (file=phenoCommonFile, phenoCommon, quote=F, row.names=F)

	return (list (genotypeFile=genoCommonFile, phenotypeFile=phenoCommonFile, trait=trait))
}
#-------------------------------------------------------------
# Filter by missing markers and samples, MAF, and HWE
# Apply filters to genotype (markers and samples) by calling external program
#-------------------------------------------------------------
filterByQCFilters <- function (genotypeFile, phenotypeFile, config) 
{
	# Format convertion from gwasp4 to plink2
	msg();msg ("Converting gwaspoly to plink formats...")
	markersIdsMap = gwaspToPlinkGenoMap (genotypeFile, "out/")
	plinkFile     = gwaspTetraGenoToPlinkPed (genotypeFile, markersIdsMap, "out/")

	cmm = paste ("plink --file", plinkFile, "--make-bed", "--out", paste0(plinkFile,"-QC"))
	if (config$filtering==T) {
		msg ("    >>>> Filtering by missingness, MAF, and HWE")
		# Filter missingness per sample (MIND)"
		if (!is.null(config$MIND)) cmm=paste (cmm, paste ("--mind", config$MIND))
		# Filter missingness per SNP    (GENO)
		if (!is.null(config$GENO)) cmm=paste (cmm, paste ("--geno", config$GENO))
		# Filter SNPs with a low minor allele frequency (MAF)
		if (!is.null(config$MAF)) cmm=paste (cmm, paste ("--maf", config$MAF))
		# Filter SNPs which are not in Hardy-Weinberg equilibrium (HWE).
		if (!is.null(config$HWE)) cmm=paste (cmm, paste ("--hwe", config$HWE))
	}

	# Recode to plink format adjusted for tassel and plink
	runCommand (cmm, "log-filtering.log" )
	cmm = sprintf ("plink --bfile %s-QC --recode tab --out %s-QC", plinkFile, plinkFile)
	runCommand (cmm, "log-filtering.log")

#	# Check for linkage disequilibrium
#	msg ("    >>>> Checking form linkage disequilibrium...")
#	cmm = sprintf ("plink --bfile %s-QC --indep-pairwise 50 10 0.9 --out %s-LD", plinkFile, plinkFile)
#	runCommand (cmm, "log-filtering.log")
#	cmm = sprintf ("plink --bfile %s-QC --exclude %s-LD.prune.out --make-bed --out %s-LDB", plinkFile, plinkFile, plinkFile)
#	runCommand (cmm, "log-filtering.log")
#	cmm = sprintf ("plink --bfile %s-LDB --recode tab --out %s-LD", plinkFile, plinkFile)
#	runCommand (cmm, "log-filtering.log")

	# Copy links of filtered plink files to main dir
	runCommand (sprintf ("ln -s %s/%s-QC.ped out/filtered-plink-genotype.ped", getwd(), plinkFile), "log-filtering.log")
	runCommand (sprintf ("ln -s %s/%s-QC.map out/filtered-plink-genotype.map", getwd(), plinkFile), "log-filtering.log")

	# Get final markers and individuals"
	msg ("Writting filtered markers and individuals...")
	filteredSamples <- as.character (read.table ("out/filtered-plink-genotype.ped")[,2])
	filteredMarkers <- as.character (read.table ("out/filtered-plink-genotype.map")[,2])

	# Filter phenotype
	#phenoAll = read.csv (phenotypeFile, header=T, check.names=F)
	phenoAll = read.csv (phenotypeFile, header=T)
	rownames (phenoAll) = phenoAll [,1]
	phenoFiltered  = phenoAll [filteredSamples,]
	trait  <- colnames (phenoAll)[2]

	# Filter genotype
	#genoAll  = read.csv (genotypeFile, header=T, check.names=F)
	genoAll  = read.csv (genotypeFile, header=T)
	rownames (genoAll) = genoAll [,1]
	genoFiltered <- genoAll [filteredMarkers,]

	msg("Writting geno/pheno filtered by MAF, Missing, HWE...")
	outGenoFile  <- "out/filtered-gwasp4-genotype.tbl"
	outPhenoFile <- "out/filtered-gwasp4-phenotype.tbl"

	write.table (file=outGenoFile, genoFiltered, row.names=F, quote=F, sep=",")
	write.table (file=outPhenoFile, phenoFiltered, row.names=F, quote=F, sep=",")

	return (list (genotypeFile=outGenoFile, phenotypeFile=outPhenoFile, trait=trait))
}

##-------------------------------------------------------------
## Calculate the inflation factor from -log10 values
##-------------------------------------------------------------
#calculateInflationFactor <- function (scores)
#{
#	remove <- which(is.na(scores))
#	if (length(remove)>0) 
#		x <- sort(scores[-remove],decreasing=TRUE)
#	else 
#		x <- sort(scores,decreasing=TRUE)
#
#	pvalues = 10^-x
#	chisq <- na.omit (qchisq(1-pvalues,1))
#	delta  = round (median(chisq)/qchisq(0.5,1), 3)
#
#	return (list(delta=delta, scores=x))
#}

#-------------------------------------------------------------
# Get params from config file and define models according to ploidy
#-------------------------------------------------------------
getConfigurationParameters <- function (configFile) 
{
	msg("Reading config file...")

	params = config::get (file=configFile) 

	message ("------------------------------------------------")
	message ("Summary of input parameters:")
	message ("------------------------------------------------")
	msg ("Genotype filename  : ", params$genotypeFile) 
	msg ("Phenotype filename : ", params$phenotypeFile) 
	msg ("Significance level : ", params$significanceLevel) 
	msg ("Correction method  : ", params$correctionMethod) 
	msg ("Trait              : ", params$trait) 
	msg ("GwAS model         : ", params$gwasModel) 
	msg ("Filtering          : ", params$filtering) 
	msg ("MIND               : ", params$MIND) 
	msg ("GENO               : ", params$GENO) 
	msg ("MAF                : ", params$MAF) 
	msg ("HWE                : ", params$HWE) 
	message ("------------------------------------------------")

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
#-------------------------------------------------------------
# Calculate threshold to decide SNPs significance
#-------------------------------------------------------------
calculateThreshold <- function (level, scores, method="FDR") 
{
	scores <- as.vector(na.omit (scores))
	m <- length(scores)
	if (method=="Bonferroni") 
		threshold <- -log10(level/m)
	else if (method=="FDR") {
		tmp <- cbind(10^(-scores),.qvalue(10^(-scores)))
		tmp <- tmp[order(tmp[,2]),]
		if (tmp[1,2] > level) {
			threshold <- -log10(tmp[1,1])*1.2
		} else {
			k <- max(which(tmp[,2] < level))
			threshold <- -log10(mean(tmp[k:(k+1),1]))
		}
	}else
		stop (paste0 ("Unknown correction method: ", method))

	return (threshold)
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
filterByMAFCommonNames <- function(geno.file, pheno.file, thresholdMAF=0.0){
	format    = "ACGT"	
	n.traits  = 1
	delim     = ","
	ploidy    = 4
	bases     = c("A","C","G","T")
	
	msg ("Reading genotype and phenotype....")
	geno      <- read.table(file=geno.file,header=T,as.is=T,check.names=T,sep=delim)
	gid.geno  <- colnames(geno)[-(1:3)]
	pheno     <- read.table(file=pheno.file,header=T,as.is=T,check.names=T,sep=delim)
	gid.pheno <- unique(pheno[,1])

	msg ("Removing duplicated markers")
	geno <- geno [!duplicated (geno[,1]),]  ### remove duplicates from geno

	map     <- data.frame(Marker=geno[,1],Chrom=factor(geno[,2],ordered=T),Position=geno[,3],stringsAsFactors=F)

	markers <- as.matrix(geno[,-(1:3)])
	rownames(markers) <- geno[,1]
	
	tmp <- apply(markers,1,get.ref)
	map$Ref <- tmp[1,]
	map$Alt <- tmp[2,]

	msg ("Calculating numeric genotype matrix...")
	M <- apply(cbind(map$Ref,markers),1,function(x){
		y <- gregexpr(pattern=x[1],text=x[-1],fixed=T)  
		ans <- as.integer(lapply(y,function(z){ifelse(z[1]<0,ploidy,ploidy-length(z))}))	
		return(ans)
	})

	msg ("Checking invalid marker calls...")
	rownames(M) <- gid.geno
	bad <- length(which(!is.element(na.omit(M),0:ploidy)))
	if (bad > 0) {stop("Invalid marker calls.")}
	
	msg ("Checking minor allele frecuency, MAF=", thresholdMAF)
	MAF <- apply(M,2,function(x){AF <- mean(x,na.rm=T)/ploidy;MAF <- ifelse(AF > 0.5,1-AF,AF)})
	polymorphic <- which(MAF>thresholdMAF)

	M <- M[,polymorphic]
	map <- map[polymorphic,]
	map <- map[order(map$Chrom,map$Position),]
	M <- M[,map$Marker]
	m <- nrow(map)
	msg("Number of polymorphic markers:",m,"\n")
	
	missing <- which(is.na(M))
	if (length(missing)>0) {
		msg("Missing marker data imputed with population mode...")
		M <- apply(M,2,impute.mode)
	}
	
	msg("Matching genotypic and phenotypic data...")
	gid <- intersect(gid.pheno, gid.geno)
	pheno <- pheno[is.element(pheno[,1],gid),]
	rownames (pheno) = pheno [,1]
	M <- M[gid,]
	N <- length(gid)
	msg ("N =",N,"individuals with phenotypic and genotypic information \n")
	
	n.fixed <- ncol(pheno) - n.traits - 1
	if (n.fixed > 0) {
		fixed <- data.frame(pheno[,(n.traits+2):ncol(pheno)],stringsAsFactors=F)
		fixed.names <- colnames(pheno)[(n.traits+2):ncol(pheno)]
		colnames(fixed) <- fixed.names
		pheno <- data.frame(pheno[,1:(1+n.traits)],stringsAsFactors=F)
		cat(paste("Detected following fixed effects:\n",paste(fixed.names,collapse="\n"),"\n",sep=""))
	} else {
		fixed <- data.frame(NULL)
	}
	trait <- colnames(pheno)[-1]
	msg("Evaluating following trait: ", trait) 
	

	msg ("Writing geno/pheno filtered by MAF, duplicated, common names")
	rownames (geno) = geno [,1]
	genoCommon      = geno [colnames(M),c(colnames(geno)[1:3], pheno[,1])]
	phenoCommon     = pheno 
	genoCommonFile  = paste0 ("out/", addLabel (geno.file, "COMMON"))
	phenoCommonFile = paste0 ("out/", addLabel (pheno.file, "COMMON"))
	write.csv (file=genoCommonFile, genoCommon, quote=F, row.names=F)
	write.csv (file=phenoCommonFile, phenoCommon, quote=F, row.names=F)
	# Write chromosome info 
	map = genoCommon [, (1:3)]
	write.table (file="out/map.tbl", map)
	# construct GWASpoly data structure
	gwaspolyData = new("GWASpoly",map=map,pheno=pheno,fixed=fixed,geno=M,ploidy=ploidy)
	return (list (genotypeFile=genoCommonFile, phenotypeFile=phenoCommonFile, data=gwaspolyData, trait=trait))
}

#-------------------------------------------------------------
# Add label to filename and new extension (optional)
#-------------------------------------------------------------
addLabel <- function (filename, label, newExt=NULL)  {
	nameext = strsplit (filename, split="[.]")
	name    = nameext [[1]][1] 
	if (is.null (newExt))
		ext     = nameext [[1]][2] 
	else
		ext     = newExt
	newName = paste0 (nameext [[1]][1], "-", label, ".", ext )
	return (newName)
}
#-------------------------------------------------------------
# Add label to filename
#-------------------------------------------------------------
old_addLabel <- function (filename, label)  {
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

#----------------------------------------------------------
# Util to print head of data
# Fast head, for debug only
#----------------------------------------------------------
hd <- function (data, m=10,n=10) {
	msg (deparse (substitute (data)),":")
	if (is.null (dim (data)))
		print (data [1:10])
	else if (ncol (data) < 10) 
		print (data[1:m,])
	else if (nrow (data) < 10)
		print (data[,1:n])
	else 
		print (data [1:m, 1:n])
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
# Run a command string using system function and writes output to log file
#-------------------------------------------------------------
runCommand <- function (command, logFile="gwas.log") 
{
	if (DEBUG==T) {
		msg (">>>> ", command)
		system (command)
	}else {
		msg ("<<<< ", command)
		#msg (">>>> ", command)
		errorsLog = paste0 (strsplit(logFile, split="[.]")[[1]], ".errors")
		system (paste0 (command, " > ", logFile," 2> ",errorsLog))
	}
}


#-------------------------------------------------------------
# Call main 
#-------------------------------------------------------------
main (args)


