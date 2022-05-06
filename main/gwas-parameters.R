#!/usr/bin/Rscript

# Module for checking input parameters before running multiGWAS
 

`%notin%` <- Negate(`%in%`)
#-------------------------------------------------------------
#-------------------------------------------------------------
main <- function () {
	library (yaml)
	params = readCheckConfigParameters ("multiGWAS.config")
	print (params)

}

#-------------------------------------------------------------
# Get params from config file and define models according to ploidy
#-------------------------------------------------------------
readCheckConfigParameters <- function (paramsFile) {
	msg ("Reading and checking config parameters...")
	
	# Check input files: config, geno, pheno, map
	params = checkInputFiles (paramsFile)

	checkPloidy (params)
	params = checkParameterR2 (params)         # default 1.0
	params = checkParameterGeneAction (params) # default "additive"
	params = checkParameterCorrectionMethod (params)
	params = checkParameterNonModelOrganism (params) 
	if (is.null (params$traitType)) params$traitType = "quantitative"

	# Change to lower case text parameters
	params$genotypeFormat   = tolower (params$genotypeFormat) 
	params$gwasModel        = tolower (params$gwasModel) 
	params$filtering        = ifelse (tolower (params$filtering)=="true", T, F) 
	params$nonModelOrganism = ifelse (tolower (params$nonModelOrganism)=="true", T, F) 
	params$tools            = tolower (params$tools) 
	params$geneAction       = tolower (params$geneAction) 

	# Create output dir, check input files, and copy files to output dir
	runningDir   = paste0 ("out-", strsplit (paramsFile, split="[.]") [[1]][1])
	createDir (runningDir)

	copyFilesToRunningDir (params, runningDir)

	# Change to the output dir and set global runningDir 
	setwd (runningDir)
	params$runningPath = getwd()
	return (params)
}

#----------------------------------------------------------
# Check non-model organisms and number of chromosomes 
#----------------------------------------------------------
checkParameterNonModelOrganism <- function (params) {
	msg ("Cheking parameters for non-model organism...")
	nonModel     = params$nonModelOrganism
	nChromosomes = params$numberOfChromosomes

	if (is.null (nonModel) || nonModel==FALSE) {# | nonModel==FALSE)
		msgmsg ("Non-model organism set to FALSE...")
		params$nonModelOrganism = FALSE
		return (params)
	}

	if (nonModel %notin% c(TRUE,FALSE)) {
		message ("------------------------------------------------------------------")
		message ("Error in nonModelOrganism parameter. Value must be TRUE or FALSE.")
		message ("------------------------------------------------------------------")
		quit()
	}

	if (is.null (nChromosomes) || !is.integer (as.integer (nChromosomes))) {
		message ("------------------------------------------------------------------")
		message ("Error in numberOfChromosomes parameter. It needs a valid number." )
		message ("------------------------------------------------------------------")
		quit()
	}

	return (params)
}

#----------------------------------------------------------
# Check possible errors in ploidy 
#----------------------------------------------------------
checkPloidy <- function (params) {
	if (params$ploidy %notin% c("2", "4"))   {
		message ("----------------------------------------------------------------------------")
		message ("Error: Ploidy: ", params$ploidy, " not supported")
		message ("----------------------------------------------------------------------------")
		quit ()
	}
}

#----------------------------------------------------------
# Check input files: params, genotype, phenotype, mapfile
#----------------------------------------------------------
checkInputFiles <- function (paramsFile) {
	msg ("Checking input files...")

	errMsg <- function (text) {
		message ("----------------------------------------------------------------------------")
		message (text)
		message ("----------------------------------------------------------------------------")
		quit() 
	}

	# Check configuration file
	if (file.exists (paramsFile)==F) 
		errMsg ("\nError: Configuration file not found\n")

	params = tryCatch (yaml.load_file (paramsFile), error=function (cond){
					errMsg (sprintf ("Error: Configuration file with invalid or repeated names!!!"))})

	# Check genotype file
	if (!file.exists (params$genotypeFile)) 
		errMsg (sprintf ("Error: Genotype file not found: '%s'", params$genotypeFile))

	# Check phenotype file
	if (!file.exists (params$phenotypeFile)) 
		errMsg (sprintf ("Error: Phenotype file not found: '%s'", params$phenotypeFile))

	# Check phenotype format
	if (ncol(read.csv(params$phenotypeFile))==1)
		errMsg ("Error: Phenotype is not in CSV format or it only has one column.")

	# Check genotype format
	format  = tolower (params$genotypeFormat)
	if (format %in% c("matrix","gwaspoly","updog"))
		if (ncol (read.csv(params$genotypeFile))==1)
			errMsg ("Error: Genotype is not in CSV format or it only has one column.")

	# Check maps file
	if (tolower (params$genotypeFormat) %in% c("matrix", "fitpoly", "updog")) 
		if (is.null (params$mapFile) || !file.exists (params$mapFile))
			errMsg ("Error: Map file not found or not specified in the config file")

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

#----------------------------------------------------------
# R2
#----------------------------------------------------------
checkParameterR2 <- function (params) {
	if (is.null (params$R2)) {
		message ("----------------------------------------------------")
		message ('WARNING!! "R2" parameter not defined, assumed R2 = 1.0')
		message ("----------------------------------------------------")
		params$R2 <- 1.0
	}
	return (params)
}

#----------------------------------------------------------
# GeneAction
#----------------------------------------------------------
checkParameterGeneAction <- function (params) {
	if (is.null (params$geneAction)) {
		message ("---------------------------------------------------------------------------")
		message ('WARNING!! "geneAction" parameter not defined, assumed geneAction = additive')
		message ("---------------------------------------------------------------------------")
		params$geneAction = "additive"
	}
	return (params)
}

#----------------------------------------------------------
# Correction method: "FDR" or "Bonferroni"
#----------------------------------------------------------
checkParameterCorrectionMethod <- function (params) {
	if (is.null (params$correctionMethod)) {
		message ("------------------------------------------------------------------------------------------")
		message ('WARNING!! "correctionMethod" parameter not defined, assumed correctionMethod = Bonferroni ')
		message ("-------------------------------------------------------------------------------------------")
		params$correctionMethod = "bonferroni"
	}

	params$correctionMethod   = tolower (params$correctionMethod) 
	if (params$correctionMethod %notin% c("bonferroni", "fdr"))
		stop (paste0 ("MG Error: Unknown correction method: ", params$correctionMethod), call.=T)

	if (params$correctionMethod=="bonferroni")
		params$correctionMethod = "Bonferroni"
	if (params$correctionMethod=="fdr")
		params$correctionMethod = "FDR"

	return (params)
}


#----------------------------------------------------------
# Call to main
#----------------------------------------------------------
#main ()
