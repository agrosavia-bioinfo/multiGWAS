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
	params = tryCatch (yaml.load_file (paramsFile), 
					   error=function (cond){
						   stop ("Some errors in configuration file. Check for valid or repeated names!!!")
					   })
	params = checkParameterR2 (params)         # default 1.0
	params = checkParameterGeneAction (params) # default "additive"
	params = checkParameterTraitType (params)  # default "quantitative"
	params = checkParameterCorrectionMethod (params)
	print (params);quit ()

	params$paramsFilename = paramsFile

	# Change to lower case text parameters
	params$genotypeFormat   = tolower (params$genotypeFormat) 
	params$gwasModel        = tolower (params$gwasModel) 
	params$filtering        = ifelse (tolower (params$filtering)=="true", T, F) 
	params$tools            = tolower (params$tools) 
	params$geneAction       = tolower (params$geneAction) 


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
# TraitType
#----------------------------------------------------------
checkParameterTraitType <- function (params) {
	if (is.null (params$traitType)) {
		message ("-----------------------------------------------------------------------------")
		message ('WARNING!! "traitType" parameter not defined, assumed traitType = quantitative')
		message ("-----------------------------------------------------------------------------")
		params$traitType = "quantitative"
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

	return (params)
}


#----------------------------------------------------------
# Call to main
#----------------------------------------------------------
#main ()
