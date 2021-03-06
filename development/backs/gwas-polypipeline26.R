#!/usr/bin/Rscript
# r2.41: Changed to parameters by config file. Beginning changes...
# r2.3: Removed annotations as it fails. Before config by file instead parameters
# r2.2: Write model type and inflation factor to significative QTLs
# r2.1: Support for multiple command line arguments (dynamic)

args = commandArgs(trailingOnly = TRUE)
#args = c("config-naive-4g-PATB96.config")

USAGE="USAGE: Rscript gwas-polypiline.R <config file>"
if (length (args) != 1) {
	message (USAGE)
	quit()
}
#-------------------------------------------------------------
#-------------------------------------------------------------

library (GWASpoly) #
suppressMessages(library (config))  # For read config file

setClass ("GWASpolyStruct", slots=c(params="list"),#testModels="character" snpModels="character", 
		  contains="GWASpoly.K")

options (width=300)
source ("gwas-formats.R")           # Module with functions to convert between different genotype formats 

#gwasModelsTypes        = c("Naive", "Kinship", "Structure", "PCs", "Kinship+Structure", "Kinship+PCs", "Structure+PCs", "Kinship+Structure+PCs")
#snpModels              = c("general","additive","1-dom", "2-dom")
multipleTestingModel    = "Bonferroni" # Bonferroni | FDR | permute
#snpModels              = c("general","additive","1-dom")
#testModels             = c("general", "additive","1-dom-alt","1-dom-ref")
#snpModels = testModels = c("general")
#-------------------------------------------------------------
# Global configs
#-------------------------------------------------------------
#genotypeFormati = "numeric"|"AB"|"ACGT"
#genotypeFormat = "numeric"
#gwasModel = "Naive"|"Structure"|"Kinship"|"Kinship+Structure|PCs"
#gwasModel      = "Kinship"
#PLOIDY         = 4
#-------------------------------------------------------------
data = data1 = data2 = data3 = data4 = NULL
phenotype = genotype = structure = NULL
dt=NULL
#-------------------------------------------------------------
# Main
#-------------------------------------------------------------
main <- function (args) 
{
	# Read and check config file arguments
	configFile =  args [1]
	params = getConfigurationParameters (configFile)
	genotypeFile  = params$genotypeFile
	phenotypeFile = params$phenotypeFile
	structureFile = params$structureFile

	msg (">>>>>>>>>>>>", params$gwasModel, "<<<<<<<<<<<") 

	# Load cache data
	if (file.exists ("gwas.RData")) load (file="gwas.RData")

	# Read, filter, and check phenotype and genotype
	data <- dataPreprocessing (genotypeFile, phenotypeFile, structureFile)

	# Read input genotype and genotype (format: "numeric", "AB", or "ACGT")
	data1 <- initGWAS (phenotypeFile, genotypeFile,
					   params$genotypePloidy, params$genotypeFormat, data1)
	# Set the kinship
	data2 <<- setKinship (data1, params$gwasModel, data2)

	# Populations structure and kinship
	data3 <- setPopulationStructure (data2, params$gwasModel, data$phenotype, 
									data$structure, data3)

	# GWAS execution
	data4 <- runGwaspoly (data3, params$gwasModel, params$snpModels, data4)

	save(data, data1, data2, data3, data4, file="gwas.RData") 

	# Plot results
	#if (params$genotypePloidy==4) ploidyLabel = "Tetra" else ploidyLabel = "Diplo"
	showResults (data4, params$testModels, data$trait, params$gwasModel, 
				 params$phenotypeFile, params$annotationsFile, params$genotypePloidy)

	runPlinkGwas (genotypeFile, phenotypeFile)
	runShesisGwas ()
}


#-------------------------------------------------------------
# 
#-------------------------------------------------------------
runPlinkGwas <- function (genotypeFile, phenotypeFile, model) 
{
	if (model=="Naive") {
		cmm = sprintf ("plink --file %s --linear --assoc --adjust --pheno %s --all-pheno 
					   --allow-no-sex --out out/out-plink-%s", genotypeFile, phenotypeFile, model)
		msg (cmm)
		system (cmm)
	}else if (model=="Kinship") {
		cmm = "plink --file filtered-plink-genotype --make-bed --out out/tmp-plinkb"
		cmm = "plink2 --bfile tmp-plinkb --make-pgen --out tmp-plinkg"
		cmm = "plink2 --pfile tmp-plinkg --king-cutoff 0.177 --out plinkg"
	}
}
	
runShesisGwas <- function () {
	inGenoPheno = "out/filtered-shesis-genopheno.tbl"
	inMarkers   = "out/filtered-shesis-markernames.tbl"
	outFile     = "out/out-shesis-naive"
	cmm=sprintf ("shesis --input %s --ploidy 4 --assoc --qtl --snpname-file %s --report-txt --adjust  --output %s", 
				 inGenoPheno, inMarkers, outFile)
	msg (cmm)
	system (cmm)
}

#-------------------------------------------------------------
# Filters the genotype by different quality control filters
# Read and check files, sample samples
#-------------------------------------------------------------
dataPreprocessing <- function (genotypeFile, phenotypeFile, structureFile) 
{
	# Format convertion from gwasp4 to plink2
	msg ("Converting gwaspoly to plink formats...")
	markersIdsMap = gwasp2plinkGenoMap (genotypeFile)
	plinkFile = gwaspTetraGenoToPlinkPed (genotypeFile, markersIdsMap)

	# Apply filters to genotype (markers and samples) by calling external program
	system2 ("gwas-filters.sh", plinkFile) 
	msg ("Reading filtered markers and individuals...")
	filteredMarkers <- read.table (file="out/filtered-names-markers.tbl", stringsAsFactors=F)[,1]
	filteredSamples <- read.table (file="out/filtered-names-samples.tbl", stringsAsFactors=F)[,1]

	# Filter phenotype
	msg ("Step 1:", phenotypeFile)
	phenoAll = read.csv (phenotypeFile, header=T, sep=",", check.names=F)
	phenoFiltered  = phenoAll [phenoAll[,1] %in% filteredSamples,]

	# Filter genotype
	msg ("Step 2:")
	genoAll  = read.csv (genotypeFile, header=T, sep=",", check.names=F)
	rownames (genoAll) = genoAll [,1]
	genoFiltered <- genoAll [filteredMarkers,]

	msg ("Step 3:")
	# Select common names geno and pheno
	samplesNamesGeno   <- colnames (genoFiltered)
	samplesNamesPheno  = phenoFiltered[,1]
	commonSamples <- Reduce (intersect, list (samplesNamesGeno, samplesNamesPheno)) 

	# Filter geno and pheno by common names
	phenoCommon = phenoFiltered [phenoFiltered[,1] %in% commonSamples,]
	genoColums  = c("Markers","Chrom","Position", commonSamples)
	genoCommon  = genoFiltered  [,colnames(genoFiltered) %in% genoColums]

	outGenoFile  = "out/filtered-gwasp4-genotype.tbl"
	outPhenoFile = "out/filtered-gwasp4-phenotype.tbl"
	write.table (file=outGenoFile, genoCommon, row.names=F, quote=F, sep=",")
	write.table (file=outPhenoFile, phenoCommon, row.names=F, quote=F, sep=",")

	trait  <- colnames (phenoCommon)[2]
	msg  (">>>> Evaluating trait ", trait)

	# Convert and write plink pheno filtered
	gwasp2plinkPhenotype  (outPhenoFile,"out/filtered-plink-phenotype.tbl") 
	gwasp2tasselPhenotype (outPhenoFile,"out/filtered-tassel-phenotype.tbl") 
	gwasp2shesisGenoPheno ("out/filtered-gwasp4-genotype.tbl", "out/filtered-gwasp4-phenotype.tbl") 

	return (list (genotypeFile=outGenoFile, phenotypeFile=outPhenoFile,
				  structureFile=structureFile, trait=trait))

}
#-------------------------------------------------------------
# Read and check command line arguments
# Currently: phenotye, genotype, structure
#-------------------------------------------------------------
old_setGwaspolyModels <- function (params) 
{
	msg("Reading Files...")

	if (params$genotypePloidy=="4") {
		#testModels  = snpModels = c ("general")
		snpModels  = c("general","additive","1-dom", "2-dom")
		testModels = c("general", "additive","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")

		params = append (params, list (snpModels  = snpModels))
		params = append (params, list (testModels = testModels))
	}else{
		params = append (params, list (snpModels  = c("general","additive","1-dom")))
		params = append (params, list (testModels = c("general", "additive","1-dom-alt","1-dom-ref")))
	}

	print (params)
	message ("------------------------------------------------")
	message ("Summary of input parameters:")
	message ("------------------------------------------------")
	msg ("Phenotype filename : ", params$phenotypeFile) 
	msg ("Genotype filename  : ", params$genotypeFile) 
	msg ("Structure filename : ", params$structureFile) 
	msg ("SNPs filename      : ", params$annotationsFile) 
	msg ("GwAS model         : ", params$gwasModel) 
	msg ("Genotype format    : ", params$genotypeFormat) 
	msg ("Genotype ploidy    : ", params$genotypePloidy)
	msg ("Kinship            : ", params$kinship) 
	message ("------------------------------------------------")

	return (params)

}

#-------------------------------------------------------------
# Set the kinship matrix using the realized relationship matrix  
#-------------------------------------------------------------
setKinship <- function (data1,  gwasModel, data2) 
{
	msg();msg("Setting kinship...");msg()
 	# Load data instead calculate it
	if (!is.null (data2)) {msg ("Loading kinship..."); return (data2) }

	#kinshipMatrix = NULL
	if (gwasModel %in% c("Kinship", "Kinship+Structure", "Kinship+PCs")) {
		msg (">>>> With default kinship... ")
		kinshipMatrix = NULL
	}else { 
		msg (">>>> Without kinship...") 
		markerNames   = data1@pheno [,1]
		n             = length (markerNames)
		kinshipMatrix = matrix (diag (n), n, n, dimnames=list (markerNames, markerNames))
	}		
	data2  = set.K (data1, K=kinshipMatrix)
	return (data2)
}
		
#-------------------------------------------------------------
# Fix Populations structure and kinship
#-------------------------------------------------------------
setPopulationStructure <- function (data2, gwasModel, phenotype, structure, data3) 
{
	st <<- structure
	msg();msg ("Setting population structure...");msg()
	#setClass ("GWASpolyStruct", slots=c(params="list"), contains="GWASpoly.K")
	data3 <- new ("GWASpolyStruct", data2)
	
	nPCs=0
	structNames = structTypes = NULL

	if (gwasModel %in% c("PCs","Kinship+PCs")) {
		nPCs=5
		msg (">>>> With nPCS=", nPCs, "...")
	}else if (gwasModel %in% c("Structure", "Kinship+Structure")) {
		msg (">>>> With population structure...")
		structNames <- colnames (structure[-1])
		structTypes <- rep ("numeric", length (structNames))

		#data3@pheno = phenoStruct
		data3@pheno = phenotype
		data3@fixed = st[,-1]
	}else 
		msg (">>>> Without populations structure")

	data3@params = set.params(n.PC=nPCs, fixed=structNames, fixed.type=structTypes)

	return (data3)
}

#-------------------------------------------------------------
# GWAS execution
#-------------------------------------------------------------
runGwaspoly <- function (data3, gwasModel, snpModels, data4) 
{
	msg();msg("Running GWASpoly...");msg()

	if (!is.null (data4)) { msg (">>>> Loading GWASpoly..."); return (data4) }

 	if (gwasModel %in% c("Naive","Kinship")) {
		msg (">>>> Without params")
		data4 = GWASpoly(data3, models=snpModels, traits=NULL, params=NULL, n.core=4)
	}else {
		msg (">>>> With params")
		data4 = GWASpoly(data3, models=snpModels, traits=NULL, params=data3@params)
	}
	
	return (data4)
}

#-------------------------------------------------------------
# Plot results
#-------------------------------------------------------------
showResults <- function (data4, testModels, trait, gwasModel, phenotypeFile, snpsAnnFile, ploidy) 
{
	msg();msg ("Plotting results...");msg()
	phenoName = strsplit (phenotypeFile, split=".tbl")[[1]][1]
	plotName = sprintf("out/out-gwasp%s-%s-plots.pdf", ploidy, gwasModel)

	pdf (file=plotName, width=11, height=7)
	n = length (testModels)
	
	# QQ-plot Output
	op <- par(mfrow = c(2,n), oma=c(0,0,3,0))
	for (i in 1:length(testModels)) {
		#par (cex.main=0.5, cex.lab=0.5, cex.axis=0.5, ann=T)
		qqPlot(data4,trait=trait, model=testModels[i], cex=0.3)
	}

	# QTL Detection
	data5 = set.threshold (data4, method=multipleTestingModel,level=0.05,n.permute=100,n.core=3)
	#significativeQTLs = get.QTL (data5)
	significativeQTLs = getQTL (data5, snpsAnnFile, gwasModel, ploidy)

	msg (">>>>", "Writing results...")
	outFile = sprintf ("out/out-gwasp%s-%s-QTLs.tbl", ploidy, gwasModel) 
	write.table (file=outFile, significativeQTLs, quote=F, sep="\t", row.names=F)

	# Manhattan plot Output
	for (i in 1:length(testModels)) {
		#par (cex=1.5)
		manhattan.plot (y.max=20,data5, trait=trait, model=testModels [i])
	}
	plotTitle = sprintf ("GWAS %s-ploidy with %s for %s trait", ploidy, gwasModel, trait)  
	mtext(plotTitle, outer=T,  cex=1.5,  line=0)
	par(op)
	dev.off()
}

#-------------------------------------------------------------
# Extracts significant QTL
#-------------------------------------------------------------
getQTL <- function(data,snpsAnnFile, gwasModel, ploidy, traits=NULL,models=NULL) 
{
	stopifnot(inherits(data,"GWASpoly.thresh"))
	sc <<-data@scores

	if (is.null(traits)) traits <- names(data@scores)
	else stopifnot(is.element(traits,names(data@scores)))

	if (is.null(models)) models <- colnames(data@scores[[1]])
	else stopifnot(is.element(models,colnames(data@scores[[1]])))

	msg ("Reading associations...")
	if (!is.null (snpsAnnFile)) snpsAnnotations <- read.csv (file=snpsAnnFile, header=T)

	n.model <- length(models)
	n.trait <- length(traits)
	output <- data.frame(NULL)
	for (j in 1:n.model) {
		ix <- which(data@scores[[traits[1]]][,models[j]] > data@threshold[traits[1],models[j]])
		markers <-  data.frame (SNP=data@map[ix,c("Marker")])
		print (markers)
		if (!is.null (snpsAnnFile)) 
			snpAnn  <- merge (markers, snpsAnnotations, by.x="SNP",by.y="SNP_id", sort=F)[,c(2,7)]
		else
			snpAnn = "None"

		scores <- data@scores[[1]][,models[j]]
		datax = calculateInflationFactor (scores)

		n.ix <- length(ix)
		
		df = data.frame(Ploidy=rep (ploidy, n.ix),
						Type=rep (gwasModel, n.ix),
						GC=rep(datax$delta,n.ix), 
					    Model=rep(models[j],n.ix),
						Score=round(data@scores[[traits[1]]][ix,models[j]],2),
					    Threshold=round(rep(data@threshold[traits[1],models[j]],n.ix),2),
						Effect=round(data@effects[[traits[1]]][ix,models[j]],2),
						data@map[ix,])
						#snpAnn) 
						#stringsAsFactors=F,check.names=F)

		output <- rbind(output, df)
	}
	#out <-cbind (Type=gwasModel, output)
	output = output [order(-output$GC,-output$Score),]
	return(output)
}
#-------------------------------------------------------------
# QQ plot
#-------------------------------------------------------------
qqPlot <- function(data,trait,model,cex=1,filename=NULL) 
{
	stopifnot(inherits(data,"GWASpoly.fitted"))
	traits <- names(data@scores)
	stopifnot(is.element(trait,traits))
	models <- colnames(data@scores[[trait]])
	stopifnot(is.element(model,models))
	scores <- data@scores[[trait]][,model]

	datax = calculateInflationFactor (scores)

	n <- length(datax$scores)
	unif.p <- -log10(ppoints(n))
	if (!is.null(filename)) {postscript(file=filename,horizontal=FALSE)}
	par(pty="s")
	plot(unif.p, datax$scores, pch=16,cex=cex,
		 xlab=expression(paste("Expected -log"[10],"(p)",sep="")),
		 ylab=expression(paste("Observed -log"[10],"(p)",sep="")),
		 main=paste(trait," (",model,") ",sep=""))

	mtext (bquote(lambda[GC] == .(datax$delta)), side=3, line=-2, cex=0.7)

	lines(c(0,max(unif.p)),c(0,max(unif.p)),lty=2)
	if (!is.null(filename)) {dev.off()}
	return(datax$delta)
}

#-------------------------------------------------------------
# Calculate the inflation factor from -log10 values
#-------------------------------------------------------------
calculateInflationFactor <- function (scores)
{
	remove <- which(is.na(scores))
	if (length(remove)>0) 
		x <- sort(scores[-remove],decreasing=TRUE)
	else 
		x <- sort(scores,decreasing=TRUE)

	pvalues = 10^-x
	chisq <- na.omit (qchisq(1-pvalues,1))
	delta  = round (median(chisq)/qchisq(0.5,1), 3)

	return (list(delta=delta, scores=x))
}

#-------------------------------------------------------------
# Read input genotype and genotype (format: "numeric", "AB", or "ACGT")
#-------------------------------------------------------------
initGWAS <- function (phenotypeFile, genotypeFile, ploidy, format, data1) 
{
	msg();msg ("Initializing GWAS...");msg()
	# When data is previously loaded
	if (!is.null (data)) {msg(">>>> Loading GWAS data..."); return (data)}

	data1 <- read.GWASpoly (ploidy = ploidy, pheno.file = phenotypeFile, 
							geno.file = genotypeFile, format = format, n.traits = 1, delim=",")

	return (data1)
}

#-------------------------------------------------------------
#-------------------------------------------------------------
getConfigurationParameters <- function (configFile) 
{
	msg("Reading config file...")

	params = config::get (file=configFile)

	if (params$genotypePloidy=="4") {
		#testModels  = snpModels = c ("general")
		snpModels = c("general","additive","1-dom", "2-dom")
		testModels = c("general", "additive","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")
	}else{
		snpModels  = c("general","additive","1-dom")
		testModels = c("general", "additive","1-dom-alt","1-dom-ref")
	}
	params = append (params, list (snpModels=snpModels))
	params = append (params, list (testModels=testModels))

	print (params)
	message ("------------------------------------------------")
	message ("Summary of input parameters:")
	message ("------------------------------------------------")
	msg ("Phenotype filename : ", params$phenotypeFile) 
	msg ("Genotype filename  : ", params$genotypeFile) 
	msg ("Structure filename : ", params$structureFile) 
	msg ("SNPs filename      : ", params$annotationsFile) 
	msg ("GwAS model         : ", params$gwasModel) 
	msg ("Genotype format    : ", params$genotypeFormat) 
	msg ("Genotype ploidy    : ", params$genotypePloidy)
	msg ("Kinship            : ", params$kinship) 
	message ("------------------------------------------------")

	return (params)
}
#-------------------------------------------------------------
# Print a log message with the parameter
#-------------------------------------------------------------
msg <- function (...) 
{
  messages = unlist (list (...))
  cat (">>>>", messages, "\n")
}

#-------------------------------------------------------------
# Call main 
#-------------------------------------------------------------
main (args)


