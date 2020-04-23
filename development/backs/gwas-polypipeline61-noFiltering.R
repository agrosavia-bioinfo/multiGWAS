#!/usr/bin/Rscript
# r5.5: Commands as bashs scripts, log files for stdout and errors
# r5.4: Tools in parallel and commands messages to log files
# r5.3: Controlling populaton structure in one function
# r5.0: Working with tassel pipeline for Naive (GLM) and K+Q (MLM)
# r4.8: With tassel pipeline for GLM 
# r4.1: Without structure file input, Working with G4,G2,PL,SH
# r2.41: Changed to parameters by config file. Beginning changes...

LOAD_DATA = FALSE
args = commandArgs(trailingOnly = TRUE)
#args = c("in/config-Gota-Naive-filtersNone-impute.config")
#args = c("in/config-test500-Naive-filtersNone.config")

USAGE="USAGE: Rscript gwas-polypiline.R <config file>"
if (length (args) != 1) {
	message (USAGE)
	quit()
}
#-------------------------------------------------------------
#-------------------------------------------------------------

library (GWASpoly) #
library (parallel) #
suppressMessages(library (config))  # For read config file

# New class for gwaspoly
setClass ("GWASpolyStruct", slots=c(params="list"), contains="GWASpoly.K")

options (width=300)
HOME = Sys.getenv ("GWASP_HOME")
source (paste0 (HOME, "/gwas-formats.R"))           # Module with functions to convert between different genotype formats 
source (paste0 (HOME, "/gwas-summary.R"))           # Module with functions to create summaries: tables and venn diagrams
source (paste0 (HOME, "/scripts/gwaspoly.R"))       # Module with gwaspoly functions

#-------------------------------------------------------------
# Global configs
#-------------------------------------------------------------
#genotypeFormats = "numeric"|"AB"|"ACGT"
#gwasModel = "Naive"|"Kinship+PCs"
#-------------------------------------------------------------
data = data1 = data2 = data3 = NULL
dt=NULL
#-------------------------------------------------------------
# Main
#-------------------------------------------------------------
main <- function (args) 
{
	# Read and check config file arguments
	configFile =  args [1]
	config <- getConfigurationParameters (configFile)
	genotypeFile  <- config$genotypeFile
	phenotypeFile <- config$phenotypeFile

	system ("rm -rf out;mkdir out")
	system (sprintf ("cp %s out/", configFile))

	msg (">>>>>>>>>>>>", config$gwasModel, "<<<<<<<<<<<") 

	# Read, filter, and check phenotype and genotype
	data <- dataPreprocessing (genotypeFile, phenotypeFile, config)
	config$genotypeFile  =data$genotypeFile
	config$phenotypeFile =data$phenotypeFile

	config$trait = data$trait

	# Run the four tools in parallel

	runGwaspGwas (config)
	#mclapply (c("Gwasp", "Plink", "Shesis", "Tassel"), runGwasTool, config, mc.cores=4)

	# Create outputs: tables, figures
	markersSummaryTable ("out/", config$gwasModel, "out/", nBEST=5)

	# Save test in new dir
	#outName = gsub ("config","test", configFile)
	#system (sprintf ("mv out %s", outName))
}

#-------------------------------------------------------------
# Used to run in parallel the other functions
#-------------------------------------------------------------
runGwasTool <- function (tool, config) 
{
	if (tool=="Gwasp")
		runGwaspGwas (config)
	else if (tool=="Plink")
		runPlinkGwas (config)
	else if (tool=="Shesis")
		runShesisGwas (config)
	else if (tool=="Tassel")
		runTasselGwas (config)
	else
		stop ("Tool not supported")
}
#-------------------------------------------------------------
#-------------------------------------------------------------
runGwaspGwas <- function (params) 
{
	genotypeFile  = params$genotypeFile
	phenotypeFile = params$phenotypeFile

	if (LOAD_DATA & file.exists ("gwas.RData")) load (file="gwas.RData")

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
				 params$phenotypeFile, params$annotationsFile, ploidy)

	if (LOAD_DATA) save(data, data1, data2, data3, file="gwas.RData") 
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
		cmm=sprintf ("scripts/plink-Naive.sh %s %s %s", inGeno, inPheno, outFile)
	else if (model=="Structure") 
		cmm=sprintf ("scripts/plink-Structure-PCs.sh %s %s %s", inGeno, inPheno, outFile)
	else
		quit (paste ("Type of GWAS:\"", model, "\", not supported"))

	runCommand (cmm, "log-Plink.log")
	runCommand (paste0 ("mv ", outPlink, " ", outFile, ".scores"), "log-Plink.log") 
}
	
#-------------------------------------------------------------
# Shesis tool
#-------------------------------------------------------------
runShesisGwas <- function (params) 
{
	msg ("Running Shesis GWAS...")
	model = params$gwasModel

	inGenoPheno = "out/filtered-shesis-genopheno.tbl"
	inMarkers   = "out/filtered-shesis-markernames.tbl"
	outFile     = paste0 ("out/out-Shesis-", model)
	outShesis   = paste0(outFile,".txt")

	cmm=sprintf ("scripts/shesis-Naive.sh %s %s %s", inGenoPheno, inMarkers, outFile)
	msg (cmm)
	runCommand (cmm, "log-Shesis.log")
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
		cmm=sprintf ("scripts/tassel-pipeline-GLM-Naive.sh %s %s %s %s",
					 inGenoPED, inGenoMAP, inPhenoTBL, outPrefix)
		runCommand (cmm, "log-tassel.log")
		outFile   = list.files("out/", pattern=sprintf("^(.*(%s).*(1).*(txt)[^$]*)$",model), full.names=T)
	}else if (model=="Structure") {
		cmm=sprintf ("scripts/tassel-pipeline-MLM-Kinship_PCs.sh %s %s %s %s",
					 inGenoPED, inGenoMAP, inPhenoTBL, outPrefix)
		runCommand (cmm, "log-tassel.log")
		outFile   = list.files("out/", pattern=sprintf("^(.*(%s).*(stats).*(txt)[^$]*)$",model), full.names=T)
	}else 
		quit (paste ("Type of GWAS:\"", model, "\", not implemented"))
	
	# Rename output file
	msg ("Tassel output file: ", outFile)
	outTassel = sprintf ("out/out-Tassel-%s.scores", model)

	# Data preprocessing: Test correction, sort and rename output file
	ts      = read.table (file=outFile, header=T, sep="\t")
	tsAdj   = ts %>% mutate (adjP=p.adjust(p,"fdr"),adjPadd=p.adjust(add_p, "fdr"),adjPdom=p.adjust(dom_p, "fdr"))
	tsMin   = tsAdj %>% rowwise %>% mutate (minP=min(adjP, adjPadd, adjPdom, na.rm=T)) 
	tsArr   = tsMin %>% arrange (minP)
	write.table (file=outTassel, tsArr, quote=F, sep="\t", row.names=F)
}

#-------------------------------------------------------------
# Filters the genotype by different quality control filters
# Read and check files, sample samples
# Convert geno an pheno to other tool formats 
#-------------------------------------------------------------
dataPreprocessing <- function (genotypeFile, phenotypeFile, config) 
{
	if (config$filtering == F) 
		return (list (genotypeFile=genotypeFile, phenotypeFile=phenotypeFile, trait="tuber_shape"))
	else {
		# First filtering: imputation, MAF, common names

		imputedFiles <-filterByImputationMAFNames (phenotypeFile, genotypeFile, config)
		#genotypeFile  = imputedFiles$genof
		#phenotypeFile = imputedFiles$phenof

		# Format convertion from gwasp4 to plink2
		msg();msg ("Converting gwaspoly to plink formats...")
		markersIdsMap = gwaspToPlinkGenoMap (genotypeFile, "out/")
		plinkFile = gwaspTetraGenoToPlinkPed (genotypeFile, markersIdsMap, "out/")

		# Apply filters to genotype (markers and samples) by calling external program
		filtered = filterByMissingHWE (plinkFile, config)

		# Filter phenotype
		phenoAll = read.csv (phenotypeFile, header=T, sep=",", check.names=F)
		write.csv (file="tmp-phenoNames.tbl", phenoAll [,1],quote=F,row.names=F)
		write.csv (file="tmp-sampleNames.tbl", filtered$Samples,quote=F,row.names=F)
		rownames (phenoAll) = phenoAll [,1]
		print (dim (phenoAll))
		print (length (filtered$Samples))
		phenoFiltered  = phenoAll [intersect (filtered$Samples, phenoAll[,1]),]
		print (dim (phenoFiltered))
		print (phenoFiltered)
		quit()

		# Filter genotype
		genoAll  = read.csv (genotypeFile, header=T, sep=",", check.names=F)
		rownames (genoAll) = genoAll [,1]
		genoFiltered <- genoAll [filtered$Markers,]

		# Select common names geno and pheno
		samplesNamesGeno   <- colnames (genoFiltered)
		samplesNamesPheno  <- phenoFiltered[,1]
		commonSamples <- Reduce (intersect, list (samplesNamesGeno, samplesNamesPheno)) 

		# Filter geno and pheno by common names
		phenoCommon <- phenoFiltered [phenoFiltered[,1] %in% commonSamples,]
		genoColums  <- c(samplesNamesGeno[1:3], commonSamples)
		genoCommon  <- genoFiltered  [,colnames(genoFiltered) %in% genoColums]

		outGenoFile  <- "out/filtered-gwasp4-genotype.tbl"
		outPhenoFile <- "out/filtered-gwasp4-phenotype.tbl"

		write.table (file=outGenoFile, genoCommon, row.names=F, quote=F, sep=",")
		write.table (file=outPhenoFile, phenoCommon, row.names=F, quote=F, sep=",")

		trait  <- colnames (phenoCommon)[2]
		msg  (">>>> Evaluating trait ", trait)

		# Convert and write plink pheno filtered
		msg (outPhenoFile)
		gwasp2plinkPhenotype  (outPhenoFile,"out/filtered-plink-phenotype.tbl") 
		gwasp2tasselPhenotype (outPhenoFile,"out/filtered-tassel-phenotype.tbl") 
		gwaspToShesisGenoPheno ("out/filtered-gwasp4-genotype.tbl", "out/filtered-gwasp4-phenotype.tbl") 

		return (list (genotypeFile=outGenoFile, phenotypeFile=outPhenoFile, trait=trait))
	}
}
#-------------------------------------------------------------
# Apply filters to genotype (markers and samples) by calling external program
#-------------------------------------------------------------
filterByMissingHWE <- function (plinkFile, config) 
{
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

	runCommand (cmm, "log-filtering.log" )
	cmm = sprintf ("plink --bfile %s-QC --recode tab --out %s-QC", plinkFile, plinkFile)
	runCommand (cmm, "log-filtering.log")

	# Copy links of filtered plink files to main dir
	runCommand (sprintf ("ln -s %s/%s-QC.ped out/filtered-plink-genotype.ped", getwd(), plinkFile), "log-filtering.log")
	runCommand (sprintf ("ln -s %s/%s-QC.map out/filtered-plink-genotype.map", getwd(), plinkFile), "log-filtering.log")

	# Get final markers and individuals"
	#sampleNames = read.table (paste0(plinkFile,"-QC.ped"))[,2]
	sampleNames = read.table ("out/filtered-plink-genotype.ped")[,2]
	write.table (file="out/filtered-names-samples.tbl", sampleNames, quote=F,row.names=F,col.names=F)
	#markerNames = read.table (paste0(plinkFile,"-QC.map"))[,2]
	markerNames = read.table ("out/filtered-plink-genotype.map")[,2]
	write.table (file="out/filtered-names-markers.tbl", markerNames, quote=F,row.names=F,col.names=F)

	msg ("Reading filtered markers and individuals...")
	filteredMarkers <- read.table (file="out/filtered-names-markers.tbl", stringsAsFactors=F)[,1]
	filteredSamples <- read.table (file="out/filtered-names-samples.tbl", stringsAsFactors=F)[,1]

	return (list (Markers=filteredMarkers, Samples=filteredSamples))
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
	msg ("GwAS model         : ", params$gwasModel) 
	msg ("Imputation         : ", params$imputation) 
	msg ("Filtering          : ", params$filtering) 
	msg ("MIND               : ", params$MIND) 
	msg ("GENO               : ", params$GENO) 
	msg ("MAF                : ", params$MAF) 
	msg ("HWE                : ", params$HWE) 
	message ("------------------------------------------------")

	runCommand (sprintf ("cp %s out/", configFile))

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
# Impute, filter by MAF, unify geno and pheno names
# Only for "ACGT" format (For other formats see Gwaspoly sources)
#-------------------------------------------------------------
filterByImputationMAFNames <- function(pheno.file, geno.file, config){
	msg ("Filtering by gwaspoly...")

	bases <- c("A","C","G","T")
	n.traits = 1
	ploidy   = 4
	delim    = ","

	msg ("Reading genotype...")
	geno = read.csv(file=geno.file,header=T,as.is=T,check.names=F,sep=delim)
	
	geno = geno [!duplicated (geno[,1]),]  ### remove duplicates from geno
	msg ("M = ", nrow (geno), " Initial markers")

	map     <- data.frame(Markers=geno[,1],Chrom=factor(geno[,2],ordered=T),Position=geno[,3],stringsAsFactors=F)
	markers <- as.matrix(geno[,-(1:3)])
	rownames(markers) <- geno[,1]
	
	msg ("Obtaining get/ref alleles...")
	tmp     <- apply(markers,1,get.ref)
	map$Ref <- tmp[1,]
	map$Alt <- tmp[2,]
	M <- apply(cbind(map$Ref,markers),1,function(x){
			y <- gregexpr(pattern=x[1],text=x[-1],fixed=T)  
			ans <- as.integer(lapply(y,function(z){ifelse(z[1]<0,ploidy,ploidy-length(z))}))	
			return(ans)
		})

	gid.geno <- colnames(geno)[-(1:3)]
	rownames(M) <- gid.geno

	bad <- length(which(!is.element(na.omit(M),0:ploidy)))
	if (bad > 0) {stop("Invalid marker calls.")}
	
	MAFTHRESHOLD=0.0
	if (config$filtering == T) MAFTHRESHOLD = config$MAF
	MAF <- apply(M,2,function(x){AF <- mean(x,na.rm=T)/ploidy;MAF <- ifelse(AF > 0.5,1-AF,AF)})
	#polymorphic <- which(MAF>MAFTHRESHOLD) ## LG
	polymorphic <- which(MAF>0.0)

	M <- M[,polymorphic]
	map <- map[polymorphic,]
	map <- map[order(map$Chrom,map$Position),]
	M <- M[,map$Marker]
	m <- nrow(map)
	cat(paste("Number of polymorphic markers:",m,"\n"))
	
	msg ("Imputing markers...")
	missing <- which(is.na(M))
	if (length(missing)>0) {
		cat("Missing marker data imputed with population mode \n")
		M <- apply(M,2,impute.mode)
	}
	
	### matching genotypic and phenotypic data 
	pheno     <- read.table(file=pheno.file,header=T,as.is=T,check.names=F,sep=delim)
	gid.pheno <- unique(pheno[,1])
	gid       <- intersect(gid.geno, gid.pheno)
	pheno     <- pheno[is.element(pheno[,1],gid),]
	M         <- M[gid,]
	N         <- length(gid)
	msg("N =",N,"individuals with phenotypic and genotypic information \n")
	
	traits <- colnames(pheno)[-1]
	msg("Detected following traits:\n",paste(traits,collapse="\n"),"\n",sep="")

	# Convert numeric to ACGT
	genoNumeric <- data.frame (Markers=colnames (M), t(M))
	SNPs        <- map [,c("Markers","Ref","Alt")]
	genoACGT    <- convertAllelesNumericToACGTFormat (genoNumeric, SNPs)
	genoImputed = data.frame (map[,-c(4,5)], genoACGT[,-1])

	# Write files
	genoFile  = paste0("out/", basename (strsplit(geno.file, split="[.]")[[1]][1]), "-IMPUTED.tbl")
	phenoFile = paste0("out/", basename (strsplit(pheno.file, split="[.]")[[1]][1]), "-IMPUTED.tbl")
	write.table (file=genoFile, genoImputed, row.names=F, quote=F, sep=",")
	write.table (file=phenoFile, pheno,  row.names=F, quote=F, sep=",")


	return (list (genof=genoFile, phenof=phenoFile, geno=M, pheno=pheno, map=map))

	# construct GWASpoly data structure
	#return(new("GWASpoly",map=map,pheno=pheno,fixed=fixed,geno=M,ploidy=ploidy))
}

read.GWASpoly2 <- function(ploidy, pheno.file, geno.file, format, n.traits, delim = ","){
	msg ("a>>>> >>>> ", geno.file)
	msg ("a>>>> >>>> ", pheno.file)
	pheno  = read.table(file=pheno.file,header=T,as.is=T,check.names=F,sep=delim)
	geno   = read.table(file=geno.file,header=T,as.is=T,check.names=F,sep=delim)
	map    = data.frame(Markers=geno[,1],Chrom=factor(geno[,2],ordered=T),Position=geno[,3],stringsAsFactors=F)

	markers <- as.matrix(geno[,-(1:3)])
	rownames(markers) <- geno[,1]
	msg ("Obtaining get/ref alleles...")
	tmp     <- apply(markers,1,get.ref)
	map$Ref <- tmp[1,]
	map$Alt <- tmp[2,]
	
	M <- apply(cbind(map$Ref,markers),1,function(x){
			y <- gregexpr(pattern=x[1],text=x[-1],fixed=T)  
			ans <- as.integer(lapply(y,function(z){ifelse(z[1]<0,ploidy,ploidy-length(z))}))	
			return(ans)
		})
	rownames (M) = colnames(geno)[-(1:3)]

	msg ("New read.GWASpoly...")
	print (M[1:10,1:10])

	fixed  = data.frame (NULL)
	# construct GWASpoly data structure
	return(new("GWASpoly",map=map,pheno=pheno,fixed=fixed,geno=M,ploidy=ploidy))
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
# Run a command string using system function and writes output to log file
#-------------------------------------------------------------
runCommand <- function (command, logFile="gwas.log") 
{
	#msg (">>>> ", command)
	errorsLog = paste0 (strsplit(logFile, split="[.]")[[1]], ".errors")
	if (logFile != "")
		system (paste0 (command, " > ", logFile," 2> ",errorsLog))
	else
		system (command)
}

#-------------------------------------------------------------
# Fast head, for debug only
#-------------------------------------------------------------
h <- function (data, n=10,m=10) {
	print (data [1:n,1:m])
}



#-------------------------------------------------------------
#-------------------------------------------------------------
read.GWASpoly <- function(ploidy, pheno.file, geno.file, format, n.traits, delim = ","){
	msg (".............read.GWaspoly...........")

	bases <- c("A","C","G","T")
	get.ref <- function(x,format) {
		y <- paste(na.omit(x),collapse="")
		ans <- apply(array(bases),1,function(z,y){length(grep(z,y,fixed=T))},y)
		if (sum(ans)>2) {stop("Error in genotype matrix: More than 2 alleles")}
		if (sum(ans)==2) {ref.alt <- bases[which(ans==1)]}
		if (sum(ans)==1) {ref.alt <- c(bases[which(ans==1)],NA)}

		return(ref.alt)
	}

	geno <- read.table(file=geno.file,header=T,as.is=T,check.names=F,sep=delim)
	map <- data.frame(Marker=geno[,1],Chrom=factor(geno[,2],ordered=T),Position=geno[,3],stringsAsFactors=F)
	markers <- as.matrix(geno[,-(1:3)])
	rownames(markers) <- geno[,1]

	tmp <- apply(markers,1,get.ref,format)
	map$Ref <- tmp[1,]
	map$Alt <- tmp[2,]
	M <- apply(cbind(map$Ref,markers),1,function(x){
		y <- gregexpr(pattern=x[1],text=x[-1],fixed=T)  
		ans <- as.integer(lapply(y,function(z){ifelse(z[1]<0,ploidy,ploidy-length(z))}))	
		return(ans)
		})

	gid.geno <- colnames(geno)[-(1:3)]
	rownames(M) <- gid.geno
	bad <- length(which(!is.element(na.omit(M),0:ploidy)))
	if (bad > 0) {stop("Invalid marker calls.")}

	impute.mode <- function(x) {
		ix <- which(is.na(x))
		if (length(ix)>0) {
			x[ix] <- as.integer(names(which.max(table(x))))
		}
		return(x)
	}

	missing <- which(is.na(M))
	if (length(missing)>0) {
		cat("Missing marker data imputed with population mode \n")
		M <- apply(M,2,impute.mode)
	}
	### matching genotypic and phenotypic data 
	pheno <- read.table(file=pheno.file,header=T,as.is=T,check.names=F,sep=delim)
	cat(paste("N =",nrow (pheno),"individuals with phenotypic and genotypic information \n"))

	fixed <- data.frame(NULL)
	traits <- colnames(pheno)[-1]
	cat(paste("Detected following traits:\n",paste(traits,collapse="\n"),"\n",sep=""))

	# construct GWASpoly data structure
	return(new("GWASpoly",map=map,pheno=pheno,fixed=fixed,geno=M,ploidy=ploidy))
}

#-------------------------------------------------------------
# Call main 
#-------------------------------------------------------------
main (args)


