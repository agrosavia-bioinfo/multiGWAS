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
#args = c("in/config-test500-Structure-noFilters.config")

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
source ("gwas-formats.R")           # Module with functions to convert between different genotype formats 
source ("gwas-summary.R")           # Module with functions to create summaries: tables and venn diagrams

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
	genotypeFile  <<- config$genotypeFile
	phenotypeFile <<- config$phenotypeFile

	system ("mkdir out")
	system (sprintf ("cp %s out/", configFile))

	msg (">>>>>>>>>>>>", config$gwasModel, "<<<<<<<<<<<") 

	# Read, filter, and check phenotype and genotype
	data <- dataPreprocessing (genotypeFile, phenotypeFile, config)

	config$trait = data$trait

	# Run the four tools in parallel

	mclapply (c("Gwasp", "Plink", "Shesis", "Tassel"), runGwasTool, config, mc.cores=4)

	# Create outputs: tables, figures
	markersSummaryTable ("out/", config$gwasModel, "out/", nBEST=5)

	# Save test in new dir
	outName = gsub ("config","test", configFile)
	system (sprintf ("mv out %s", outName))
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
	data1 <- initGWAS (phenotypeFile, genotypeFile, ploidy, params$genotypeFormat, data1)

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
	# First filtering: imputation, MAF, common names

	imputedFiles <-filterByImputationMAFNames (phenotypeFile, genotypeFile, config)
	genotypeFile  = imputedFiles$geno
	phenotypeFile = imputedFiles$pheno

	# Format convertion from gwasp4 to plink2
	msg();msg ("Converting gwaspoly to plink formats...")
	markersIdsMap = gwaspToPlinkGenoMap (genotypeFile)
	plinkFile = gwaspTetraGenoToPlinkPed (genotypeFile, markersIdsMap)

	# Apply filters to genotype (markers and samples) by calling external program
	filtered = filterByMissingHWE (plinkFile, config)

	# Filter phenotype
	phenoAll = read.csv (phenotypeFile, header=T, sep=",", check.names=F)
	phenoFiltered  = phenoAll [phenoAll[,1] %in% filtered$Samples,]

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
	gwasp2shesisGenoPheno ("out/filtered-gwasp4-genotype.tbl", "out/filtered-gwasp4-phenotype.tbl") 

	return (list (genotypeFile=outGenoFile, phenotypeFile=outPhenoFile, trait=trait))

}
#-------------------------------------------------------------
# Apply filters to genotype (markers and samples) by calling external program
#-------------------------------------------------------------
filterByMissingHWE <- function (plinkFile, config) 
{
	cmm = paste ("plink --file", plinkFile, "--make-bed", "--out", paste0(plinkFile,"-QC"))
	if (config$filtering==T) {
		# Filter missingness per sample (MIND)"
		if (!is.null(config$MIND)) cmm=paste (cmm, paste ("--mind", config$MIND))
		# Filter missingness per SNP    (GENO)
		if (!is.null(config$GENO)) cmm=paste (cmm, paste ("--geno", config$GENO))
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
# Control population structure using default Kinship and PCs
#-------------------------------------------------------------
controlPopulationStratification <- function (data1, gwasModel, data2) 
{
	msg();msg("Controlling populations structure...")

 	# Load data instead calculate it
	if (!is.null (data2)) {msg (">>>> Loading kinship..."); return (data2) }

	if (gwasModel=="Naive") {
		msg ("    >>>> Without any correction") 
		markerNames  = data1@pheno [,1]
		dataTmp      = set.K (data1, K=NULL)
		data2 <- new ("GWASpolyStruct", dataTmp)
	}else if (gwasModel == "Structure") {
		msg ("    >>>> Using default Kinship and PCs=5 ")
		kinshipMatrix = NULL
		dataTmp       = set.K (data1)
		data2         = new ("GWASpolyStruct", dataTmp)
		data2@params  = set.params (n.PC=5, fixed=NULL, fixed.type=NULL)
	}else 
		stop ("unknown ploidy number (only 2 or 4)")

	return (data2)
}

#-------------------------------------------------------------
# GWAS execution
#-------------------------------------------------------------
runGwaspoly <- function (data2, gwasModel, snpModels, data3) 
{
	msg();msg("Running GWASpoly...")

	if (!is.null (data3)) { msg (">>>> Loading GWASpoly..."); return (data3) }

 	if (gwasModel %in% c("Naive","Kinship")) {
		msg (">>>> Without params")
		data3 = GWASpoly(data2, models=snpModels, traits=NULL, params=NULL, n.core=4)
	}else {
		msg (">>>> With params")
		data3 = GWASpoly(data2, models=snpModels, traits=NULL, params=data2@params)
	}
	
	return (data3)
}

#-------------------------------------------------------------
# Plot results
#-------------------------------------------------------------
showResults <- function (data3, testModels, trait, gwasModel, phenotypeFile, snpsAnnFile, ploidy) 
{
	msg();msg ("Writing results...", trait)

	msg (">>>> Plotting results...")
	phenoName = strsplit (phenotypeFile, split=".scores")[[1]][1]
	plotName = sprintf("out/out-Gwasp%s-%s-plots.pdf", ploidy, gwasModel)

	n = length (testModels)
	
	# QTL Detection
	data5 = set.threshold (data3, method="Bonferroni",level=0.05,n.permute=100,n.core=4)

	# Plots
	pdf (file=plotName, width=11, height=7)
	# QQ plot 
	op <- par(mfrow = c(2,n), oma=c(0,0,3,0))
	for (i in 1:length(testModels)) {
		#par (cex.main=0.5, cex.lab=0.5, cex.axis=0.5, ann=T)
		qqPlot(data3,trait=trait, model=testModels[i], cex=0.3)
	}

	# Manhattan plot 
	for (i in 1:length(testModels)) {
		#par (cex=1.5)
		#manhattan.plot (y.max=20,data5, trait=trait, model=testModels [i])
		manhattan.plot (data5, trait=trait, model=testModels [i])
	}
	plotTitle = sprintf ("%s gwas %s-ploidy for %s trait", gwasModel, ploidy, trait)  
	mtext(plotTitle, outer=T,  cex=1.5,  line=0)
	par(op)
	dev.off()

	msg (">>>> Writing QTLs...")
	significativeQTLs = getQTL (data5, snpsAnnFile, gwasModel, ploidy)
	outFile = sprintf ("out/out-Gwasp%s-%s-QTLs.scores", ploidy, gwasModel) 
	write.table (file=outFile, significativeQTLs, quote=F, sep="\t", row.names=F)
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

	if (!is.null (snpsAnnFile)) {
		msg (">>>> Reading associations...")
		snpsAnnotations <- read.csv (file=snpsAnnFile, header=T)
	}

	n.model <- length(models)
	n.trait <- length(traits)
	output <- data.frame(NULL)
	for (j in 1:n.model) {
		ix <- which(data@scores[[traits[1]]][,models[j]] > data@threshold[traits[1],models[j]])
		markers <-  data.frame (SNP=data@map[ix,c("Marker")])
		msg ("         ", models [j],": ", as.character (markers[,1]))

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
	sc <<- data@scores
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
# Get params from config file and define models according to ploidy
#-------------------------------------------------------------
getConfigurationParameters <- function (configFile) 
{
	msg("Reading config file...")

	params = config::get (file=configFile) 

	message ("------------------------------------------------")
	message ("Summary of input parameters:")
	message ("------------------------------------------------")
	msg ("Phenotype filename : ", params$phenotypeFile) 
	msg ("Genotype filename  : ", params$genotypeFile) 
	msg ("GwAS model         : ", params$gwasModel) 
	msg ("Genotype format    : ", params$format) 
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

	geno    <- read.table(file=geno.file,header=T,as.is=T,check.names=F,sep=delim)
	print (dim(geno))
	map     <- data.frame(Markers=geno[,1],Chrom=factor(geno[,2],ordered=T),Position=geno[,3],stringsAsFactors=F)
	markers <- as.matrix(geno[,-(1:3)])
	rownames(markers) <- geno[,1]
	
	tmp     <- apply(markers,1,get.ref)
	map$Ref <- tmp[1,]
	map$Alt <- tmp[2,]
	M <- apply(cbind(map$Ref,markers),1,function(x){
			y <- gregexpr(pattern=x[1],text=x[-1],fixed=T)  
			ans <- as.integer(lapply(y,function(z){ifelse(z[1]<0,ploidy,ploidy-length(z))}))	
			return(ans)
		})


	gid.geno <<- colnames(geno)[-(1:3)]
	rownames(M) <- gid.geno


	bad <- length(which(!is.element(na.omit(M),0:ploidy)))
	if (bad > 0) {stop("Invalid marker calls.")}
	
	MAFTHRESHOLD=0.0
	if (config$filtering == T) MAFTHRESHOLD = config$MAF
	MAF <- apply(M,2,function(x){AF <- mean(x,na.rm=T)/ploidy;MAF <- ifelse(AF > 0.5,1-AF,AF)})
	polymorphic <- which(MAF>=MAFTHRESHOLD)

	M <- M[,polymorphic]
	map <- map[polymorphic,]
	#map <- map[order(map$Chrom,map$Position),]
	M <- M[,map$Marker]
	m <- nrow(map)
	cat(paste("Number of polymorphic markers:",m,"\n"))
	
	if (config$imputation==T) {
		missing <- which(is.na(M))
		if (length(missing)>0) {
			cat("Missing marker data imputed with population mode \n")
			M <- apply(M,2,impute.mode)
		}
	}
	
	### matching genotypic and phenotypic data 
	pheno <- read.table(file=pheno.file,header=T,as.is=T,check.names=F,sep=delim)
	gid.pheno <<- unique(pheno[,1])
	gid <<- intersect(gid.geno, gid.pheno)
	pheno <- pheno[is.element(pheno[,1],gid),]
	M <- M[gid,]
	N <- length(gid)
	cat(paste("N =",N,"individuals with phenotypic and genotypic information \n"))
	
	n.fixed <- ncol(pheno) - n.traits - 1
	if (n.fixed > 0) {
		fixed <- data.frame(pheno[,(n.traits+2):ncol(pheno)],stringsAsFactors=F)
		fixed.names <- colnames(pheno)[(n.traits+2):ncol(pheno)]
		colnames(fixed) <- fixed.names
		pheno <- data.frame(pheno[,1:(1+n.traits)],stringsAsFactors=F)
		cat(paste("Detected following fixed effects:\n",paste(fixed.names,collapse="\n"),"\n",sep=""))
	} else 
		fixed <- data.frame(NULL)
	
	traits <- colnames(pheno)[-1]
	cat(paste("Detected following traits:\n",paste(traits,collapse="\n"),"\n",sep=""))
	

	# Convert numeric to ACGT
	genoNumeric <<- data.frame (Markers=colnames (M), t(M))
	SNPs        <- map [,c("Markers","Ref","Alt")]
	genoACGT    <<- numericGenotypeToACGTFormatFromData (genoNumeric, SNPs)
	genoImputed = data.frame (map[,-c(4,5)], genoACGT[,-1])

	# Write files
	genoFile  = paste0("out/", basename (strsplit(geno.file, split="[.]")[[1]][1]), "-IMPUTED.tbl")
	phenoFile = paste0("out/", basename (strsplit(pheno.file, split="[.]")[[1]][1]), "-IMPUTED.tbl")
	write.table (file=genoFile, genoImputed, row.names=F, quote=F, sep=",")
	write.table (file=phenoFile, pheno,  row.names=F, quote=F, sep=",")

	return (list (geno=genoFile, pheno=phenoFile))

	# construct GWASpoly data structure
	#return(new("GWASpoly",map=map,pheno=pheno,fixed=fixed,geno=M,ploidy=ploidy))
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
# Call main 
#-------------------------------------------------------------
main (args)


