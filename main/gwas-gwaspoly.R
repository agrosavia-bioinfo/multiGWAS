#!/usr/bin/Rscript
# INFO   : Script to run GWASpoly for tetraploides (modified from GWASpoly main)
# AUTHOR : Luis Garreta (lgarreta@agrosavia.co)
# DATA   : Feb/2020
# LOG    :
	# r1.1:  Added main. Changed naive kinship matrix to NULL
	# r1.02:  Hidden warnigns qchisq
	# r1.01:  Removed annotations from functions. Output results to "out/" dir 
#-------------------------------------------------------------
main <- function () {
	library (GWASpoly)
	# New class for gwaspoly
	setClass ("GWASpolyStruct", slots=c(params="list"), contains="GWASpoly.K")

	args = commandArgs(trailingOnly = TRUE)
	params = list()
	params$genotypeFile      = args [1] 
	params$phenotypeFile     = args [2] 
	params$trait             = colnames (read.csv (params$phenotypeFile))[2]
	params$ploidy            = 4
	params$gwasModel         = "naive"
	params$snpModels         = c("general","additive","1-dom", "2-dom", "diplo-general", "diplo-additive")
	params$correctionMethod  = "Bonferroni"
	params$significanceLevel = 0.05

	msgmsg("Running GWASpoly...")

	# Read input genotype and genotype (format: "numeric", "AB", or "ACGT")
	data1 <- initGWAS (params$phenotypeFile, params$genotypeFile, params$ploidy, "ACGT")

	# Control population structure
	data2 = controlPopulationStratification (data1, params$gwasModel)

	# GWAS execution
	data3 <- runGwaspoly (data2, params, 8)
	showResults (data3, params$snpModels, params$trait, params$gwasModel, 
				 params$phenotypeFile, params$ploidy)
 }

#-------------------------------------------------------------
# Set parameters and run GWASpoly tool
#-------------------------------------------------------------
runToolGwaspoly <- function (params) {
	msgmsg("Running GWASpoly GWAS...")

	scoresFile = paste0 ("out/tool-GWASpoly-scores-", params$gwasModel, ".csv")

	# Set gene action models (automatic or specific ones)
	modelsDiplo  = c ("general","additive","diplo-general", "diplo-additive", "1-dom")
	modelsTetra  = c (modelsDiplo, "2-dom")

	# Default models whether diplo or tetra 
	if (params$geneAction %in% c("all","automatic"))
		if (params$ploidy == 4) 
			snpModels = modelsTetra
		else
			snpModels = modelsDiplo
	else  # But if user specifided just one model
		if (params$geneAction == "dominant" & params$ploidy==4)
			snpModels = c ("1-dom", "2-dom")
		else if (params$geneAction == "dominant" & params$ploidy==2)
			snpModels = c ("1-dom")
		else
			snpModels = c(params$geneAction)

	msg ("SNPModels: ", snpModels)

	params = append (params, list (snpModels=snpModels))

	# Read input genotype and genotype (format: "numeric", "AB", or "ACGT")

	data1 <- initGWAS (params$phenotypeFile, params$genotypeFile, params$ploidy, "ACGT")

	# Control population structure
	data2 = controlPopulationStratification (data1, params$gwasModel)

	# GWAS execution
	data3 <- runGwaspoly (data2, params, NCORES) 
						  #params$gwasModel, params$snpModels,params$correctionMethod )

	# Show results in GWASpoly way
	msgmsg ("    >>>> GWASpoly showResults...")
	showResults (data3, params$snpModels, params$trait, params$gwasModel,params$phenotypeFile, params$ploidy)
	msg ("...Ending GWASpoly")

	# Get SNP associations
	scores  = getQTLGWASpoly (data3, params$gwasModel, params$ploidy)
	colnames (scores)[colnames(scores) %in% c("Chrom","Position")] = c ("CHR","POS")

	scoresColumns = c("MODEL", "GC", "Marker", "CHR", "POS", "P", "SCORE", "THRESHOLD", "DIFF")
	scores <- data.frame (scores[,scoresColumns], scores [,setdiff (colnames(scores), scoresColumns)])

	write.table (file=scoresFile, scores, quote=F, sep="\t", row.names=F)

	return (list (tool="GWASpoly", scoresFile=scoresFile, scores=scores))
}

#-------------------------------------------------------------
# Control population structure using default Kinship and PCs
#-------------------------------------------------------------
controlPopulationStratification <- function (data1, gwasModel) {
	msgmsg ();msgmsg("Controlling populations structure...")

	if (gwasModel=="naive") {
		msgmsg("    >>>> Without any correction") 
		markers       = data1@pheno [,1]
		n             = length (markers)
		kinshipMatrix = matrix (diag (n), n, n, dimnames=list (markers, markers))
		#dataTmp      <- set.K (data1, K=kinshipMatrix)
		dataTmp       = set.K (data1, K=NULL)
		data2         = new ("GWASpolyStruct", dataTmp)
	}else if (gwasModel == "full") {
		msgmsg("    >>>> Using default Kinship and PCs=5 ")
		kinshipMatrix = NULL
		dataTmp       = set.K (data1)
		write.csv (dataTmp@K,"out/GWASpoly-kinship-matrix.csv", quote=F, row.names=T)
		data2         = new ("GWASpolyStruct", dataTmp)
		data2@params  = set.params (n.PC=5, fixed=NULL, fixed.type=NULL)
	}else 
		stop ("Unknown ploidy number (only 2 or 4)")

	return (data2)
}

#-------------------------------------------------------------
# GWAS execution
#-------------------------------------------------------------
runGwaspoly <- function (data2, params, NCORES) {
	gwasModel        = params$gwasModel
	snpModels        = params$snpModels
	correctionMethod = params$correctionMethod
	signLevel        = params$significanceLevel

 	if (gwasModel %in% c("naive")) {
		msgmsg(">>>> Without params")
		data3 = GWASpoly(data2, models=snpModels, traits=NULL, params=NULL, n.core=NCORES)
	}else {
		msgmsg(">>>> With params")
		data3 = GWASpoly(data2, models=snpModels, traits=NULL, params=data2@params)
	}
	
	# QTL Detection
	data4 = set.threshold (data3, method=correctionMethod,level=signLevel,n.core=NCORES)

	return (data4)
}

#-------------------------------------------------------------
# Plot results
#-------------------------------------------------------------
showResults <- function (data3, models, trait, gwasModel, phenotypeFile, ploidy) {
	msgmsg ();msgmsg("Writing GWASpoly results...")
	#outFile       = paste0 ("out/tool-GWASpoly-scores-", gwasModel)
	#scoresFile    = paste0 (outFile,".csv")
	#scoresFileAll = paste0 (outFile,"-all.csv")
	plotFile      = paste0 ("out/out-GWASpoly-", gwasModel, "-plots.pdf") 

	# QTL Detection
	#data5 = set.threshold (data3, method=correctionMethod,level=0.05,n.core=4)

	# Plot results	
	plotMahattanQQ (plotFile, models, data3, trait, data3, gwasModel, ploidy) 

	#msgmsg(">>>> Writing QTLs to file: ", scoresFile, "...")
	#write.GWASpoly (data5, trait, paste0(scoresFile,".qtls"), "scores", delim="\t")

	#scoresTableAll  = getQTLGWASpoly (data3, gwasModel, ploidy)
	#write.table (file=scoresFileAll, scoresTableAll, quote=F, sep="\t", row.names=F)

	#scoresTableSorted = scoresTableAll [order (scoresTableAll$GC, scoresTableAll$DIFF, decreasing=T),]
	#scoresTable       = scoresTableSorted [!duplicated (scoresTableSorted$Marker, fromLast=F),]
	#write.table (file=scoresFile, scoresTable, quote=F, sep="\t", row.names=F)
}

#-------------------------------------------------------------
# Manhattan and QQ plots
#-------------------------------------------------------------
plotMahattanQQ <- function (plotFile, models, data5, trait, data3, gwasModel, ploidy) {
	# Create test models for each ref|alt allele if dominant models present
	testModels = c()
	for (m in models)
		if (m=="1-dom") testModels = c (testModels, "1-dom-alt", "1-dom-ref")
		else if (m=="2-dom") testModels = c (testModels, "2-dom-alt", "2-dom-ref")
		else testModels = c (testModels, m)

	n = length (testModels)

	pdf (file=plotFile, width=11, height=15)
	op <- par(mfrow = c(n,2), mar=c(3.5,3.5,2,1), oma=c(0,0,0,0), mgp = c(2.2,1,0)) #MultiGWAS tools
	for (i in 1:n) {
		manhattan.plot (data5, trait=trait, model=testModels [i])
		qqPlot(data3,trait=trait, model=testModels[i], cex=0.3)
	}

	plotTitle = sprintf ("%s gwas %s-ploidy for %s trait", gwasModel, ploidy, trait)  
	mtext(plotTitle, outer=T,  cex=1.5,  line=0)
	par(op)
	dev.off()
}

#-------------------------------------------------------------
# Extracts significant QTL
#-------------------------------------------------------------
getQTLGWASpoly <- function(data,gwasModel, ploidy, traits=NULL,models=NULL) {
	stopifnot(inherits(data,"GWASpoly.thresh"))

	if (is.null(traits)) traits <- names(data@scores)
	else stopifnot(is.element(traits,names(data@scores)))

	if (is.null(models)) models <- colnames(data@scores[[1]])
	else stopifnot(is.element(models,colnames(data@scores[[1]])))

	n.model <- length(models)
	n.trait <- length(traits)
	output <- data.frame(NULL)
	for (j in 1:n.model) {
		#ix <- which(data@scores[[traits[1]]][,models[j]] > (data@threshold[traits[1],models[j]]) - 1)
		ix <- which (data@scores[[traits[1]]][,models[j]] != 0)
		markers <-  data.frame (SNP=data@map[ix,c("Marker")])

		scores <- data@scores[[1]][,models[j]]
		datax = calculateInflationFactor (scores)

		n.ix <- length(ix)
		
		gc=rep(datax$delta,n.ix) 
		scores=round(data@scores[[traits[1]]][ix,models[j]],2)
		thresholds=round(rep(data@threshold[traits[1],models[j]],n.ix),2)
		diffs = (scores - thresholds)
		pvalues = 10^(-scores)
		df = data.frame(Ploidy=rep (ploidy, n.ix), Type=rep (gwasModel, n.ix),
						data@map[ix,], GC=gc, MODEL=rep(models[j],n.ix),
						P=pvalues,SCORE=scores, THRESHOLD=thresholds, DIFF=diffs,
						Effect=round(data@effects[[traits[1]]][ix,models[j]],2))
						#stringsAsFactors=F,check.names=F)

		output <- rbind(output, df)
	}
	#out <-cbind (Type=gwasModel, output)
	#output <- output [order(-output$GC, -output$DIFF),]
	output = output [order(-output$DIFF),]
	#output = output [!duplicated (output$Marker),]
	#outputPositives = output [output$DIFF > 0,]
	#outputNegatives = output [output$DIFF <= 0,]

	#outQTLsAllSNPs = rbind (outputPositives, outputNegatives)

	return(output)
}

#-------------------------------------------------------------
# Calculate the inflation factor from -log10 values
# It can fire warning, here they are hidign
#-------------------------------------------------------------
calculateInflationFactor <- function (scores) {
	oldw <- getOption("warn")
	options(warn = -1)

	remove <- which(is.na(scores))
	if (length(remove)>0) 
		x <- sort(scores[-remove],decreasing=TRUE)
	else 
		x <- sort(scores,decreasing=TRUE)

	pvalues = 10^-x
	chisq <- na.omit (qchisq(1-pvalues,1))
	delta  = round (median(chisq)/qchisq(0.5,1), 3)

	options (warn = oldw)

	return (list(delta=delta, scores=x))
}

#-------------------------------------------------------------
# QQ plot
#-------------------------------------------------------------
qqPlot <- function(data,trait,model,cex=1,filename=NULL) {
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
# Read input genotype and genotype (format: "numeric", "AB", or "ACGT")
#-------------------------------------------------------------
initGWAS <- function (phenotypeFile, genotypeFile, ploidy, format="ACGT") {
	msgmsg ();msgmsg("Initializing GWAS...");msgmsg ()

	data1 <- read.GWASpoly (ploidy = ploidy, pheno.file = phenotypeFile, 
							geno.file = genotypeFile, format = "ACGT", n.traits = 1, delim=",")

	return (data1)
}
#-------------------------------------------------------------
# Add label to filename
#-------------------------------------------------------------
addLabel <- function (filename, label)  {
	nameext = strsplit (filename, split="[.]")
	newName = paste0 (nameext [[1]][1], "-", label, ".", nameext [[1]][2])
	return (newName)
}
#-------------------------------------------------------------
#-------------------------------------------------------------
msgmsg <- function (...) {
  messages = unlist (list (...))
  cat ("\t>>>>", messages, "\n")
}

#-------------------------------------------------------------
#main ()
