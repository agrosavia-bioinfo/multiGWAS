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
# Control population structure using default Kinship and PCs
#-------------------------------------------------------------
controlPopulationStratification <- function (data1, gwasModel, data2) 
{
	msg();msg("Controlling populations structure...")

 	# Load data instead calculate it
	if (!is.null (data2)) {msg (">>>> Loading kinship..."); return (data2) }

	if (gwasModel=="Naive") {
		msg ("    >>>> Without any correction") 
		markers = data1@pheno [,1]
		n       = length (markers)
		kinshipMatrix = matrix (diag (n), n, n, dimnames=list (markers, markers))
		dataTmp      <- set.K (data1, K=kinshipMatrix)
		#dataTmp      <- set.K (data1, K=NULL)
		data2        = new ("GWASpolyStruct", dataTmp)
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
showResults <- function (data3, testModels, trait, gwasModel, correctionMethod, phenotypeFile, snpsAnnFile, ploidy, qtlsFile) 
{
	msg();msg ("Writing results...", trait)

	msg (">>>> Plotting results...")
	phenoName = strsplit (phenotypeFile, split=".scores")[[1]][1]
	plotName = sprintf("out/out-Gwasp%s-%s-plots.pdf", ploidy, gwasModel)

	n = length (testModels)
	
	# QTL Detection
	data5 = set.threshold (data3, method=correctionMethod,level=0.05,n.core=4)

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
	#write.GWASpoly (data5, trait, qtlsFile, "scores", delim="\t")

	significativeQTLs  = getQTL (data5, snpsAnnFile, gwasModel, ploidy)
	#significativesFile = addLabel (qtlsFile, "SIGNIFICATIVES")
	write.table (file=qtlsFile, significativeQTLs, quote=F, sep="\t", row.names=F)

}


#-------------------------------------------------------------
# Extracts significant QTL
#-------------------------------------------------------------
getQTL <- function(data,snpsAnnFile, gwasModel, ploidy, traits=NULL,models=NULL) 
{
	stopifnot(inherits(data,"GWASpoly.thresh"))

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
		#ix <- which(data@scores[[traits[1]]][,models[j]] > (data@threshold[traits[1],models[j]]) - 1)
		ix <- which (data@scores[[traits[1]]][,models[j]] != 0)
		markers <-  data.frame (SNP=data@map[ix,c("Marker")])

		if (!is.null (snpsAnnFile)) 
			snpAnn  <- merge (markers, snpsAnnotations, by.x="SNP",by.y="SNP_id", sort=F)[,c(2,7)]
		else
			snpAnn = "None"

		scores <- data@scores[[1]][,models[j]]
		datax = calculateInflationFactor (scores)

		n.ix <- length(ix)
		
		gc=rep(datax$delta,n.ix) 
		scores=round(data@scores[[traits[1]]][ix,models[j]],2)
		thresholds=round(rep(data@threshold[traits[1],models[j]],n.ix),2)
		diffs = (scores - thresholds)
		pvalues = 10^(-scores)
		df = data.frame(Ploidy=rep (ploidy, n.ix),
						Type=rep (gwasModel, n.ix),
						data@map[ix,],
						GC=gc,
						Model=rep(models[j],n.ix),
						P=pvalues,SCORE=scores, THRESHOLD=thresholds, DIFF=diffs,
						Effect=round(data@effects[[traits[1]]][ix,models[j]],2))
						#snpAnn) 
						#stringsAsFactors=F,check.names=F)

		output <- rbind(output, df)
	}
	#out <-cbind (Type=gwasModel, output)
	output <- output [order(-output$GC, -output$DIFF),]
	output = output [!duplicated (output$Marker),]
	outputPositives = output [output$DIFF > 0,]
	outputNegatives = output [output$DIFF <= 0,]

	outputTotal = rbind (outputPositives, outputNegatives)
	return(outputTotal)
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
# Read input genotype and genotype (format: "numeric", "AB", or "ACGT")
#-------------------------------------------------------------
initGWAS <- function (phenotypeFile, genotypeFile, ploidy, format="ACGT", data1) 
{
	msg();msg ("Initializing GWAS...");msg()
	# When data is previously loaded
	if (!is.null (data)) {msg(">>>> Loading GWAS data..."); return (data)}

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
