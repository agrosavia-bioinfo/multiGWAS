#!/usr/bin/Rscript

# INFO   : Create different summarized reports (tables, plots, Venn diagrams) for multiGWAS tool analysis
# AUTHOR : Luis Garreta (lgarreta@agrosavia.co)
# DATA   : feb/2020
# LOG: 
#	r4.1: Fixed transparency and axis error (options (bitmapType="cairo"))
#	r4.0: Modified to work with markdown, but better to only report outputs (PNGs) to be included by markdown
#	r3.0: Manhattan and QQ plots. Formated to create a Markdown report (not yet)
#	r2.0: Improve selection of data from tables. Titles to graphics and files
#	r1.0: Message to log files
#	r0.9: Full working with funtions: create snpTables and ven diagrams using parameters
#	r0.8: Create venn diagrams, summary table of first Ns"
#     


suppressMessages (library (dplyr))
suppressMessages (library (qqman))
suppressMessages (library (VennDiagram))
suppressMessages (library (config))  # For read config file

options (bitmapType="cairo", width=300)
#options(scipen=999)

#-------------------------------------------------------------
# Main function
# Input files are taken from input dir
# Outupt are written to output dir
#-------------------------------------------------------------
main <- function () {

	msg ("Main...")

	inputDir    = "out/"
	outputDir   = "report/"
	gwasModel    = "Full"
	reportTitle = "GWAS Report"
	nBest = 6

	createReports (inputDir, gwasModel, outputDir, outputDir)
}

#-------------------------------------------------------------
# Function to create different reports:
#	1- 1 table of best SNPs
#	2- 1 table of significative SNPs
#   3- 1 Venn diagram of best SNPs
#   4- 1 Venn diagram of significative SNPs
#	5- 1 multiplot of 4x4 manhattan and QQ plots
#-------------------------------------------------------------
createReports <- function (inputDir, genotypeFile, phenotypeFile, gwasModel, outputDir, nBest=7) 
{
	msg ("Creating reports for ", gwasModel, "...")
	createDir (outputDir)

	msg ("Writing table input config parameters...")
	writeConfigurationParameters (inputDir, outputDir)

	msg ("Writing table with summary results...")
	snpTables = markersSummaryTable (inputDir, gwasModel, outputDir,  nBest)

	msg ("Writing table with ", nBest, " best ranked SNPs Table...")
	outName = paste0(outputDir, "/out-multiGWAS-scoresTable-best.scores")
	write.table (file=outName, snpTables$best, row.names=F,quote=F, sep="\t")

	msg ("Writing table with significative SNPs...")
	outName = paste0(outputDir, "/out-multiGWAS-scoresTable-significatives.scores")
	write.table (file=outName, snpTables$significatives, row.names=F,quote=F, sep="\t")

	msg ("Writing Venn diagram with best SNPs...")
	png (paste0 (outputDir,"/out-multiGWAS-vennDiagram-best.png"), res=72)
	commonBest = markersVennDiagrams (snpTables$best, gwasModel, "Best")
	dev.off()
	pdf (paste0 (outputDir,"/out-multiGWAS-vennDiagram-best.pdf"))
	commonBest = markersVennDiagrams (snpTables$best, gwasModel, "Best")
	dev.off()

	msg ("Writing Venn diagram with significative SNPs...")
	outFilename = 
	png (paste0 (outputDir,"/out-multiGWAS-vennDiagram-significatives.png"), res=72)
	commonSign = markersVennDiagrams (snpTables$significatives, gwasModel, "Significatives")
	dev.off ()
	pdf (paste0 (outputDir,"/out-multiGWAS-vennDiagram-significatives.pdf"))
	commonSign = markersVennDiagrams (snpTables$significatives, gwasModel, "Significatives")
	dev.off ()

	msg ("Writing Manhattan and QQ plots...")
	png (paste0 (outputDir, "/out-multiGWAS-manhattanQQ-plots.png"), width=11, height=15, units="in", res=120)
	op=markersManhattanPlots (inputDir, gwasModel, commonBest, commonSign, snpTables$significatives, outputDir)
	dev.off()
	pdf (paste0 (outputDir, "/out-multiGWAS-manhattanQQ-plots.pdf"), width=11, height=15)
	op=markersManhattanPlots (inputDir, gwasModel, commonBest, commonSign, snpTables$significatives, outputDir)
	par (op)
	dev.off()

	msg ("Creating SNP heatmaps for the 4 best ranked SNPs...")
	genoNumericFilename = ACGTToNumericGenotypeFormat (genotypeFile)
	snpList = snpTables$best$SNP [1:4]
	createHeatmapForSNPList (outputDir, genotypeFile, genoNumericFilename, phenotypeFile, commonSign)
}

#-------------------------------------------------------------
# Calculate the inflation factor from -log10 values
#-------------------------------------------------------------
## @knitr calculateInflationFactor
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

#------------------------------------------------------------------------
#------------------------------------------------------------------------
## @knitr markersManhattanPlots
markersManhattanPlots <- function (inputDir, gwasModel, commonBest, commonSign, summarySignificatives, outputDir, nBest=8) {
	files =  list.files(inputDir, pattern=paste0("^(.*(",gwasModel,").*(scores)[^$]*)$"), full.names=T)
	#pdf (paste0 (outputDir, "/out-summary.manhattan-qq-plots.pdf"), width=11, height=7)
	op <- par(mfrow = c(4,2), mar=c(3.5,3.5,3,1), oma=c(0,0,0,0), mgp = c(2.2,1,0))
	for (filename in files) {
		data           = na.omit (read.table (file=filename, header=T))

		bestThresholdScore = data [nBest,c("SCORE")]
		bestThreshold      = 10^-bestThresholdScore

		names          = unlist (strsplit (basename (filename), "[-|.]"))
		mainTitle      = paste0 (names[2],"-", names [3])

		if (grepl ("GWASpoly", filename)) {
			tool = "GWASpoly"
			data = data [data[,"Model"] %in% "additive",]
			gwasResults = data.frame (SNP=data$Marker, CHR=data$Chrom, BP=data$Position, P=10^-data$SCORE)
		}
		else if (grepl ("SHEsis", filename)) {
			tool = "SHEsis"
			gwasResults = data.frame (SNP=data$SNP, CHR=data$CHR, BP=data$POS, P=10^-data$SCORE)
		}
		else if (grepl ("PLINK", filename)) {
			tool = "PLINK"
			gwasResults = data.frame (SNP=data$SNP, CHR=data$CHR, BP=data$POS, P=10^-data$SCORE)
		}
		else if (grepl ("TASSEL", filename)) {
			tool = "TASSEL"
			gwasResults = data.frame (SNP=data$Marker, CHR=data$Chr, BP=data$Pos, P=10^-data$SCORE)
		}
		ss = summarySignificatives
		if (tool %in% ss$TOOL)
			signThresholdScore = min (ss [ss$TOOL==tool,"SCORE"])
		else
			signThresholdScore = 0.95*ceiling (data[1, "SCORE"]) 

		colorsBlueOrange = c("blue4", "orange3")
		manhattan(gwasResults,col = c("gray10", "gray60"), highlight=intersect (commonBest, gwasResults$SNP), annotatePval=bestThreshold, annotateTop=F,
				  suggestiveline=bestThresholdScore, genomewideline=signThresholdScore, main=mainTitle, logp=T, cex=2)

		text (x=0, y=signThresholdScore*1.02, "Significative",, col="red", pos=4)
		text (x=0, y=bestThresholdScore*0.92, "Best",, col="blue", pos=4)

		datax = calculateInflationFactor (-log10 (gwasResults$P))
		qq (gwasResults$P)
		mtext (bquote(lambda[GC] == .(datax$delta)), side=3, line=-2, cex=0.7)
		#title (bquote(lambda[GC] == .(datax$delta)))
	}
	#par (op)
	#dev.off()
	return (op)
}

#------------------------------------------------------------------------
# Create Venn diagram of common markers using info from summary table
#------------------------------------------------------------------------
## @knitr markersVennDiagrams
markersVennDiagrams <- function (summaryTable, gwasModel, scoresType){
	flog.threshold(ERROR)
	x <- list()
	x$GWASpoly = summaryTable %>% filter (TOOL %in% "GWASpoly") %>% select (SNP) %>% .$SNP
	x$SHEsis   = summaryTable %>% filter (TOOL %in% "SHEsis") %>% select (SNP) %>% .$SNP
	x$PLINK    = summaryTable %>% filter (TOOL %in% "PLINK")  %>% select (SNP) %>% .$SNP
	x$TASSEL   = summaryTable %>% filter (TOOL %in% "TASSEL") %>% select (SNP) %>% .$SNP

	a = intersect (x$GWASpoly, x$SHEsis)
	b = intersect (x$GWASpoly, x$PLINK)
	c = intersect (x$GWASpoly, x$TASSEL)
	d = intersect (x$SHEsis,   x$PLINK)
	e = intersect (x$SHEsis,   x$TASSEL)
	f = intersect (x$PLINK,    x$TASSEL)
	#commonSNPs = intersect (intersect (x$GWASpoly, x$SHEsis), intersect(x$Plink, x$Tassel))
	commonSNPs = union (union (a,b),union (union (c,d),union (e,f)))

	mainTitle = paste0(gwasModel, "-", scoresType)
	msg();msg (mainTitle);msg()
	v0 <- venn.diagram(x, height=12000, width=12000, alpha = 0.5, filename = NULL, # main=mainTitle,
						col = c("red", "blue", "green", "yellow"), cex=0.9, margin=0.0,
						fill = c("red", "blue", "green", "yellow")) 

	overlaps <- calculate.overlap(x)
	overlaps <- rev(overlaps)

	posOverlap = as.numeric (gsub ("a","", (names (overlaps))))
	for (i in 1:length(overlaps)){
		pos = posOverlap [i]
		v0[[pos+8]]$label <- paste(overlaps[[i]], collapse = "\n")
 	}
 
	#pdf(paste0(outputDir,"/out-summary-gwas-venndiagram-", scoresType, ".pdf"))
	grid.draw(v0)
	#dev.off()

	return (commonSNPs)
}

#------------------------------------------------------------------------
# Create a summary table of best and significative markers
#------------------------------------------------------------------------
markersSummaryTable <- function (inputDir, gwasModel, outputDir="out", nBest=5) {

	map = read.table (file=paste0(inputDir,"/map.tbl"))
	rownames (map) = map [,1]

	files =  list.files(inputDir, pattern=paste0("^(.*(",gwasModel,").*(scores)[^$]*)$"), full.names=T)
	#msg ("Files: ", inputDir); print (files)
	msg ("CREATING THE SUMMARY...")
	summaryTable = data.frame ()

	tool=""
	for (f in files) {
		msg ("Processing file: ", f)
		data <- read.table (file=f, header=T)
		#if (nrow(data)>nBest) data=data [1:nBest,] 
		pVal	<- data$P
		pscores <- data$SCORE
		tscores <- data$THRESHOLD
		signf   = pscores >= tscores

		flagNewData = F
		if (grepl("GWASpoly", f)) {
			tool    = "GWASpoly"
			snps    <- data$Marker
			chrom   <- data$Chrom
			pos	    <- data$Position
			flagNewData = T
		}else if (grepl ("PLINK", f)) {
			tool    = "PLINK"
			snps    = data$SNP
			chrom   = data$CHR
			pos	    = map [snps, "Position"]
			flagNewData = T
		}else if (grepl ("TASSEL", f)) {
			tool    = "TASSEL"
			snps    = data$Marker
			chrom   = data$Chr
			pos		= data$Pos
			flagNewData = T
		}else if (grepl ("SHEsis", f)) {
			tool    = "SHEsis"
			snps    = data$SNP
			chrom	= map [snps, "Chrom"]
			pos	    = map [snps, "Position"]
			flagNewData = T
		}
		if (flagNewData==T) {
			dfm = data.frame (TOOL=tool, MODEL=gwasModel, CHROM=chrom, POSITION=pos, SNP=snps, 
							  PVALUE = round (pVal,6), SCORE=pscores, THRESHOLD=tscores, SIGNIFICANCE=signf )
			dfm = dfm %>% distinct (SNP, .keep_all=T)
			if (nrow(dfm)>nBest) dfm=dfm [1:nBest,] 
			summaryTable <- rbind (summaryTable, dfm)
			flagNewData = F
		}
	}

	summaryTable = summaryTable [which(!is.na(summaryTable$SIGNIFICANCE)),]
	summarySignificatives = summaryTable %>% filter (SIGNIFICANCE%in%T) 

	return (list (best=summaryTable, significatives=summarySignificatives))
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
#-------------------------------------------------------------
# Print a log message with the parameter
#-------------------------------------------------------------
## @knitr msg
msg <- function (...) 
{
		messages = unlist (list (...))
		cat (">>>>", messages, "\n")
}

#-------------------------------------------------------------
# Add label to filename
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
# Get params from config file and define models according to ploidy
#-------------------------------------------------------------
## @knitr writeConfigurationParameters
writeConfigurationParameters <- function (inputDir, outputDir) 
{
	msg("Reading config file...", inputDir)
	configFile = paste0(inputDir, list.files (inputDir, pattern="config")[1])
	msg (configFile)

	params = config::get (file=configFile) 

	paramsDF = data.frame (PARAMETER=character(), VALUE=character ())
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="Genotype filename", VALUE=toString (params$genotypeFile)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="Phenotype filename", VALUE=toString (params$phenotypeFile)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="Significance level (Genome-wide significance level)", VALUE=toString (params$significanceLevel)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="Correction method (Bonferroni or FDR)", VALUE=toString (params$correctionMethod)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="GWAS model (Full or Naive)", VALUE=toString (params$gwasModel)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="nBest (Number of best-ranked SNPs to be reported)", VALUE=toString (params$nBest)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="Filtering (TRUE or FALSE)", VALUE=toString (params$filtering)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="MIND Filter (Individual with missing genotype)", VALUE=toString (params$MIND)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="GENO Filter (SNPs with missing genotype)", VALUE=toString (params$GENO)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="MAF Filter (Minor allele frequency)", VALUE=toString (params$MAF)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="HWE Filter (Hardy-Weinberg test)", VALUE=toString (params$HWE)))

	msg("Writing config file...", inputDir)
	outName = paste0(outputDir, "/out-multiGWAS-inputParameters.tbl")
	write.table (file=outName, paramsDF, quote=F, sep="\t", row.names=F)
	return (paramsDF)
}
#-------------------------------------------------------------
# Get alternate allele from a row of alleles
#-------------------------------------------------------------
#-------------------------------------------------------------
# Call to main function (first lines)
#-------------------------------------------------------------
#main ()

