#!/usr/bin/Rscript

# INFO   : Create different summarized reports (tables, plots, Venn diagrams) for multiGWAS tool analysis
# AUTHOR : Luis Garreta (lgarreta@agrosavia.co)
# DATA   : feb/2020
# LOG: 
#	r3.0: Manhattan and QQ plots. Formated to create a Markdown report (not yet)
#	r2.0: Improve selection of data from tables. Titles to graphics and files
#	r1.0: Message to log files
#	r0.9: Full working with funtions: create snpsTables and ven diagrams using parameters
#	r0.8: Create venn diagrams, summary table of first Ns"
#     

## @knitr loadLibraries
library (stringr)
library (dplyr)
library (qqman)
library (VennDiagram)
library (config)  # For read config file
KNITFLAG=FALSE

options (width=300)
#options(scipen=999)

## @knitr initParams
inputDir    = "in/"
outputDir   = "out/"
gwasType    = "Naive"
reportTitle = "GWAS Report"
nBEST = 6

## kinitr main
#-------------------------------------------------------------
# Main function
# Input files are taken from input dir
# Outupt are written to output dir
#-------------------------------------------------------------
main <- function () {
	msg ("Main...")

	createDir (outputDir)
	createReports (inputDir, gwasType, reportTitle, outputDir)
}

#-------------------------------------------------------------
# Function to create different reports:
#	1- 1 table of best SNPs
#	2- 1 table of significative SNPs
#   3- 1 Venn diagram of best SNPs
#   4- 1 Venn diagram of significative SNPs
#	5- 1 multiplot of 4x4 manhattan and QQ plots
#-------------------------------------------------------------
createReports <- function (inputDir, gwasType, title, outputDir, nBEST=7) 
{
	# Get the two summary snpsTables: best SNPs, and significative SNPs
	snpsTables = markersSummaryTable ("in/", "Structure", "Tittle", "out/",  nBEST=7)

	# Title: Report from MultiGWAS tool
	## The report contains the following elements:
	## 1. Table of N best ranged SNPs
	## 2. Table of significative SNPs
	## 3. Venn diagrams of N best ranged SNPs
	## 4. Venn diagrams of sinificative SNPs
	## 5. Manhattan and QQ plots

	# 1. Table of N best ranged SNPs
	## Table with the first N best ranged SNPs from the GWAS analysis
	msg ("Writing Best SNPS Table...")
	outName = paste0(outputDir, "/out-summary-gwas-best", nBEST)
	write.table (file=paste0(outName,".scores"), snpsTables$best, row.names=F,quote=F, sep="\t")

	# 2. Table of significative SNPs
	## Table with the significative SNPs from the GWAS analysis
	msg ("Writing Significative SNPS Table...")
	outName = paste0(outputDir, "/out-summary-gwas-signficatives.scores")
	write.table (file=outName, snpsTables$significatives, row.names=F,quote=F, sep="\t")

	# 3. Venn diagrams of N best ranged SNPs
	## Venn diagram showing the N best ranged SNPs 
	scoresType = paste0("Best",nBEST)
	outFilename = paste0(outputDir,"/out-summary-gwas-venndiagram-", scoresType, ".pdf")
	pdf (outFilename)
	commonBest = markersVennDiagrams (snpsTables$best, title, outputDir)
	dev.off()

	# 4. Venn diagrams of sinificative SNPs
	## Venn diagram showing the significative SNPs from the GWAS analysis
	scoresType = "Significatives"
	outFilename = paste0(outputDir,"/out-summary-gwas-venndiagram-", scoresType, ".pdf")
	pdf (outFilename)
	commonSign = markersVennDiagrams (snpsTables$significatives, "Significatives", title, outputDir)
	dev.off ()

	# 5. Manhattan and QQ plots
	## Manhattan and QQ plots for the GWAS analysis of the four tools

	pdf (paste0 (outputDir, "/out-summary.manhattan-qq-plots.pdf"), width=11, height=7)
	markersManhattanPlots (inputDir, gwasType, commonBest, commonSign, snpsTables$significatives, outputDir)
	dev.off()

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
markersManhattanPlots <- function (inputDir, gwasType, commonBest, commonSign, summarySignificatives, outputDir, nBest=8) {
	files =  list.files(inputDir, pattern=paste0("^(.*(",gwasType,").*(scores)[^$]*)$"), full.names=T)
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
		else if (grepl ("Plink", filename)) {
			tool = "Plink"
			gwasResults = data.frame (SNP=data$SNP, CHR=data$CHR, BP=data$POS, P=10^-data$SCORE)
		}
		else if (grepl ("Tassel", filename)) {
			tool = "Tassel"
			gwasResults = data.frame (SNP=data$Marker, CHR=data$Chr, BP=data$Pos, P=10^-data$SCORE)
		}
		ss = summarySignificatives
		if (tool %in% ss$TOOL)
			signThresholdScore = min (ss [ss$TOOL==tool,"SCORE"])
		else
			signThresholdScore = 0.95*ceiling (data[1, "SCORE"]) 

		colorsBlueOrange = c("blue4", "orange3")
		manhattan(gwasResults,col = c("gray10", "gray60"), highlight=commonBest, annotatePval=bestThreshold, annotateTop=F,
				  suggestiveline=bestThresholdScore, genomewideline=signThresholdScore, main=mainTitle, logp=T, cex=2)

		text (x=0, y=signThresholdScore*1.02, "Significative",, col="red", pos=4)
		text (x=0, y=bestThresholdScore*0.92, "Best",, col="blue", pos=4)

		datax = calculateInflationFactor (-log10 (gwasResults$P))
		qq (gwasResults$P)
		mtext (bquote(lambda[GC] == .(datax$delta)), side=3, line=-2, cex=0.7)
		#title (bquote(lambda[GC] == .(datax$delta)))
	}
	par (op)
	#dev.off()
}

#------------------------------------------------------------------------
# Create Venn diagram of common markers using info from summary table
#------------------------------------------------------------------------
## @knitr markersVennDiagrams
markersVennDiagrams <- function (summaryTable, scoresType, title="", outputDir="out"){
	flog.threshold(ERROR)
	x <- list()
	x$GWASpoly = summaryTable %>% filter (TOOL %in% "GWASpoly") %>% select (SNP) %>% .$SNP
	x$SHEsis   = summaryTable %>% filter (TOOL %in% "SHEsis") %>% select (SNP) %>% .$SNP
	x$Plink    = summaryTable %>% filter (TOOL %in% "Plink")  %>% select (SNP) %>% .$SNP
	x$Tassel   = summaryTable %>% filter (TOOL %in% "Tassel") %>% select (SNP) %>% .$SNP

	a = intersect (x$GWASpoly, x$SHEsis)
	b = intersect (x$GWASpoly, x$Plink)
	c = intersect (x$GWASpoly, x$Tassel)
	d = intersect (x$SHEsis,   x$Plink)
	e = intersect (x$SHEsis,   x$Tassel)
	f = intersect (x$Plink,    x$Tassel)
	#commonSNPs = intersect (intersect (x$GWASpoly, x$SHEsis), intersect(x$Plink, x$Tassel))
	commonSNPs = union (union (a,b),union (union (c,d),union (e,f)))

	mainTitle = paste0(title, "-", scoresType)
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
## @knitr markersSummaryTable
markersSummaryTable <- function (inputDir, gwasType, title="", outputDir="out", nBEST=5, significanceLevel=0.05, correctionMethod="FDR") {
	map = read.table (file=paste0(inputDir,"/map.tbl"))
	rownames (map) = map [,1]

	files =  list.files(inputDir, pattern=paste0("^(.*(",gwasType,").*(scores)[^$]*)$"), full.names=T)
	msg ("CREATING THE SUMMARY...")
	summaryTable = data.frame ()

	tool=""
	for (f in files) {
		msg ("Processing file: ", f)
		data <- read.table (file=f, header=T)
		if (nrow(data)>nBEST) data=data [1:nBEST,] 
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
		}else if (grepl ("Plink", f)) {
			tool    = "Plink"
			snps    = data$SNP
			chrom   = data$CHR
			pos	    = map [snps, "Position"]
			flagNewData = T
		}else if (grepl ("Tassel", f)) {
			tool    = "Tassel"
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
			dfm = data.frame (TOOL=tool, MODEL=gwasType, CHR=chrom, POS=pos, SNP=snps, P = round (pVal,6), SCORE=pscores, THRESHOLD=tscores, SIGNF=signf )
			dfm = dfm %>% distinct (SNP, .keep_all=T)
			summaryTable <- rbind (summaryTable, dfm)
			flagNewData = F
		}
	}

	summaryTable = summaryTable [which(!is.na(summaryTable$SIGNF)),]
	#outName = paste0(outputDir, "/out-summary-gwas-best", nBEST)
	#msg ("Writting summary results to ", outName, "...")
	#write.table (file=paste0(outName,".scores"), summaryTable, row.names=F,quote=F, sep="\t")
	#commonBest = markersVennDiagrams (summaryTable, paste0("Best",nBEST), title, outputDir)

	summarySignificatives = summaryTable %>% filter (SIGNF%in%T) 
	#outName = paste0(outputDir, "/out-summary-gwas-signficatives.scores")
	#write.table (file=outName, summarySignificatives, row.names=F,quote=F, sep="\t")
	#commonSign = markersVennDiagrams (summarySignificatives, "Significatives", title, outputDir)

	#markersManhattanPlots (inputDir, gwasType, commonBest, commonSign, summarySignificatives)

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
	if (KNITFLAG==FALSE) { 
		print (KNITFLAG)
		messages = unlist (list (...))
		cat (">>>>", messages, "\n")
	}
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
## @knitr getConfigurationParameters
getConfigurationParameters <- function (inputDir) 
{
	msg("Reading config file...")
	configFile = paste0 (inputDir,"/", list.files (inputDir, pattern="config")[1])
	print (configFile)

	paramsDF = data.frame (PARAMETER=character(), VALUE=character ())

	params = config::get (file=configFile) 
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="Genotype filename", VALUE=toString (params$genotypeFile)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="Phenotype filename", VALUE=toString (params$phenotypeFile)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="Significance level", VALUE=toString (params$significanceLevel)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="Correction method", VALUE=toString (params$correctionMethod)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="Trait", VALUE=toString (params$trait)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="GWAS model", VALUE=toString (params$gwasModel)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="Filtering", VALUE=toString (params$filtering)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="MIND Filter (Individual with missing genotype)", VALUE=toString (params$MIND)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="GENO Filter (SNPs with missing genotype)", VALUE=toString (params$GENO)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="MAF Filter (Minor allele frequency)", VALUE=toString (params$MAF)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="HWE Filter (Hardy-Weinberg test)", VALUE=toString (params$HWE)))

	#row.names = c("PARAMETER", "VALUE"))
	return (paramsDF)
}
#-------------------------------------------------------------
# Get alternate allele from a row of alleles
#-------------------------------------------------------------
#-------------------------------------------------------------
# Call to main function (first lines)
#-------------------------------------------------------------
main ()

