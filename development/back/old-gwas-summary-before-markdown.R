#!/usr/bin/Rscript
## Log: 
##      r2.0: Improve selection of data from tables. Titles to graphics and files
##      r1.0: Message to log files
##      r0.9: Full working with funtions: create tables and ven diagrams using parameters
##      r0.8: Create venn diagrams, summary table of first Ns"
##     

library (stringr)
library (dplyr)
library(qqman)

options (width=300)
#options(scipen=999)

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

#------------------------------------------------------------------------
#------------------------------------------------------------------------
markersManhattanPlots <- function (inputDir, gwasType, commonBest, commonSign, summarySignificatives, nBest=8) {
	files =  list.files(inputDir, pattern=paste0("^(.*(",gwasType,").*(scores)[^$]*)$"), full.names=T)
	pdf ("out-qqman.pdf", width=11, height=7)
	op <- par(mfcol = c(2,4), mar=c(3.5,3.5,3,1), oma=c(0,0,0,0), mgp = c(2.2,1,0))
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
				  suggestiveline=bestThresholdScore, genomewideline=signThresholdScore, main=mainTitle, logp=T)

		msg ("THRESHOLD: ", bestThresholdScore, "thScr: ",  signThresholdScore, "file: ", filename)

		text (x=0, y=signThresholdScore*1.02, "Significative",, col="red", pos=4)
		text (x=0, y=bestThresholdScore*0.92, "Best",, col="blue", pos=4)

		datax = calculateInflationFactor (-log10 (gwasResults$P))
		qq (gwasResults$P)
		mtext (bquote(lambda[GC] == .(datax$delta)), side=3, line=-2, cex=0.7)
		#title (bquote(lambda[GC] == .(datax$delta)))
	}
	par (op)
	dev.off()
}

#------------------------------------------------------------------------
# Create Venn diagram of common markers using info from summary table
#------------------------------------------------------------------------
markersVennDiagrams <- function (summaryTable, scoresType, title="", outDir="out"){
	require(VennDiagram)
	st <<-summaryTable
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
	v0 <- venn.diagram(x, height=12000, width=12000, alpha = 0.5, filename = NULL, main=mainTitle,
						col = c("red", "blue", "green", "yellow"), cex=0.9,
						fill = c("red", "blue", "green", "yellow")) 

	overlaps <- calculate.overlap(x)
	overlaps <- rev(overlaps)

	posOverlap = as.numeric (gsub ("a","", (names (overlaps))))
	for (i in 1:length(overlaps)){
		pos = posOverlap [i]
		v0[[pos+8]]$label <- paste(overlaps[[i]], collapse = "\n")
 	}
 
	pdf(paste0(outDir,"/out-summary-gwas-venndiagram-", scoresType, ".pdf"))
	grid.draw(v0)
	dev.off()

	return (commonSNPs)
}

#------------------------------------------------------------------------
# Create a summary table of best and significative markers
#------------------------------------------------------------------------
markersSummaryTable <- function (inputDir, gwasType, title="", outDir="out", nBEST=5, significanceLevel=0.05, correctionMethod="FDR") {
	print (getwd())
	map = read.table (file=paste0(outDir,"/map.tbl"))
	rownames (map) = map [,1]

	files =  list.files(inputDir, pattern=paste0("^(.*(",gwasType,").*(scores)[^$]*)$"), full.names=T)
	msg ("CREATING THE SUMMARY...")
	summaryTable = data.frame ()

	tool=""
	msg ("Max number of best scored SNPs:", nBEST)
	msg ("Tools:")
	print (files)
	for (f in files) {
		msg ("Tool file: ", f)
		data <- read.table (file=f, header=T)
		if (nrow(data)>nBEST) data=data [1:nBEST,] 
		pVal	<- data$P
		pscores <- data$SCORE
		tscores <- data$THRESHOLD
		signf   = pscores >= tscores

		flagNewData = F
		if (str_detect(f, "GWASpoly")) {
			tool    = "GWASpoly"
			snps    <- data$Marker
			chrom   <- data$Chrom
			pos	    <- data$Position
			flagNewData = T
		}else if (str_detect (f, "Plink")) {
			tool    = "Plink"
			snps    = data$SNP
			chrom   = data$CHR
			pos	    = map [snps, "Position"]
			flagNewData = T
		}else if (str_detect (f, "Tassel")) {
			tool    = "Tassel"
			snps    = data$Marker
			chrom   = data$Chr
			pos		= data$Pos
			flagNewData = T
		}else if (str_detect (f, "SHEsis")) {
			tool    = "SHEsis"
			snps    = data$SNP
			chrom	= map [snps, "Chrom"]
			pos	    = map [snps, "Position"]
			flagNewData = T
		}
		if (flagNewData==T) {
			msg ("    ", tool)
			dfm = data.frame (TOOL=tool, MODEL=gwasType, CHR=chrom, POS=pos, SNP=snps, P = round (pVal,6), SCORE=pscores, THRESHOLD=tscores, SIGNF=signf )
			dfm = dfm %>% distinct (SNP, .keep_all=T)
			summaryTable <- rbind (summaryTable, dfm)
			flagNewData = F
		}
	}

	summaryTable = summaryTable [which(!is.na(summaryTable$SIGNF)),]
	outName = paste0(outDir, "/out-summary-gwas-best", nBEST)
	msg ("Writting summary results to ", outName, "...")
	write.table (file=paste0(outName,".scores"), summaryTable, row.names=F,quote=F, sep="\t")
	commonBest = markersVennDiagrams (summaryTable, paste0("Best",nBEST), title, outDir)

	summarySignificatives = summaryTable %>% filter (SIGNF%in%T) 
	outName = paste0(outDir, "/out-summary-gwas-signficatives.scores")
	write.table (file=outName, summarySignificatives, row.names=F,quote=F, sep="\t")
	commonSign = markersVennDiagrams (summarySignificatives, "Significatives", title, outDir)

	markersManhattanPlots (inputDir, gwasType, commonBest, commonSign, summarySignificatives)

	return (summaryTable)
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
#-------------------------------------------------------------
# Main
#-------------------------------------------------------------

# Create Venn diagram of common markers
#markersSummaryTable ("./", "Structure", "Tittle", "out/",  nBEST=7)

