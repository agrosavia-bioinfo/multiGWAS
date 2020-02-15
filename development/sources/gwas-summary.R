#!/usr/bin/Rscript
# INFO   : Functios to summarize results from the different tools
# AUTHOR : Luis Garreta (lgarreta@agrosavia.co)
# DATE   : Feb/2020
# LOG: 
#      r2.0: Improve selection of data from tables. Titles to graphics and files
#      r1.0: Message to log files
#      r0.9: Full working with funtions: create tables and ven diagrams using parameters
#      r0.8: Create venn diagrams, summary table of first Ns"

library (stringr)
library (dplyr)
options (width=300)
#options(scipen=999)

#------------------------------------------------------------------------
# Create Venn diagram of common markers using info from summary table
#------------------------------------------------------------------------
markersVennDiagrams <- function (summaryTable, scoresType, title="", outDir="out"){
	require(VennDiagram)
	st <<-summaryTable
	flog.threshold(ERROR)
	x <- list()
	x$Gwaspoly = summaryTable %>% filter (TOOL %in% "Gwaspoly") %>% select (SNP) %>% .$SNP
	x$Shesis = summaryTable %>% filter (TOOL %in% "Shesis") %>% select (SNP) %>% .$SNP
	x$Plink  = summaryTable %>% filter (TOOL %in% "Plink")  %>% select (SNP) %>% .$SNP
	x$Tassel = summaryTable %>% filter (TOOL %in% "Tassel") %>% select (SNP) %>% .$SNP

	mainTitle = paste0(title, "-", scoresType)
	msg();msg (mainTitle);msg()
	v0 <-venn.diagram(x, height=12000, width=12000, alpha = 0.5, filename = NULL, main=mainTitle,
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
}

#------------------------------------------------------------------------
# Create a summary table of best and significative markers
#------------------------------------------------------------------------
markersSummaryTable <- function (inputDir, gwasType, title="", outDir="out", nBEST=5, significanceLevel=0.05, correctionMethod="FDR") {
	map = read.table (file="out/map.tbl")
	rownames (map) = map [,1]

	files =  list.files(inputDir, pattern=paste0("^(.*(",gwasType,").*(scores)[^$]*)$"), full.names=T)
	msg ("CREATING THE SUMMARY...")
	summaryTable = data.frame ()

	tool=""
	msg ("Max number of best scored SNPs:", nBEST)
	msg ("Tools:")
	for (f in files) {
		data <- read.table (file=f, header=T)
		if (nrow(data)>nBEST) data=data [1:nBEST,] 
		pVal	<- data$P
		pscores <- data$SCORE
		tscores <- data$THRESHOLD
		signf   = pscores >= tscores

		flagNewData = F
		if (str_detect(f, "Gwaspoly")) {
			tool    = "Gwaspoly"
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
		}else if (str_detect (f, "Shesis")) {
			tool    = "Shesis"
			snps    = data$SNP
			chrom	= map [snps, "Chrom"]
			pos	    = map [snps, "Position"]
			flagNewData = T
		}
		if (flagNewData==T) {
			msg ("    ", tool)
			dfm = data.frame (TOOL=tool, MODEL=gwasType, CHR=chrom, POS=pos, SNP=snps, P = round (pVal,6), FDRSCORE=pscores, THRESHOLD=tscores, SIGNF=signf )
			dfm = dfm %>% distinct (SNP, .keep_all=T)
			summaryTable <- rbind (summaryTable, dfm)
			flagNewData = F
		}
	}

	summaryTable = summaryTable [which(!is.na(summaryTable$SIGNF)),]
	outName = paste0(outDir, "/out-summary-gwas-best", nBEST)
	msg ("Writting summary results to ", outName, "...")
	write.table (file=paste0(outName,".scores"), summaryTable, row.names=F,quote=F, sep="\t")
	markersVennDiagrams (summaryTable, paste0("Best",nBEST), title, outDir)
	summarySignificatives = summaryTable %>% filter (SIGNF%in%T) 
	outName = paste0(outDir, "/out-summary-gwas-signficatives.scores")
	write.table (file=outName, summarySignificatives, row.names=F,quote=F, sep="\t")
	markersVennDiagrams (summarySignificatives, "Significatives", title, outDir)

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

# Create Venn diagram of common markers
#markersSummaryTable ("out/", "Naive", "Tittle", "out/",  nBEST=7)

