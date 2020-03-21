#!/usr/bin/Rscript
## Log: 
##      r1.0: Message to log files
##      r0.9: Full working with funtions: create tables and ven diagrams using parameters
##      r0.8: Create venn diagrams, summary table of first Ns"
##     

library (stringr)
library (dplyr)
options (width=300)
#options(scipen=999)

#------------------------------------------------------------------------
# Create Venn diagram of common markers using info from summary table
#------------------------------------------------------------------------
markersVennDiagrams <- function (summaryTable, scoresType, outDir="out"){
	require(VennDiagram)
	flog.threshold(ERROR)
	x <- list()
	x$Gwasp4 = summaryTable %>% filter (TOOL %in% "Gwasp4") %>% select (SNP) %>% .$SNP
	x$Shesis = summaryTable %>% filter (TOOL %in% "Shesis") %>% select (SNP) %>% .$SNP
	x$Plink  = summaryTable %>% filter (TOOL %in% "Plink")  %>% select (SNP) %>% .$SNP
	x$Tassel = summaryTable %>% filter (TOOL %in% "Tassel") %>% select (SNP) %>% .$SNP

	v0 <-venn.diagram(x, height=9000, width=9000, alpha = 0.5, filename = NULL,
						col = c("red", "blue", "green", "yellow"),
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
markersSummaryTable <- function (inputDir, gwasType, outDir="out", nBEST=5, SIGNIFICANCE=0.05) {
	nMARKERS = 2017
	THRESHOLD = round (-log10 (SIGNIFICANCE/nMARKERS),4)

	files =  list.files(inputDir, pattern=paste0("^(.*(",gwasType,").*(scores)[^$]*)$"), full.names=T)
	message (">>>> SUMMARY FOR THE FILES: ")
	summaryTable = data.frame ()

	tool=""
	for (f in files) {
		if (str_detect(f, "Gwasp4")) {
			#if (str_detect(f, "Gwasp4")) tool <- "Gwasp4" else tool <- "Gwasp2"
			data = read.table (file=f, header=T)
			tool    = "Gwasp4"
			data    = data %>% add_count (Marker, sort=T, name="N1") %>% arrange (desc(N1), desc(Score))
			data    = data %>% distinct (Marker, .keep_all=T)
			data    = if (ncol(data)>nBEST) data [1:nBEST,] else data
			snps    <- data$Marker
			pVal	<- round (10^(-data$Score),10)
			pscores <- data$Score
			tscores <- data$Threshold
			chrom   <- data$Chrom
			pos	    <- data$Position
			signf   = pscores >= tscores
		}else if (str_detect (f, "Plink")) {
			data = read.table (file=f, header=T)
			data    = if (ncol(data)>nBEST) data [1:nBEST,] else data
			tool    = "Plink"
			snps    = data$SNP
			pVal    = p.adjust (data$UNADJ, "fdr")
			pscores = round (-log10 (pVal), 4)
			tscores = THRESHOLD
			chrom   = data$CHR
			pos	    = NA
			signf   = pscores >= tscores
			
		}else if (str_detect (f, "Tassel")) {
			data    = read.table (file=f, header=T)
			data    = if (ncol(data)>nBEST) data [1:nBEST,] else data
			tool    = "Tassel"
			snps    = data$Marker
			#pVal    = unlist (data %>% rowwise %>% mutate (minP=min(p, add_p, dom_p, na.rm=T)) %>% select (minP))
			#pVal    = p.adjust (pVal, "fdr")
			pVal    = data$minP
			pscores = unlist (round (-log10 (pVal),4))
			tscores = THRESHOLD
			chrom   = data$Chr
			pos		= data$Pos
			signf   = pscores >= tscores
		}else if (str_detect (f, "Shesis")) {
			data    = read.table (file=f, header=T, sep="\t")
			data    = data %>% arrange (P.value) %>% head (nBEST) #top_n (-1*nBEST,P.value)
			tool    = "Shesis"
			snps    = data$SNP
			pVal    = data$P.value
			pVal    = p.adjust (pVal, "fdr")
			pscores = unlist (round (-log10 (pVal),4))
			tscores = THRESHOLD
			chrom   = NA
			pos		= NA
			signf   = pscores >= tscores
		}
		message ("    ", tool, ": ", f)
		dfm = data.frame (TOOL=tool, MODEL=gwasType, CHR=chrom, POS=pos, SNP=snps, P = round (pVal,6), FDRSCORE=pscores, THRESHOLD=tscores, SIGNF=signf )
		dfm = dfm %>% distinct (SNP, .keep_all=T)
		summaryTable <- rbind (summaryTable, dfm)
	}
	summaryTable = summaryTable [which(!is.na(summaryTable$SIGNF)),]
	outName = paste0(outDir, "/out-summary-gwas-best", nBEST)
	msg ("Writing summary results to ", outName, "...")
	write.table (file=paste0(outName,".scores"), summaryTable, row.names=F,quote=F, sep="\t")
	markersVennDiagrams (summaryTable, paste0("best",nBEST), outDir)
	summarySignificatives = summaryTable %>% filter (SIGNF%in%T) 
	outName = paste0(outDir, "/out-summary-gwas-signficatives.scores")
	write.table (file=outName, summarySignificatives, row.names=F,quote=F, sep="\t")
	markersVennDiagrams (summarySignificatives, "significatives", outDir)

	return (summaryTable)
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
#markersSummaryTable ("out/", "Naive", "out/")

