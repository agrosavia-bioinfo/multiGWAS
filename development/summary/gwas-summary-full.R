#!/usr/bin/Rscript
## Log: 
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
markersVennDiagrams <- function (summaryTable, typeScores){
	require(VennDiagram)
	flog.threshold(ERROR)
	x <- list()
	x$Gwasp4 = summaryTable %>% filter (TOOL %in% "Gwasp4") %>% select (SNP) %>% .$SNP
	x$Gwasp2 = summaryTable %>% filter (TOOL %in% "Gwasp2") %>% select (SNP) %>% .$SNP
	x$Plink  = summaryTable %>% filter (TOOL %in% "Plink")  %>% select (SNP) %>% .$SNP
	x$Tassel = summaryTable %>% filter (TOOL %in% "Tassel") %>% select (SNP) %>% .$SNP

	v0 <<-venn.diagram(x, height=9000, width=9000,
										col = c("red", "blue", "green", "yellow"),
										fill = c("red", "blue", "green", "yellow"), 
										alpha = 0.5, filename = NULL)

	overlaps <- calculate.overlap(x)
	overlaps <- rev(overlaps)

	posOverlap = as.numeric (gsub ("a","", (names (overlaps))))
	for (i in 1:length(overlaps)){
		pos = posOverlap [i]
		v0[[pos+8]]$label <- paste(overlaps[[i]], collapse = "\n")
	}

	pdf(paste0("out-venn-diagram-", typeScores, ".pdf"))
	grid.draw(v0)
	dev.off()
}
#------------------------------------------------------------------------
# Create a summary table of best and significative markers
#------------------------------------------------------------------------
markersSummaryTable <- function (model, nBEST=10, SIGNIFICANCE=0.05) {
	nMARKERS = 2017
	THRESHOLD = round (-log10 (SIGNIFICANCE/nMARKERS),4)

	model<<-"Naive"
	files =  list.files(".", pattern="^(.*(Naive).*(tbl)[^$]*)$")
	files = c("out-Gwasp4-Naive-significativeQTLs.tbl", "out-Gwasp2-Naive-significativeQTLs.tbl", 
			  "out-Plink-Naive-assoc.linear.adjusted.tbl", "out-Tassel-Naive-GLM_Stats_geno+pheno.tbl")
	#files = c("out-Gwasp4-Naive-significativeQTLs.tbl", 
	summaryTable = data.frame ()

	for (f in files) {
		data = read.table (file=f, header=T)
		print (f)
		if (str_detect(f, "Gwasp")) {
			if (str_detect(f, "Gwasp4")) tool <- "Gwasp4" else tool <- "Gwasp2"
			data    = data %>% add_count (Marker, sort=T, name="N1") %>% arrange (desc(N1), desc(Score))
			data    = data %>% distinct (Marker, .keep_all=T)
			#write.table (file=paste0("data.scores-",tool,".scores"), data, quote=F, sep="\t")
			data    = data [1:nBEST,]
			snps    <- data$Marker
			pVal	<- round (10^(-data$Score),10)
			pscores <- data$Score
			tscores <- data$Threshold
			chrom   <- data$Chrom
			pos	    <- data$Position
		}else if (str_detect (f, "Plink")) {
			data    = data [1:nBEST,]
			tool    = "Plink"
			snps    = data$SNP
			pVal    = p.adjust (data$UNADJ, "fdr")
			pscores = round (-log10 (pVal), 4)
			tscores = THRESHOLD
			chrom   = data$CHR
			pos	    = NA
			
		}else if (str_detect (f, "Tassel")) {
			data    = data %>% top_n (-1*nBEST,p)
			tool    = "Tassel"
			snps    = data$Marker
			pVal    = unlist (data %>% rowwise %>% mutate (minP=min(p, add_p, dom_p, na.rm=T)) %>% select (minP))
			pVal    = p.adjust (pVal, "fdr")
			pscores = unlist (round (-log10 (pVal),4))
			tscores = THRESHOLD
			chrom   = data$Chr
			pos		= data$Pos
		}
		signf   = pscores >= tscores
		#data    = data %>% top_n (-1*nBEST,pscores)
		print (tool)
		print (length (pVal))
		dfm <- data.frame (TOOL=tool, MODEL=model, CHR=chrom, POS=pos, SNP=snps, P = round (pVal,6), FDRSCORE=pscores, THRESHOLD=tscores, SIGNF=signf )
		dfm = dfm %>% distinct (SNP, .keep_all=T)
		summaryTable <- rbind (summaryTable, dfm)
	}
	summaryTable = summaryTable [which(!is.na(summaryTable$SIGNF)),]
	message ("Writing...")
	outName = paste0("out-summary-gwas-best", nBEST)
	write.table (file=paste0(outName,".scores"), summaryTable, row.names=F,quote=F, sep="\t")
	markersVennDiagrams (summaryTable, paste0("best",nBEST))
	summarySignificatives = summaryTable %>% filter (SIGNF%in%T) 
	write.table (file="out-summary-gwas-significatives.scores", summarySignificatives, row.names=F,quote=F, sep="\t")
	markersVennDiagrams (summarySignificatives, "significatives")

	return (summaryTable)
}

# Create Venn diagram of common markers
summaryTable <<- markersSummaryTable ("Naive")
