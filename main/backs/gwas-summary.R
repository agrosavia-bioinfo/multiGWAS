#!/usr/bin/Rscript

# INFO   : Create different summarized reports (tables, plots, Venn diagrams) for multiGWAS tool analysis
# AUTHOR : Luis Garreta (lgarreta@agrosavia.co)
# DATA   : feb/2020
# LOG: 
#	r5.2: Fixed unused factors in gwasResults
#	r5.1: Fixed chord diagram when no shared SNPs
#	r5.0: Heuristic to select best gene action model
#	r4.1: Fixed transparency and axis error (options (bitmapType="cairo"))

HOME = Sys.getenv ("MULTIGWAS_HOME")
.libPaths (paste0(HOME, "/opt/Rlibs"))

suppressMessages (library (dplyr))
suppressMessages (library (qqman))
suppressMessages (library (VennDiagram))
suppressMessages (library (yaml))          # For read config file
suppressMessages (library ("RColorBrewer"))  # For chord diagrams
suppressMessages (library(circlize))         # For chord diagrams
suppressMessages (library (doParallel))
suppressMessages (library(ldsep))         # For linkage disequilibrium in matchSNPs:w

options (bitmapType="cairo", width=300)
#options(scipen=999)

#-------------------------------------------------------------
# Main function
# Input files are taken from input dir
# Outupt are written to output dir
#-------------------------------------------------------------
main <- function () {
	options (warn=1)
	source (paste0(HOME,"/sources/gwas-heatmap.R"))            # Module with functions to create summaries: tables and venn diagrams
	source (paste0(HOME,"/sources/gwas-preprocessing.R"))      # Module with functions to convert between different genotype formats 

	msgmsg ("Main...")

	inputDir        = "out/"
	genotypeFile    = "out/filtered-gwasp4-genotype.tbl"
	phenotypeFile   = "out/filtered-gwasp4-phenotype.tbl"
	genotypeNumFile = "out/filtered-gwasp4-genotype-NUM.csv"
	outputDir       = "report/"
	gwasModel       = "naive"
	nBest           = 10
	ploidy          = 4
	geneAction      = "additive"

	tools1 = list(list (tool="GWASpoly", scoresFile="out/tool-GWASpoly-scores-full.csv")) 

	tools2 = list(list (tool="GWASpoly", scoresFile="out/tool-GWASpoly-scores-full.csv"), 
				  list (tool="SHEsis", scoresFile="out/tool-SHEsis-scores-full.csv"))

	tools3 = list(list (tool="GWASpoly", scoresFile="out/tool-GWASpoly-scores-full.csv"), 
				  list (tool="SHEsis", scoresFile="out/tool-SHEsis-scores-full.csv"),
				  list (tool="TASSEL", scoresFile="out/tool-TASSEL-scores-full.csv"))

	tools4 = list(list (tool="GWASpoly", scoresFile="out/tool-GWASpoly-scores-full.csv"), 
				  list (tool="SHEsis", scoresFile="out/tool-SHEsis-scores-full.csv"),
				  list (tool="TASSEL", scoresFile="out/tool-TASSEL-scores-full.csv"),
				  list (tool="PLINK", scoresFile="out/tool-PLINK-scores-additive-full.csv"))

	listOfResultsFile = tools1

	createReports (inputDir, genotypeFile, phenotypeFile, genotypeNumFile,
				   ploidy, gwasModel, outputDir, nBest, geneAction, listOfResultsFile)
}

#-------------------------------------------------------------
# Function to create different reports:
#	1- 1 table of best SNPs
#	2- 1 table of significative SNPs
#   3- 1 Venn diagram of best SNPs
#   4- 1 Venn diagram of significative SNPs
#	5- 1 multiplot of 4x4 manhattan and QQ plots
#-------------------------------------------------------------
#-------------------------------------------------------------
createReports <- function (listOfResultsFile, params) {
	inputDir  = params$outputDir
	outputDir = params$reportDir

	msgmsg ("Creating tables and plots with summary results for ", params$gwasModel, " model...")
	createDir (outputDir)

	# Define filenams for outputs
	fileBestScores               = paste0 (outputDir,  "/out-multiGWAS-scoresTable-best.scores")
	fileSignificativeScores      = paste0 (outputDir,  "/out-multiGWAS-scoresTable-significatives.scores")
	fileBestVennDiagram          = paste0 (outputDir, "/out-multiGWAS-vennDiagram-best")
	fileSignificativeVennDiagram = paste0 (outputDir, "/out-multiGWAS-vennDiagram-significatives")
	fileManhattanPlotPNG         = paste0 (outputDir, "/out-multiGWAS-manhattanQQ-plots.png")
	fileManhattanPlotPDF         = paste0 (outputDir, "/out-multiGWAS-manhattanQQ-plots.pdf")

	msgmsg ("Writing input config parameters...")
	config = writeConfigurationParameters (inputDir, outputDir)

	if (length (listOfResultsFile)==0) {
		msgError ("WARNING: No result files for any tool. Check config file parameters (e.g. tools, geneAction, gwasModel)")
		quit ()
	}

	snpTables = markersSummaryTable (listOfResultsFile, params$gwasModel, params$nBest, params$geneAction, params$genotypeNumFile)

	msgmsg ("Writing tables with best ranked and signficative SNPs... ")
	write.table (file=fileBestScores, snpTables$best, row.names=F,quote=F, sep="\t")
	write.table (file=fileSignificativeScores, snpTables$significatives, row.names=F,quote=F, sep="\t")

	msgmsg ("Writing Venn diagram for best and significative SNPs...")
	commonBest = markersVennDiagrams (listOfResultsFile, snpTables$best, params$gwasModel, "Best", fileBestVennDiagram)
	commonSign = markersVennDiagrams (listOfResultsFile, snpTables$significatives, params$gwasModel, "Significatives", fileSignificativeVennDiagram)

	msgmsg ("Writing Manhattan and QQ plots...")
	png (fileManhattanPlotPNG, width=11, height=15, units="in", res=90)
	op=markersManhattanPlots (listOfResultsFile, commonBest, snpTables, params$nBest, params$geneAction)
	dev.off()

	pdf (fileManhattanPlotPDF, width=11, height=15)
	op=markersManhattanPlots (listOfResultsFile, commonBest, snpTables, params$nBest, params$geneAction)
	par (op)
	dev.off()

	# Create heat maps
	msgmsg ("Creating heatmaps for best ranked SNPs...")
	createHeatmapForSNPList (outputDir, params$genotypeFile, params$genotypeNumFile, params$phenotypeFile, commonBest, params$ploidy)

	# Create chord diagrams
	msgmsg ("Creating chord diagrams for chromosome vs SNPs...")
	createChordDiagramSharedSNPs (fileBestScores)

	# Call to rmarkdown report
	createMarkdownReport (config)
}

#-------------------------------------------------------------
# Create Rmkardown report
#-------------------------------------------------------------
createMarkdownReport  <- function (config) {
	msg ("Creating html rmarkdown report...")
	outputFile = paste0 (config$workingDir, "/multiGWAS-report.html") 
	title      = paste0 ("MultiGWAS report for ", config$gwasModel, " GWAS model")


	# Create html with embbeded images (using javascrips) for Java WebView
	rmarkdown::render (paste0(HOME,"/sources/gwas-markdown.Rmd"), output_file=outputFile, output_format="html_document", 
					   params=list (workingDir=config$workingDir, reportTitle=title, nBest=config$nBest), quiet=T)

	# Create html with external images (using subfolders) for Java JEditorPane
	#rmarkdown::render (paste0(HOME,"/sources/gwas-markdown.Rmd"), output_file=outputFile, 
	#				   output_format="html_document", output_options=list(self_contained=F),
	#				   params=list (workingDir=config$workingDir, reportTitle=title, nBest=config$nBest), quiet=T)
}


#------------------------------------------------------------------------
#------------------------------------------------------------------------
markersManhattanPlots <- function (listOfResultsFile, commonBest, snpTables, nBest, geneAction) {
	op <- par(mfrow = c(4,2), mar=c(3.5,3.5,3,1), oma=c(0,0,0,0), mgp = c(2.2,1,0))
	#op <- layout (matrix (c(1,1,2,3, 4:16), 4,4, byrow=T))#, widths=c(2,1), heights=c(2,1)))
	for (res in listOfResultsFile) {
		tool       = res$tool
		scoresFile = res$scoresFile

		data = read.table (file=scoresFile, header=T)
		data = data [!is.na (data$P),]
		#data = selectBestModel (data, nBest, tool, geneAction)
		gwasResults = data.frame (SNP=data$Marker, CHR=data$CHR, BP=data$POS, P=10^-data$SCORE, MODEL=data$MODEL)

		# Remove old column factors
		gwasResults[] <- lapply(gwasResults, function(x) if(is.factor(x)) factor(x) else x)

		# Check for few rows
		nRows = nrow (data)
		nBest = ifelse (nBest>nRows, nRows, nBest )

		# Select limiting lines for Manhattan (best and significant)
		bestThresholdScore = data [nBest,c("SCORE")]
		bestThreshold      = 10^-bestThresholdScore
		signThresholdScore = data [1, "THRESHOLD"]
		signThreshold      = 10^-signThresholdScore

		# 
		ss = snpTables$significatives
		if (tool %in% ss$TOOL)
			signThresholdScore = min (ss [ss$TOOL==tool,"SCORE"])
		else
			signThresholdScore = ceiling (data[1, "SCORE"]) 

		bestSNPsTool     = unlist (dplyr::select (filter (snpTables$best, TOOL==tool), "SNP"))
		sharedSNPs       = intersect (commonBest, bestSNPsTool)
		colorsBlueOrange = c("blue4", "orange3")
		ylims   = c (0, ceiling (signThresholdScore))

		# Check if all chromosome names are numeric
		anyNonNumericChrom <- function (chrs) {
			suppressWarnings (any (is.na (as.numeric (chrs))))
		}

		# if non-numeric Chromosome names, convert to numeric using factors
		chrs = as.character (gwasResults$CHR)
		if  (anyNonNumericChrom (chrs)==TRUE) {
			msgmsg ("!!!Mapping chromosome names to numbers (see 'out-mapped-chromosome-names.csv') file...")
			chrs            = as.factor (chrs)
			levels (chrs)   = 1:length (levels (chrs))
			write.csv (data.frame (ORIGINAL_CHROM=gwasResults$CHR, NEW_CHROM=chrs), "out-mapped-chromosome-names.csv", quote=F, row.names=F)
			gwasResults$CHR = as.numeric (chrs)
		}

		names      = unlist (strsplit (basename (scoresFile), "[-|.]"))
		mainTitle  = paste0 (names[2],"-", names [3])

		manhattan(gwasResults,col = c("orange", "midnightblue"), highlight=sharedSNPs, annotatePval=bestThreshold, annotateTop=F,
				  suggestiveline=bestThresholdScore, genomewideline=signThresholdScore, main=mainTitle, logp=T, cex=2)

		text (x=0, y=signThresholdScore*0.92, "               Significants",, col="red", pos=4)
		text (x=0, y=bestThresholdScore*0.92, "Best",, col="blue", pos=4)

		#datax = calculateInflationFactor (-log10 (gwasResults$P))
		qqMGWAS (gwasResults, geneAction)
		#qq (gwasResults$P)
		#mtext (bquote(lambda[GC] == .(datax$delta)), side=3, line=-2, cex=0.7)
		#title (bquote(lambda[GC] == .(data$delta)))
	}
	#par (op)
	#dev.off()
	return (op)
}

#-------------------------------------------------------------
# QQ plot
#-------------------------------------------------------------
qqMGWAS <- function(gwasResults, geneAction) {
	models = levels (gwasResults$MODEL)
	#---- local fun -----
	qqValues <- function (pValues) {
		scores = -log10 (pValues)
		datax  = calculateInflationFactor (scores)
		n      = length(datax$scores)
		unif.p = -log10(ppoints(n))
		return (list (scores=scores, unif.p=unif.p, gc=datax$delta))
	}
	#---- local fun -----
	values = qqValues (gwasResults$P)
	xMax = max (values$unif.p)
	yMax = max (values$scores)

	#par(pty="s")
	#plot.new ()
	plot(1, type="n", xlim = c(0, xMax), ylim = c(0, yMax),
		 xlab=expression(paste("Expected -log"[10],"(p)",sep="")),
		 ylab=expression(paste("Observed -log"[10],"(p)",sep="")),
		 main=paste (models,collapse=" / "))
	lines(c(0,max(values$unif.p)),c(0,max(values$unif.p)),lty=2, col=c("red"))

	colors = c("black", "magenta3", "cyan")
	LEGEND = list()
	for (i in 1:length(models)) {
		pValues  = gwasResults [gwasResults$MODEL==models[i], c("P")]
		values=qqValues (pValues)
		points (values$unif.p, values$scores, pch=16, col=colors[i])
		GC = as.expression (bquote (.(models[i]) ~ lambda[GC] == .(values$gc)))
		LEGEND = append(LEGEND, GC) 
	}

	legend ("bottomright", lty=1, legend=LEGEND, col=colors, pch=16) 
}

#------------------------------------------------------------------------
# Create Venn diagram of common markers using info from summary table
#------------------------------------------------------------------------
markersVennDiagrams <- function (listOfResultsFile, summaryTable, gwasModel, scoresType, outFile){
	WIDTH  = 6; HEIGHT = 7

	if (nrow (summaryTable) == 0) {
		v0 = grid.text ("No Venn Diagram  (without significative SNPs)")
		sharedSNPs = NULL
	}else { 
		# Params for figure shape and fonts
		CEXLABELS = 0.6; CEXTITLES = 1.0

		flog.threshold(ERROR)

		x = list()
		toolNames = c()
		for (resultsFile in listOfResultsFile) {
			toolNames  = c (toolNames, resultsFile$tool)
			markers = dplyr::select (filter (summaryTable, TOOL %in% resultsFile$tool), SNP) %>% .$SNP %>% as.character
			x          = append (x, list (markers))
		}
		names (x) = toolNames
		nTools = length (x)

		# Create Venn diagram
		mainTitle = paste0(gwasModel, "-", scoresType)
		COLORS= c("red", "blue", "yellow", "green")[1:nTools]
		v0 <- venn.diagram(x, height=1000, width=2000, alpha = 0.5, filename=NULL, col=COLORS, 
						   cex=CEXLABELS, cat.cex=CEXTITLES, margin=0.001, fill=COLORS, euler.d=F, scaled=F)
		overlaps   = rev (calculate.overlap(x))
		# Check if only one tool
		if (length (overlaps)==1) posOverlap = NA
		else                      posOverlap = as.numeric (gsub ("a","", (names (overlaps))))

		# Set items for areas, for two and one area are special cases
		if (nTools==2) {
			v0[[5]]$label <- paste0 ("\n", paste(setdiff(x[[1]], x[[2]]), collapse="\n"))  
			v0[[6]]$label <- paste0 ("\n", paste(setdiff(x[[2]], x[[1]])  , collapse="\n"))  
			v0[[7]]$label <- paste0 ("\n", paste(intersect(x[[1]], x[[2]]), collapse="\n"))  
		}else for (i in 1:length(overlaps)){
			pos = if (length (posOverlap)==1) 1 else  posOverlap [i] 
			v0[[pos+2*nTools]]$label <- paste0("\n", paste(overlaps[[i]], collapse = "\n"))
		}

		# Get shared SNPs
		message (" Getting shared SNPs...")
		dataSNPsNs     = data.frame (add_count (summaryTable, SNP, sort=T)); 
		dataSNPsShared = dataSNPsNs[dataSNPsNs$n > 1,]
		dataSNPsNoDups = dataSNPsShared [!duplicated (dataSNPsShared$SNP),]
		sharedSNPs     = dataSNPsNoDups$SNP
		message (" >>> Getting shared SNPs...")
	}

	png (paste0 (outFile,".png"), width=WIDTH, height=HEIGHT, units="in", res=120)
	grid.draw(v0); dev.off()

	pdf (paste0 (outFile,".pdf"), width=WIDTH,height=HEIGHT)
	grid.draw(v0); dev.off()
	
	return (sharedSNPs)
}

#------------------------------------------------------------------------
# Create a summary table of best and significative markers
#------------------------------------------------------------------------
markersSummaryTable <- function (listOfResultsFile, gwasModel, nBest, geneAction, genotypeFile) {
	summaryTable = data.frame (stringsAsFactors=F)

	for (resultsFile in listOfResultsFile) {
		TOOL       = resultsFile$tool
		scoresFile = resultsFile$scoresFile
		data       = read.table (file=scoresFile, header=T, stringsAsFactors=F)
		#data       = selectBestModel (data, nBest, TOOL, geneAction)

		MODEL        = data$MODEL
		GC           = data$GC
		SNP          = data$Marker
		CHROM        = data$CHR
		POSITION	 = data$POS
		PVALUE	     = round (data$P, 6)
		SCORE        = round (data$SCORE, 4)
		THRESHOLD    = round (data$THRESHOLD, 4)
		SIGNIFICANCE = SCORE >= THRESHOLD

		dfm = data.frame (TOOL, MODEL, GC, SNP, CHROM, POSITION, PVALUE, SCORE, THRESHOLD, SIGNIFICANCE)
		dfm = dfm [!duplicated (dfm$SNP),]
		if (nrow(dfm)>nBest) 
			dfm=dfm [1:nBest,] 

		summaryTable <- rbind (summaryTable, dfm)
	}
	#summaryTable = matchSNPsByLDAllTools (genotypeFile, summaryTable, 0.99)

	summaryTable = summaryTable [which(!is.na(summaryTable$SIGNIFICANCE)),]
	summarySignificatives = summaryTable %>% filter (SIGNIFICANCE%in%T) 
	return (list (best=summaryTable, significatives=summarySignificatives))
}


#------------------------------------------------------------------------
# Match SNPs by LD and select the best tagSNP when r2=1
#------------------------------------------------------------------------
matchSNPsByLDAllTools <- function (genotypeFile, scores, maxLD) {
	msgmsg ("Matching LD SNPs for ", genotypeFile, "...")
	#scores = read.table (scoresFile, sep="\t", header=T)
	geno   = read.csv (genotypeFile, row.names=1)

	# Create hash list of SNPs
	snpList = list()
	snps    = as.character (scores [, "SNP"])
	snpList = sapply  (snps, function (x) append (snpList,x))

	# Create genotype matrix from SNPs
	genoSNPs    = as.matrix (geno [snps,-1:-2])
	ldmat       = mldest(geno = genoSNPs, K = 4, nc = 7, type = "comp", se=F);
	ldSNPs      = ldmat [, c(3,4,7)]
	ldSNPs$snpi = sapply (ldSNPs$snpi, function (x) strsplit (x, "[.]")[[1]][1])
	ldSNPs$snpj = sapply (ldSNPs$snpj, function (x) strsplit (x, "[.]")[[1]][1])

	# Filter SNPs by R2
	ldSNPsFiltered = ldSNPs [ldSNPs$snpi!=ldSNPs$snpj,]
	ldSNPsFiltered = ldSNPsFiltered [!duplicated (ldSNPsFiltered[c(1,2)]),]
	ldSNPsR2       = ldSNPsFiltered [ldSNPsFiltered$r2 > maxLD,] ; ldSNPsR2

	# Match SNPs
	n = nrow (ldSNPsR2)
	for (i in n:1) {
		snpi = ldSNPsR2 [i, "snpi"]
		snpj = ldSNPsR2 [i, "snpj"]
		snpList [snpList %in% snpj] = snpi
		message (">>> snpi: ", snpi, " snpj: ", snpj)
	}
	#scores = data.frame (SNPLD=as.character (snpList), scores)
	scores$SNP = as.character (snpList)
	#outFile = addLabel (scoresFile, "LD")
	#write.table (scores, outFile, sep="\t", col.names=T, row.names=F, quote=F)
	return (scores)
}


#-------------------------------------------------------------
# Get params from config file and define models according to ploidy
#-------------------------------------------------------------
writeConfigurationParameters <- function (inputDir, outputDir) {
	paramsFile = paste0(inputDir, list.files (inputDir, pattern="config")[1])
	#params     = params::get (file=paramsFile) 
	params     = yaml.load_file (paramsFile, merge.precedence="order") 

	configDF = data.frame (PARAMETER=character(), VALUE=character ())
	configDF = rbind  (configDF, data.frame (PARAMETER="Ploidy (4 or 2)", VALUE=toString (params$ploidy)))
	configDF = rbind  (configDF, data.frame (PARAMETER="Genotype filename", VALUE=toString (params$genotypeFile)))
	configDF = rbind  (configDF, data.frame (PARAMETER="Phenotype filename", VALUE=toString (params$phenotypeFile)))
	configDF = rbind  (configDF, data.frame (PARAMETER="Significance level (Genome-wide significance level)", VALUE=toString (params$significanceLevel)))
	configDF = rbind  (configDF, data.frame (PARAMETER="Correction method (Bonferroni or FDR)", VALUE=toString (params$correctionMethod)))
	configDF = rbind  (configDF, data.frame (PARAMETER="GWAS model (Full or Naive)", VALUE=toString (params$gwasModel)))
	configDF = rbind  (configDF, data.frame (PARAMETER="nBest (Number of best-ranked SNPs to be reported)", VALUE=toString (params$nBest)))
	configDF = rbind  (configDF, data.frame (PARAMETER="Filtering (TRUE or FALSE)", VALUE=toString (params$filtering)))
	configDF = rbind  (configDF, data.frame (PARAMETER="MIND Filter (Individual with missing genotype)", VALUE=toString (params$MIND)))
	configDF = rbind  (configDF, data.frame (PARAMETER="GENO Filter (SNPs with missing genotype)", VALUE=toString (params$GENO)))
	configDF = rbind  (configDF, data.frame (PARAMETER="MAF Filter (Minor allele frequency)", VALUE=toString (params$MAF)))
	configDF = rbind  (configDF, data.frame (PARAMETER="HWE Filter (Hardy-Weinberg test)", VALUE=toString (params$HWE)))
	configDF = rbind  (configDF, data.frame (PARAMETER="GWAS Tools", VALUE=toString (params$tools)))

	outName = paste0(outputDir, "/out-multiGWAS-inputParameters.tbl")
	write.table (file=outName, configDF, quote=F, sep="\t", row.names=F)
	params$workingDir = getwd ()
	return (params)
}

#-----------------------------------------------------------
# Create a chord diagram for SNPs vs Chromosomes from
# summary table of best scores
#-----------------------------------------------------------
createChordDiagramSharedSNPs <- function (scoresFile) {
	# >>>>> local function: matrix and colors 
	createChord <- function (matrixChord=NULL, colorsChord=NULL) {
		if (is.null (matrixChord)){
			plot.new()
			mtext ("No chord diagram (without shared SNPs)")
		}else {
			chordDiagram(matrixChord, annotationTrack = "grid", directional = -1, direction.type = c("arrows"), # c("diffHeight", "arrows"),
						 grid.col = colorsChord, preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(matrixChord))))))

			# we go back to the first track and customize sector labels
			circos.track(track.index = 1, panel.fun = function(x, y) {
				xlim = get.cell.meta.data("xlim")
				xplot = get.cell.meta.data("xplot")
				ylim = get.cell.meta.data("ylim")
				sector.name = get.cell.meta.data("sector.index")							 

				#if(abs(xplot[2] - xplot[1]) < 20) {
				# Check if top or botton for Markers of Chromosomes
				if(abs(xplot[1]) < 180) 
					circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), col = "red")
				 else 
					circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), col = "blue")
				}, bg.border = NA) # here set bg.border to NA is important

				mtext ("Markers", side=3, col="red", cex=1.5)
				mtext ("Chromosomes", side=1, col="blue", cex=1.5)
		}
	} # >>>>> local function
	
	">>>>>> get shared SNPs <<<<<<>"
	getSharedSNPsFromFile <- function (scoresFile, N) {
		scores = read.table (file=scoresFile, header=T, sep="\t"); 
		summary = data.frame (add_count (scores, SNP, sort=T)); 
		sharedDups = summary [summary$n > 1,]
		shared = sharedDups [!duplicated (sharedDups$SNP),]
		return (shared)
	}
	">>>>>> get shared SNPs <<<<<<>"

	outFile = paste0(strsplit (scoresFile, split="[.]")[[1]][1], "-chordDiagram") 
	scores  = getSharedSNPsFromFile (scoresFile)

	# Check for shared SNPs
	matrixChord = NULL
	colorsChord = NULL
	if (nrow (scores) > 0) {
		# With shared SNPs
		tbl       = scores [,c("TOOL","CHROM","SNP")]

		# Group by TOOL and select 3 SNPs for each one
		#tblr = Reduce (rbind, by (tbl, tbl["TOOL"], head, n=2))
		tblm = tbl [,c(2,3)]
		chrs = sort (tblm [!duplicated (tblm [,1]), 1])
		snps = sort (as.character (tblm [!duplicated (tblm [,2]), 2]))

		# Create matrix Chroms X SNPs
		nChrs   = length (chrs)
		nSNPs   = length (snps)
		mat = as.data.frame (matrix (rep (0,nChrs*nSNPs),nrow=nChrs, ncol=nSNPs), stringAsFactor=F )
		rownames (mat) = chrs
		colnames (mat) = snps

		# Fill the matrix

		dmat = as.data.frame (mat)
		for (i in 1:nrow (tblm)) {
			chr = as.character (tblm [i, 1])
			snp = as.character (tblm [i, 2])
			dmat [chr,snp] = dmat [chr,snp] + 1
		}

		# Params for chordDiagram: mat and colors
		matrixChord = as.matrix (dmat) 

		colorSNPs <- colorRampPalette(brewer.pal(8, "Set2"))(length(chrs))

		#colorsChrs = setNames (brewer.pal (n=nChrs+3, name="RdBu"), c(chrs, "chXX", "chYY", "chZZ"))
		colorsChrs = setNames (colorSNPs, c(chrs))
		colorsSNPs = setNames (rep ("grey", nSNPs), snps)
		colorsChord = c(colorsChrs, colorsSNPs)
	}

	funCreateChords <- function () {
		createChord (matrixChord, colorsChord)
	}

	# PDF
	pdf (file=paste0 (outFile, ".pdf"), width=7, height=7)
		funCreateChords ()
	dev.off()
	#PNG
	png (file=paste0 (outFile, ".png"), width=7, height=7	, units="in", res=90)
		funCreateChords ()
	dev.off()
}

#----------------------------------------------------------
# Util to print head of data
# Fast head, for debug only
#----------------------------------------------------------
view <- function (data, m=10,n=10, tool="") {
	filename = paste0 ("x",tool,"-", deparse (substitute (data)),".csv")
	msgmsg (deparse (substitute (data)),": ", dim (data))
	if (is.null (dim (data)))
		print (data [1:10])
	else if (ncol (data) < 10) 
		print (data[1:m,])
	else if (nrow (data) < 10)
		print (data[,1:n])
	else 
		print (data [1:m, 1:n])

	write.table (file=filename, data, quote=F, sep="\t", row.names=F)
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
# Call to main function (first lines)
#-------------------------------------------------------------

#source ("lglib06.R")
#main ()

