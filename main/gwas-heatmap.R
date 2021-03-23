#!/usr/bin/Rscript
# LOG: 
#     Jul 11: Removed NAs with complete.cases
#     Jun 06: Fixed, not include NAs in both geno and pheno

suppressMessages (library (gplots)) 
#----------------------------------------------------------
# Main
#----------------------------------------------------------
main <- function () {
	library ("parallel")
	args = commandArgs (trailingOnly=T)

	genoFileACGT <- "out/filtered-gwasp4-genotype.tbl"
	genoFileNUM  <- "out/filtered-gwasp4-genotype-NUM.tbl"
	phenoFile    <- "out/filtered-gwasp4-phenotype.tbl"
	snpList      <- c("c2_51250","c1_7770", "c1_2978")
	ploidy       <- 4
	createHeatmapForSNPList ("./",genoFileACGT, genoFileNUM, phenoFile, snpList, ploidy) 
}

#----------------------------------------------------------
# create SNP profiles for a list of SNPS
#----------------------------------------------------------
createHeatmapForSNPList <- function (outputDir, genoFileACGT, genoFileNUM, phenoFile, snpList, ploidy) {
	outName = paste0 (outputDir, "/out-SNPProfile") 

	genotypeACGT    <- read.csv (genoFileACGT, na.strings = "NA", dec = ".", strip.white = TRUE, check.names=F )
	genotypeNUMERIC <- read.csv (genoFileNUM,  na.strings = "NA", dec = ".", strip.white = TRUE, check.names=F)
	phenotype       <- read.csv (phenoFile,    na.strings = "NA", dec = ".", strip.white = TRUE, check.names=F)

	pdfHeatMap <- function (snp) {
		msgmsg  ("\tHeatmap for snp: ", snp)
		createHeatmapForSNP (outputDir, genotypeACGT, genotypeNUMERIC, phenotype, snp, ploidy)
	}

	# Not works in parallel
	NCORES = detectCores ()
	for (i in snpList) 
		pdfHeatMap (i)
	#res=mclapply (snpList, pdfHeatMap, mc.cores=1)
}

#----------------------------------------------------------
# create a SNP profile for a snpId
#----------------------------------------------------------
createHeatmapForSNP <- function (outputDir, genotypeACGT, genotypeNUMERIC, phenotype, snpId, ploidy) {

	# Get names and values for phenotype
	phenoNames  = phenotype [,1]
	phenoValues = phenotype [,2]
	trait       = colnames (phenotype) [2]

	# Get samples for marker in both matrices: Numeric and ACGT
	samplesMarker          = t(genotypeNUMERIC[genotypeNUMERIC[,1]==snpId,])
	samplesMarkerMatrixNUM = as.matrix(samplesMarker[-1:-3,])

	samplesMarker           = t(genotypeACGT[genotypeACGT[,1]==snpId,])
	samplesMarkerMatrixACGT = as.matrix(samplesMarker[-1:-3,])

	# Count genotype types
	one.hot<-function(sqnce, alphabet){
	  y<-unlist(strsplit(sqnce,""))
	  sapply(y,function(x){match(alphabet,x,nomatch=0)})
	}

	alphabet_AA = as.character (0:ploidy)

	c           <- t(one.hot (samplesMarkerMatrixNUM, alphabet = alphabet_AA))
	names(c)    <- names(samplesMarkerMatrixNUM)
	e           <- as.data.frame(c,row.names = samplesMarkerMatrixNUM[1,])
	e$Name      <- row.names(samplesMarkerMatrixNUM)
	genoxfeno_2 <- subset(e,e$Name%in%phenoNames)

	genoxfeno_3 <- genoxfeno_2[,1:(ploidy+1)]*phenoValues

	# Remove rows with NAs from genoxfeno_3 
	g3NAs                 = genoxfeno_3
	rownames (g3NAs)      = phenoNames 
	g3noNAs               = g3NAs[complete.cases(g3NAs),]
	genoxfenov4           = as.matrix(g3noNAs)
	rownames(genoxfenov4) = rownames (g3noNAs)

	data_geno_NUM_ACGT<-data.frame(samplesMarkerMatrixNUM,samplesMarkerMatrixACGT)


	# Get genotype names (e.g AAAA,..,AAGG,...,GGGG)
	namesGen = c()
	for (i in 0:ploidy) {
		gen = subset (data_geno_NUM_ACGT, samplesMarkerMatrixNUM==i,select = samplesMarkerMatrixACGT)
		g   = as.character (gen$samplesMarkerMatrixACGT[1])    
		namesGen = c(namesGen, g)
	}
	colnames(genoxfenov4) = namesGen

	#my_palette <- colorRampPalette(c("white", "blue", "darkblue", "darkred"))(n = 100)
	#my_palette <- colorRampPalette(c("white", "blue", "darkblue", "darkred"))(n = 100)
	my_palette <- colorRampPalette(c("black", "blue", "red"))(n = 50)
	lmat <- rbind(c(5,4), c(2,3), c(2,1))
	lhei <- c(12,0.1,32)
	lwid <- c(3,9)

	myplot<- function() {
	  oldpar <- par("mar")
	  hist(phenoValues, main = "Histogram", xlab=trait)
	}

	# Function that call heatmap.2 (only to save code)
	fun_heatmap <- function () {
		hmap = heatmap.2(genoxfenov4,
						 #cellnote=genoxfenov4,
						 dendrogram = "row",
						 reorderfun=function(snpId, w) reorder(snpId, w, agglo.FUN = mean),
						 Colv=FALSE,
						 adjCol = c(NA,0),
						 key=TRUE,
						 srtCol=360,
						 col=my_palette,
						 cexCol = 1,
						 lmat=lmat,
						 lhei=lhei, 
						 lwid=lwid, 
						 extrafun=myplot,
						 key.xlab=paste0 ("Value of ", trait), 
						 xlab=snpId, key.title="Color Key",
						 #rowsep=c(0:nrow(genoxfenov4)),
		                 #sepcolor = "red",
						 #sepwidth = c(0.001,0.001)
						 colsep=c(0:ploidy+1)
						 )
	}

	outName = paste0 (outputDir, "/out-SNPProfile") 
	pdf (paste0(outName,"-", snpId, ".pdf"), width=7, height=7)
		fun_heatmap ()
	dev.off()

	png (paste0(outName,"-",snpId, ".png"), width=7, height=7, units="in", res=72)
		fun_heatmap ()
	dev.off()
}

#-------------------------------------------------------------
# Print a log message with the parameter
#-------------------------------------------------------------
msg <- function (...) 
{
		messages = unlist (list (...))
		cat (">>>>", messages, "\n")
}

#----------------------------------------------------------
#----------------------------------------------------------
#main()
