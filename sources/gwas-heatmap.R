#!/usr/bin/Rscript
suppressMessages (library (gplots)) 
#----------------------------------------------------------
# Main
#----------------------------------------------------------
main <- function () {
	args = commandArgs (trailingOnly=T)

	genoFileACGT <- "filtered-gwaspoly-genotype-ACGT.tbl"
	genoFileNUM  <- "filtered-gwaspoly-genotype-NUMERIC.tbl"
	phenoFile    <- "filtered-gwaspoly-phenotype.tbl"
	snpList      <- c ("c1_8019", "c2_21314")         # SNPs a visualizar
	createHeatmapForSNPList ("./",genoFileACGT, genoFileNUM, phenoFile, snpList) 
}

#----------------------------------------------------------
# createHeatmapForSNP
#----------------------------------------------------------
createHeatmapForSNPList <- function (outputDir, genoFileACGT, genoFileNUM, phenoFile, snpList) {
	outName = paste0 (outputDir, "/out-SNPProfile") 
	for (snp in snpList) {
		#msgmsg ("    >>>> Plotting SNP profile for SNP: ", snp)
		pdf (paste0(outName,"-",snp, ".pdf"), width=7, height=7)
		createHeatmapForSNP (outputDir, genoFileACGT, genoFileNUM, phenoFile, snp)
		dev.off()

		png (paste0(outName,"-",snp, ".png"), width=7, height=7, units="in", res=72)
		createHeatmapForSNP (outputDir, genoFileACGT, genoFileNUM, phenoFile, snp)
		dev.off()
	}
}

createHeatmapForSNP <- function (outputDir, genoFileACGT, genoFileNUM, phenoFile, snpId) {
	genotypeACGT    <- read.table(genoFileACGT, header = TRUE, sep = ",", na.strings = "NA", dec = ".", strip.white = TRUE)
	genotypeNUMERIC <- read.table(genoFileNUM, header = TRUE, sep = ",", na.strings = "NA", dec = ".", strip.white = TRUE)
	phenotype       <- read.table(phenoFile, header = TRUE, sep = ",", na.strings = "NA", dec = ".", strip.white = TRUE)

	phenoNames  = phenotype [,1]
	phenoValues = phenotype [,2]
	trait = colnames (phenotype) [2]


	marker=t(genotypeNUMERIC[genotypeNUMERIC$Marker==snpId,])
	marker_2<-as.matrix(marker[-1:-3,])
	head(marker_2)
	alphabet_AA<-c("0","1","2","3","4")

	one.hot<-function(sqnce, alphabet){
	  y<-unlist(strsplit(sqnce,""))
	  sapply(y,function(x){match(alphabet,x,nomatch=0)})
	}

	a <- marker_2
	c <- t(one.hot(a, alphabet = alphabet_AA))
	names(c)<-names(marker_2)
	e<-as.data.frame(c,row.names = marker_2[1,])
	e$Name<-row.names(marker_2)
	e
	genoxfeno_2<-subset(e,e$Name%in%phenoNames)
	genoxfeno_3<-genoxfeno_2[,1:5]*phenoValues
	genoxfeno_3
	head(phenoValues)
	genoxfenov4<-as.matrix(genoxfeno_3[complete.cases(genoxfeno_3),])
	marker=t(genotypeACGT[genotypeACGT$Marker==snpId,])
	marker_3<-as.matrix(marker[-1:-3,])
	head(genotypeACGT)
	marker_3
	z<-data.frame(marker_3)
	data_geno_num_ACGT<-data.frame(marker_2,marker_3)
	data_geno_num_ACGT
	gen_0<-subset(data_geno_num_ACGT,marker_2=="0",select = marker_3)
	o<-as.character(gen_0$marker_3[1])    
	gen_4<-subset(data_geno_num_ACGT,marker_2=="4",select = marker_3)
	a<-as.character(gen_4$marker_3[1])
	gen_1<-subset(data_geno_num_ACGT,marker_2=="1",select = marker_3)
	l<-as.character(gen_1$marker_3[1])
	gen_2<-subset(data_geno_num_ACGT,marker_2=="2",select = marker_3)
	s<-as.character(gen_2$marker_3[1])
	gen_3<-subset(data_geno_num_ACGT,marker_2=="3",select = marker_3)
	E<-as.character(gen_3$marker_3[1])
	colnames(genoxfenov4)<-c(o,l,s,E,a)
	rownames(genoxfenov4)<-phenoNames
	my_palette <- colorRampPalette(c("white", "black"))(n = 30)
	lmat <- rbind(c(5,4), c(2,3), c(2,1))
	lhei <- c(12,0.1,32)
	lwid <- c(3,9)

	myplot<- function() {
	  oldpar <- par("mar")
	  hist(phenoValues, main = "Histogram", xlab=trait)
	}
	my_palette <- colorRampPalette(c("white", "black"))(n = 30)

	heatmap.2(genoxfenov4,dendrogram = "row",reorderfun=function(snpId, w) reorder(snpId, w, agglo.FUN = mean),
			  Colv=FALSE,adjCol = c(NA,0),key=TRUE,srtCol=360,
			  col=my_palette,cexCol = 1,lmat=lmat,lhei=lhei, lwid=lwid, extrafun=myplot, 
			  key.xlab=paste0 ("Value of ", trait), xlab=snpId, key.title="Color Key")
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
