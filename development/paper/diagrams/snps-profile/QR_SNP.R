#!/usr/bin/Rscript

library (gplots)

#setwd("D:/Mis Documentos/Escritorio/papa/gwas")
geno_8019ACGT <- read.table("filtered-gwaspoly-genotype-ACGT.tbl", header = TRUE, sep = ",", na.strings = "NA", dec = ".", strip.white = TRUE)
geno_8019 <- read.table("filtered-gwaspoly-genotype-NUMERIC.tbl", header = TRUE, sep = ",", na.strings = "NA", dec = ".", strip.white = TRUE)
feno_8019 <- read.table("filtered-gwaspoly-phenotype.tbl", header = TRUE, sep = ",", na.strings = "NA", dec = ".", strip.white = TRUE)



marker=t(geno_8019[geno_8019$Marker=="c1_8019",])
#este archivo toca generarlo el cano 
#;0 -> 1;0;0;0;0 
#;1 -> 0;1;0;0;0
#;2 -> 0;0;1;0;0
#;3 -> 0;0;0;1;0
#;4 -> 0;0;0;0;1
#PAULA: por hacer - este archivo cano de manera automática

genocano8019<-read.csv("c1_8019_geno_cano.csv",sep=";")
genoxfeno<-genocano8019[,2:6]*feno_8019$tuber_shape
rownames(genoxfeno)<-genocano8019$Name
genoxfenov3<-as.matrix(genoxfeno[complete.cases(genoxfeno),])

#PAULA: por hacer
#Depende del SNP
# identificar caso 0 y caso 4 para ver casos en el medio y dejar en texto las opciones no 0 -4

nombres<-c("AAAA","AAAG","AAGG","AGGG","GGGG")
colnames(genoxfenov3)<-nombres
my_palette <- colorRampPalette(c("white", "black"))(n = 30)

heatmap.2(genoxfenov3,dendrogram = "row",trace=c("none"),Colv = FALSE,col=my_palette,cexCol = 1,cexRow = 0.6,
		  colsep=1:4, sepcolor="red", sepwidth=c(0.01,0.01), lhei=c(0.2,0.8), lwid=c(0.2,0.8), lmat=rbind (c(3,4),c(2,1)), key.title=NA )

