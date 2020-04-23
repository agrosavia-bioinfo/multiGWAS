#!/usr/bin/Rscript

#----------------------------------------------------------
# Util to print head of data
# Fast head, for debug only
#----------------------------------------------------------
hd <- function (data, m=10,n=10) {
	message (deparse (substitute (data)),":")
	if (is.null (dim (data)))
		print (data [1:10])
	else if (ncol (data) < 10)
		print (data[1:m,])
	else if (nrow (data) < 10)
		print (data[,1:n])
	else
		print (data [1:m, 1:n])
}

#

args = commandArgs(trailingOnly = TRUE)
args = c("Positions_markers_chip_8K_papa_jberdugo.txt")

SNPsFile = args [1]

SNPs = read.table (file=SNPsFile, header=T)
hd (SNPs[,2])

allelesList = strsplit (as.character (SNPs[,2]), split="/")
splitf  <- funcion (lst) {
	return (lst[[1]][1]
print (allelesList)
alleles = list(REF=allelesList[[1]][1], ALT=allelesList[[1]][2])
print (alleles)



