#!/usr/bin/Rscript
library(qqman)

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
#-------------------------------------------------------------
manhattanSimple <- function (filename) {
	msg ("File:", filename)
	newFile =  gsub (".scores", "-manhattan.pdf", filename)
	msg ("newFile:", newFile)
	data <- read.table (file=filename, header=T)
	hd (data)
	pdf (newFile)
		names = unlist (strsplit (basename (filename), "[-|.]"))
		mainTitle = paste0 (names[2],"-", names [3], " manhattan plot")

		if (grepl ("Gwaspoly", filename)) 
			manhattan(data,col = c("blue4", "orange3"),  snp="Marker", chr="Chrom", bp="Position", p="P", 
				  suggestiveline=F,genomewideline=F, main=mainTitle)
		else
			manhattan(gwasResults,col = c("blue4", "orange3"),  snp="SNP", chr="CHR", bp="POS", p="P", 
				  suggestiveline=F,genomewideline=F, main=mainTitle)
	dev.off()
}
#-------------------------------------------------------------
#-------------------------------------------------------------
inputDir  = "."
gwasType  = "Naive"
files     =  list.files(inputDir, pattern=paste0("^(.*(",gwasType,").*(scores)[^$]*)$"), full.names=T)
modelsLst = list ("general"=0,"additive"=10,"1-dom-alt"=20, "1-dom-ref"=30, "2-dom-alt"=40, "2-dom-ref"=50)
models    = c("general", "additive", "1-dom-alt", "1-dom-ref", "2-dom-alt", "2-dom-ref")
#data.frame (c("general", "additive, "1-dom-alt", "1-dom-ref", "2-dom-alt", "2-dom-ref"), 0, 10, 20, 30, 40, 50)

for (filename in files) {
	manhattanSimple (filename)
}
