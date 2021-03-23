
args = commandArgs(trailingOnly = TRUE)
options (width=300)
#----------------------------------------------------------
# Util to print head of data
# Fast head, for debug only
#----------------------------------------------------------
view <- function (data, n=5,m=6) {
	filename = deparse (substitute (data))
	name = paste (deparse (substitute (data)),":  ")
	if (is.null (dim (data))) {
		dimensions = paste (length (data))
		message (name, class(data), " : (", paste0 (dimensions),")")
		if (length (data) < 6) n = length(data)
		print (data[1:n])
	}else {
		dimensions = paste0 (unlist (dim (data)),sep=c(" x ",""))
		message (name, class(data), " : (", paste0 (dimensions),")")
		if (n==0 | nrow (data) < 5) n = nrow(data)
		if (m==0 | ncol (data) < 6) m = ncol(data)
		print (data[1:n,1:m])
	}
	write.csv (data, paste0("x-", filename, ".csv"), quote=F, row.names=F)
}

#-------------------------------------------------------------
# Add label to filename and new extension (optional)
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
hd=view

#--------------------------------------------------------------
# Convert a dataframe to list of vectors
#--------------------------------------------------------------
dataFrameToList <- function (dataFrame) {
	dataList = split(dataFrame,seq(nrow(dataFrame)))
	dataList = setNames(dataList, rownames(dataFrame))
	return (dataList)
}
