	sep <- function (allele) {
		s="";
		for (i in 1:4) s=paste0(s, substr (allele,start=i,stop=i)," ");
		return (s)
	}

	genotypeFile  = "fltGeno.tbl"
	phenotypeFile = "fltPheno.tbl"
	geno    = read.csv (file=genotypeFile)
	pheno   = read.csv (file=phenotypeFile)
	rownames (pheno) <- pheno [,1]
	rownames (geno)  <- geno [,1]


	#markers = as.character (geno  [,1])
	#samples = as.character (pheno [,1])

	alleles    <- geno[,-c(1,2,3)]
	allelesMat <- t(sapply (alleles, sep))
	#allelesMat <- matrix (allelesSep, nrow=nrow(alleles), ncol=ncol(alleles), byrow=F)

	samples = rownames (allelesMat)
	markers = colnames (allelesMat)
	pheno   <- pheno [samples,]
	genoPhenoShesis    = data.frame (Sample=pheno[,1], Trait=pheno[,2],  allelesMat)

	outFile = "out/filtered-shesis-genopheno.tbl"
	write.table (file=outFile, genoPhenoShesis, quote=F,row.names=F,col.names=F, sep="\t")

	outFile = "out/filtered-shesis-markernames.tbl"
	write.table (file=outFile, rownames(geno), quote=F,row.names=F,col.names=F, sep="\t")
