#!/usr/bin/Rscript
# Recode a tetraploid genotipe (0,1,2,3) to diploid (0,1,2)
# also recode to ACGT (diplo) using info in solcap scafolds

# Remove duplicated SNPs

library(parallel)

options(stringsAsFactors = FALSE)

message (">>> Reading genotype...")
genotypeDups = read.csv ("agrosavia-genotype.tbl", header=T)
genotype = genotypeDups [!duplicated(genotypeDups$Markers),]
dossages = genotype [,-c(1,2,3)]

message (">>> Reading SNPs")
SNPs     = read.table ("SNPs-refs-alts-AGs.tbl", header=T, row.names=1)
#----------------------------------------------------------
# To numeric format (0,1,2)
#----------------------------------------------------------
dossages [dossages==2] = 1
dossages [dossages==3] = 1
dossages [dossages==4] = 2

newGenotype = cbind (genotype [,c(1,2,3)], dossages)

message (">>> Writing genotype diplo numeric")
write.csv (file="agrosavia-genotype-diplo-NUM.tbl", newGenotype, quote=F,row.names=F)

#----------------------------------------------------------
# To alleles format (AA, AB, BB)
#----------------------------------------------------------
dossages = genotype [,-c(1,2,3)]

dossages [dossages==0] = "AA"
dossages [dossages==1] = "AB"
dossages [dossages==2] = "AB"
dossages [dossages==3] = "AB"
dossages [dossages==4] = "BB"

newGenotype = cbind (genotype [,c(1,2,3)], dossages)

message (">>> Writing genotype diplo ABs")
write.table (file="agrosavia-genotype-diplo-ABs.tbl", newGenotype, quote=F,row.names=F,sep="\t")

#----------------------------------------------------------
# To ACGTs format using reference/alternate info chip
#----------------------------------------------------------
genoAGs = read.table (file="agrosavia-genotype-diplo-ABs.tbl", header=T, row.names=1)

n = nrow (genoAGs)
#genoAGs [,which (genoAGs [1,]=="AB")] = SNPs[genoAGs[,1], "AB"]

message (">>> Writing genotype diplo AGs")
#changeRow <- function (id) {
	#genoAGs [id,which (genoAGs [id,]=="AA")] <<- SNPs[id, "AA"]
	#genoAGs [id,which (genoAGs [id,]=="AB")] <<- SNPs[id, "AB"]
	#genoAGs [id,which (genoAGs [id,]=="BB")] <<- SNPs[id, "BB"]
#}
#mclapply (rownames(genoAGs), changeRow, mc.cores=4)

for (id in rownames (genoAGs)) {
	genoAGs [id,which (genoAGs [id,]=="AA")] = SNPs[id, "AA"]
	genoAGs [id,which (genoAGs [id,]=="AB")] = SNPs[id, "AB"]
	genoAGs [id,which (genoAGs [id,]=="BB")] = SNPs[id, "BB"]
}

genoAGsDF = data.frame ("Markers"=rownames(genoAGs), genoAGs)
write.table (file="agrosavia-genotype-diplo-AGs.tbl", genoAGsDF, quote=F,row.names=F, sep="\t")


