#!/usr/bin/Rscript
# Recode a tetraploid genotipe (0,1,2,3) to diploid (0,1,2)
# also recode to ACGT (diplo) using info in solcap scafolds

# Remove duplicated SNPs

options (width=300)
library(parallel)

args = commandArgs(trailingOnly = TRUE)

args = c(genotype="agrosavia-genotype-CCC-CLEANED", snps="solcap-SNPs-refs-alts-AGs.tbl")
genotypeFile = args [1]
snpsFile     = args [2]
outName      = strsplit (genotypeFile, split="[.]")[[1]][1]

options(stringsAsFactors = FALSE)

#----------------------------------------------------------
# Reading files
#----------------------------------------------------------
message (">>> Reading genotype...")
genotypeDups = read.csv (genotypeFile, header=T)
genotype = genotypeDups [!duplicated(genotypeDups$Markers),]
dossages = genotype [,-c(1,2,3)]

message (">>> Reading SNPs")
SNPs     = read.table (snpsFile, header=T, row.names=1)
#----------------------------------------------------------
# To numeric diplo format (0,1,2)
#----------------------------------------------------------
dossages [dossages==2] = 1
dossages [dossages==3] = 1
dossages [dossages==4] = 2

newGenotype = cbind (genotype [,c(1,2,3)], dossages)

message (">>> Writing genotype diplo numeric")
write.csv (file="agrosavia-genotype-diplo-NUM.tbl", newGenotype, quote=F,row.names=F)

#----------------------------------------------------------
# To alleles diplo format (AA, AB, BB)
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
# Diplo ACGTs format using reference/alternate info chip
#----------------------------------------------------------
genoAGs = read.table (file="agrosavia-genotype-diplo-ABs.tbl", header=T, row.names=1)

n = nrow (genoAGs)
#genoAGs [,which (genoAGs [1,]=="AB")] = SNPs[genoAGs[,1], "AB"]

message (">>> Writing genotype diplo AGs")

for (id in rownames (genoAGs)) {
	genoAGs [id,which (genoAGs [id,]=="AA")] = SNPs[id, "AA"]
	genoAGs [id,which (genoAGs [id,]=="AB")] = SNPs[id, "AB"]
	genoAGs [id,which (genoAGs [id,]=="BB")] = SNPs[id, "BB"]
}

genoAGsDF = data.frame ("Markers"=rownames(genoAGs), genoAGs)
write.table (file="agrosavia-genotype-diplo-ACGT.tbl", genoAGsDF, quote=F,row.names=F, sep="\t")


#----------------------------------------------------------
# Tetra ACGTs format using reference/alternate info chip
#----------------------------------------------------------
changeAllele <- function (row, id, gAllele, sAllele) {
	geno [id,which (geno [id,]==gAllele)] = SNPs[id, sAllele]
}

dossages = genotype [,-c(1,2,3)]
rownames (dossages) = genotype [,1] 
message (">>> Writing genotype tetra AGs")
for (id in rownames (dossages)) {
	cat(".")
	dossages [id,which (dossages [id,]==0)] = SNPs[id, "AAAA"]
	dossages [id,which (dossages [id,]==1)] = SNPs[id, "AAAB"]
	dossages [id,which (dossages [id,]==2)] = SNPs[id, "AABB"]
	dossages [id,which (dossages [id,]==3)] = SNPs[id, "ABBB"]
	dossages [id,which (dossages [id,]==4)] = SNPs[id, "BBBB"]
}

newGenotype = cbind (genotype [,c(1,2,3)], dossages)

outFile = paste0 (outName,"-tetra-ACGT.tbl")
write.table (file=outFile, newGenotype, quote=F,row.names=F, sep="\t")

