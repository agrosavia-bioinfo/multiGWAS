#!/usr/bin/Rscript

options(java.parameters = c("-Xmx1g", "-Xms1g"))
library(rTASSEL)

rTASSEL::startLogger(fullPath = NULL, fileName = NULL)

# Load in hapmap file
message ("")
message (">>>>> Genotype")
geno <- "tassel-geno.vcf"
tasGenoHMP <- rTASSEL::readGenotypeTableFromPath(path=geno)
print (tasGenoHMP)

# Load into rTASSEL "GenotypePhenotype" object
message ("")
message (">>>>> Phenotype")
pheno="phenot.tbl"
tasPheno <- rTASSEL::readPhenotypeFromPath(path=pheno)
print (tasPheno)

message ("")
message (">>>> Geno/Pheno")
tasGenoPheno <- rTASSEL::readGenotypePhenotype(
	genoPathOrObj = tasGenoHMP,
    phenoPathDFOrObj = tasPheno
)

print (tasGenoPheno)


message ("")
trait = "tuber_shape"
message ("Calculating GLM...")
tasGLM <- rTASSEL::assocModelFitter(tasObj=tasGenoPheno, formula=list(tuber_shape)~., fitMarkers=T, kinship=NULL, fastAssociation=F)
#
print (tasGLM)
write.table (tasGLM$GLM_Stats, file="tasse_glm.scores", quote=F, col.names=T, row.names=F, sep="\t")
