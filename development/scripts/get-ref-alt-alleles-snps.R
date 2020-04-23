#!/usr/bin/Rscript

# Get a table with Reference/Alternate alleles for each SNP
# Remove duplicated SNPs

args = commandArgs (trailingOnly=T)

args = c ("potato_infinium_8303_map_context_DM_v3_superscaffolds.txt")

inFile = args [1]
name = strsplit (inFile, split="\\.")[[1]][1]
#outFile = "snps-ref-alt-nucleotides.tbl"

snpsAll   = read.table (inFile, header=T)
snps      = snpsAll [!duplicated (snpsAll$Name),]
ids       = snps [1]
seq       = snps [4]
alleles   = unlist(strsplit (substr (seq[,1], 52,54), split="/"))
allelesDF = data.frame (matrix (alleles, ncol=2, byrow=T))
ref      = as.character (allelesDF [,1])
alt      = as.character (allelesDF [,2])

AA       = paste0 (ref, ref)
AB       = paste0 (ref, alt)
BB       = paste0 (alt, alt)

AAAA     = paste0 (ref,ref,ref,ref)
AAAB     = paste0 (ref,ref,ref,alt)
AABB     = paste0 (ref,ref,alt,alt)
ABBB     = paste0 (ref,alt,alt,alt)
BBBB     = paste0 (alt,alt,alt,alt)


out = cbind (snp=ids, ref=ref, alt=alt, AA=AA, AB=AB, BB=BB, AAAA=AAAA, AAAB=AAAB, AABB=AABB, ABBB=ABBB, BBBB=BBBB)
outFilename = sprintf ("%s-ref-alt-snps.tbl", name)
write.table (file=outFilename, out, quote=F, row.names=F, sep="\t")


