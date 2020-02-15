#!/usr/bin/Rscript

# Filter and write plink phenotype according to .ped filtered file

library (dplyr)
options (width=300, stringsAsFactors=T)
args = commandArgs(trailingOnly=T)

args = c("pheno.tbl", "agrosavia-genotype-FLT-plink.ped")
phenoFile = args [1]
pedFile   = args [2]

pheno = read.table (file=phenoFile, header=T)
ped   = read.table (file=pedFile, header=F)

phenoFLT = pheno %>% filter (IID %in% ped[,2]) %>% arrange (match (IID, ped[,2]))

write.table (file="pheno-FLT.tbl", phenoFLT, row.names=F, quote=F, sep="\t")


