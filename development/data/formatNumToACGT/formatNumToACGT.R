#!/usr/bin/Rscript

source ("gwas-formats.R")
#----------------------------------------------------------
# Main
#----------------------------------------------------------
#args = c ("agrosavia-genotype-checked-tetra-NUM.tbl", "potato_infinium_8303_map_context_DM_v3_superscaffolds-ref-alt-snps.tbl")
args = c ("agrosavia-genotype-checked-tetra-NUM.tbl", "agrosavia-genotype-map-altref.tbl")

genotypeFile  = args [1]
mapSNPs  = args [2]

numericToACGTFormatGenotype (genotypeFile, mapSNPs)

