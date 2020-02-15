#!/usr/bin/Rscript

source ("gwas-formats.R")
#----------------------------------------------------------
# Main
#----------------------------------------------------------
main <- function () 
{
	args = c ("agrosavia-genotype-checked.tbl")

	genotypeFile  = args [1]

	numericTetraToNumericDiploGenotype (genotypeFile)
}

main()

