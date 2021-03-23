#!/usr/bin/Rscript

parametersLines = readLines ("full.config");parametersLines

params=list()
for (line in  parametersLines) {
	if (line=="") next
	parVal = strsplit (line, ":")[[1]]
	par = trimws (parVal[1], whitespace="[ \"\t\r\n]")
	val = trimws (parVal[2], whitespace="[ \"\t\r\n]")
	names (val) = par
	params = append (params, val)
}

# Check possible errors
`%notin%` = Negate (`%in%`)
if (params$ploidy %notin% c("2" or "4"))
	stop ("Ploidy not supported")

if (!file.exists (params$genotypeFile))
	stop ("Genotype file not found")

if (!file.exists (params$phenotypeFile))
	stop ("Phenotype file not found")

if (params$genotypeFormat %in% c("kmatrix", "fitpoly") && !file.exists (param$mapFile))
	stop ("Map file not found")


		


