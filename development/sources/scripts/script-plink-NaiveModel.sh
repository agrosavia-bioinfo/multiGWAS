#!/bin/bash

# AUTHOR: Luis Garreta (lgarreta@agrosavia.co)
# Plink GWAS without correction for population structure 

GENOPLINK=$1
PHENOTBL=$2
OUTFILE=$3

plink --file $GENOPLINK --linear --adjust --pheno $PHENOTBL --all-pheno --allow-no-sex --out $OUTFILE

#plink --file $GENOPLINK --linear --adjust --pfilter 0.001 --pheno $PHENOTBL --all-pheno --allow-no-sex --out $OUTFILE
