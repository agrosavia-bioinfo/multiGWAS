#!/bin/bash

# AUTHOR: Luis Garreta (lgarreta@agrosavia.co)
# Plink GWAS without correction for population structure 

GENOPLINK=$1
PHENOTBL=$2
OUTFILE=$3
GENEACTION=$4

#plink --file $GENOPLINK --allow-extra-chr --linear --adjust --pheno $PHENOTBL --all-pheno --allow-no-sex --out $OUTFILE --linear $GENEACTION 
plink --file $GENOPLINK --allow-extra-chr --pheno $PHENOTBL --all-pheno --allow-no-sex --out $OUTFILE --linear $GENEACTION 

#cmm="plink --file $GENOPLINK ${TRAITTYPE} --allow-extra-chr --pheno $PHENOTBL --all-pheno --allow-no-sex --covar out-PCs.eigenvec --out $OUTFILE --linear $GENEACTION" 
#plink --file $GENOPLINK --linear --adjust --pfilter 0.001 --pheno $PHENOTBL --all-pheno --allow-no-sex --out $OUTFILE
