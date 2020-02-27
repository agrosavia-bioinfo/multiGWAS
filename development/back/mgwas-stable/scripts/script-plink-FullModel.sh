#!/bin/bash

# AUTHOR: Luis Garreta (lgarreta@agrosavia.co)
# Plink GWAS with population structure using PCs (N=10)

GENOPLINK=$1
PHENOTBL=$2
OUTFILE=$3

plink --file $GENOPLINK --pca 10 --out out-PCs
plink --file $GENOPLINK --linear --adjust --pheno $PHENOTBL --all-pheno --allow-no-sex --covar out-PCs.eigenvec --out $OUTFILE
