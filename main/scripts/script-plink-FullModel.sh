#!/bin/bash

# AUTHOR: Luis Garreta (lgarreta@agrosavia.co)
# Plink GWAS with population structure using PCs (N=10)

GENOPLINK=$1
TRAITTYPE=$2
PHENOTBL=$3
OUTFILE=$4
GENEACTION=$5

cmm="plink --file $GENOPLINK --allow-extra-chr --pca 5 --out out-PCs"
echo ">>>> cmm1: " $cmm
eval $cmm
cmm="plink --file $GENOPLINK ${TRAITTYPE} --allow-extra-chr --pheno $PHENOTBL --all-pheno --allow-no-sex --covar out-PCs.eigenvec --out $OUTFILE --linear $GENEACTION" 
echo ">>>> cmm2: " $cmm
eval $cmm
