#!/bin/bash

# AUTHOR: Luis Garreta (lgarreta@agrosavia.co)
# Execute tassel pipeline for MLM model with default kinship and structure correction
# using principal components (N=5)

GENOPED=$1
GENOMAP=$2
PHENOTBL=$3
OUTFILE=$4

$TASSEL_HOME/run_pipeline.pl \
	-fork1 -plink -ped $GENOPED -map $GENOMAP \
	-fork2 -r $PHENOTBL \
	-fork4 -plink -ped $GENOPED -map $GENOMAP -PrincipalComponentsPlugin -covariance -reportEigenvalues false -reportEigenvectors false -endPlugin \
    -combine3 -input1 -input2 -input4 -intersect \
	-fork5 -plink -ped $GENOPED -map $GENOMAP -KinshipPlugin -method Dominance_Centered_IBS -endPlugin \
	-combine6 -input3 -input5 \
	-mlm -mlmVarCompEst P3D \
	-mlmOutputFile $OUTFILE

