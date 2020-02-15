#!/bin/bash

# AUTHOR: Luis Garreta (lgarreta@agrosavia.co)
# Execute tassel pipeline for MLM model with default kinship and structure correction
# using principal components (N=5)

GENOPED=$1
GENOMAP=$2
PHENOTBL=$3
OUTFILE=$4

$MULTIGWAS_HOME/tools/run_pipeline.pl \
	-fork1 -plink -ped $GENOPED -map $GENOMAP \
	-fork2 -r $PHENOTBL \
	-fork3 -plink -ped $GENOPED -map $GENOMAP -PrincipalComponentsPlugin -covariance -reportEigenvalues false -reportEigenvectors false -endPlugin \
    -combine4 -input1 -input2 -input3 -intersect \
	-fork5 -plink -ped $GENOPED -map $GENOMAP -KinshipPlugin -method Centered_IBS -endPlugin \
    -combine6 -input4 -input5 \
	-mlm -mlmVarCompEst P3D \
	-mlmOutputFile $OUTFILE

