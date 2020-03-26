#!/bin/bash

# AUTHOR: Luis Garreta (lgarreta@agrosavia.co)
# Execute tassel pipeline for MLM model with default kinship and structure correction
# using principal components (N=5)

GENOVCF=$1
PHENOTBL=$2
OUTFILE=$3

$MULTIGWAS_HOME/tools/run_pipeline.pl \
	-fork1 -vcf $GENOVCF \
	-fork2 -r $PHENOTBL \
	-fork3 -vcf $GENOVCF -PrincipalComponentsPlugin -covariance -reportEigenvalues false -reportEigenvectors false -endPlugin \
    -combine4 -input1 -input2 -input3 -intersect \
	-fork5 -vcf $GENOVCF -KinshipPlugin -method Centered_IBS -endPlugin \
    -combine6 -input4 -input5 \
	-mlm -mlmVarCompEst P3D \
	-mlmOutputFile $OUTFILE

