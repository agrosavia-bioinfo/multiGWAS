#!/bin/bash

# AUTHOR: Luis Garreta (lgarreta@agrosavia.co)
# Execute tassel pipeline for GLM model only with genotype and phenotype
# No population structure correction (Kinship or PCs)

GENOPED=$1
GENOMAP=$2
PHENOTBL=$3
OUTFILE=$4
$MULTIGWAS_HOME/tools/run_pipeline.pl \
	-fork1 -plink -ped $GENOPED -map $GENOMAP \
	-fork2 -r $PHENOTBL \
	-combine3 -input1 -input2 -intersect \
	-FixedEffectLMPlugin -endPlugin \
	-export $OUTFILE-
