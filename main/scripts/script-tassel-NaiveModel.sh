#!/bin/bash

# AUTHOR: Luis Garreta (lgarreta@agrosavia.co)
# Execute tassel pipeline for GLM model only with genotype and phenotype
# No population structure correction (Kinship or PCs)

GENOVCF=$1
PHENOTBL=$2
OUTFILE=$3
$MULTIGWAS_HOME/opt/tools/run_pipeline.pl \
	-fork1 -vcf $GENOVCF \
	-fork2 -r $PHENOTBL \
	-combine3 -input1 -input2 -intersect \
	-FixedEffectLMPlugin -appendAddDom -endPlugin \
	-export $OUTFILE-
