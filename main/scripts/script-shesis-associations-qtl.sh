#!/bin/bash

# AUTHOR: Luis Garreta (lgarreta@agrosavia.co)
# Shesis with default parameters for quantitative traits

SHESISGENOPHENO=$1
PLOIDY=$2
SHESISMARKERSNAMES=$3
OUTFILE=$4
FLAGQTL=$5

SHEsis --input $SHESISGENOPHENO --ploidy $PLOIDY --assoc --snpname-file $SHESISMARKERSNAMES --report-txt --adjust  --output $OUTFILE $FLAGQTL 

# Format output file by removing initial and end lines and replacing multiples tabs
#tail -n +5 $OUTFILE.txt | head -n-1 > $OUTFILE.tmp
#sed -r "s:\t+:\t:g" $OUTFILE.tmp > $OUTFILE.csv
