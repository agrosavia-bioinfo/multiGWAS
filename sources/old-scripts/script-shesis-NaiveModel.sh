#!/bin/bash

# AUTHOR: Luis Garreta (lgarreta@agrosavia.co)
# Shesis with default parameters for quantitative traits

SHESISGENOPHENO=$1
SHESISMARKERSNAMES=$2
OUTFILE=$3

SHEsis --input $SHESISGENOPHENO --ploidy 4 --assoc --qtl --snpname-file $SHESISMARKERSNAMES --report-txt --adjust  --output $OUTFILE 

# Format output file by removing initial and end lines and replacing multiples tabs
tail -n +5 $OUTFILE.txt | head -n-1 > $OUTFILE.tmp
sed -r "s:\t+:\t:g" $OUTFILE.tmp > $OUTFILE.scores
