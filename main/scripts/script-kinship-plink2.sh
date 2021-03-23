#!/bin/bash

# Calculate the related individuals until 2nd degree relationships
# It uses the kin software

PLINKPREFIXFILE=$1
OUTFILE=$2

# Call plink2 to filter related individuals 
# It is conventional to use a cutoff of ~0.354 (the geometric mean of 0.5 and 0.25) to screen for monozygotic twins and duplicate samples
# https://www.cog-genomics.org/plink/2.0/distance
cmm="plink2 --bfile $PLINKPREFIXFILE --allow-extra-chr --king-cutoff 0.177 --out $OUTFILE"
echo -e "\n>>>>" $cmm "\n"
eval $cmm

# Out filenames
outRelatedFile=$OUTFILE.king.cutoff.in.id
outUnRelatedFile=$OUTFILE.king.cutoff.out.id

# Remove lines from PED file
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
cmm="cp $PLINKPREFIXFILE.ped $OUTFILE.ped"
echo -e "\n>>>>" $cmm "\n"
eval $cmm
cmm="ln -s ../$PLINKPREFIXFILE.map $OUTFILE.map"
echo -e "\n>>>>" $cmm "\n"
eval $cmm

cmm="cat $outRelatedFile|xargs -I % sh -c  \"sed -i '/%/d' $OUTFILE.ped\""
echo -e "\n>>>>" $cmm "\n"
eval $cmm

