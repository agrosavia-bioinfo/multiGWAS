#!/bin/bash

# Calculate the related individuals until 2nd degree relationships
# It uses the kin software

PLINKPREFIXFILE=$1
OUTFILE=$2

# Convert plink file to binary
cmm="plink --file $PLINKPREFIXFILE --make-bed --out $PLINKPREFIXFILE"
echo -e "\n>>>>" $cmm "\n"
eval $cmm

# Call plink2 to filter related individuals 
# It is conventional to use a cutoff of ~0.354 (the geometric mean of 0.5 and 0.25) to screen for monozygotic twins and duplicate samples
# https://www.cog-genomics.org/plink/2.0/distance
cmm="plink2 --bfile $PLINKPREFIXFILE --king-cutoff 0.354 --out $OUTFILE"
echo -e "\n>>>>" $cmm "\n"
eval $cmm

## Call kin to search for related individuals
## king -b $PLINKPREFIXFILE.bed --unrelated --degree 2 --prefix $OUTFILE-

# Out filenames
outRelatedFile=$OUTFILE.king.cutoff.out.id
outUnRelatedFile=$OUTFILE.king.cutoff.in.id

# Remove lines from PED file
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
cmm="cp $PLINKPREFIXFILE.ped $PLINKPREFIXFILE-kinship.ped"
echo -e "\n>>>>" $cmm "\n"
eval $cmm
cmm="ln -s ../$PLINKPREFIXFILE.map $PLINKPREFIXFILE-kinship.map"
echo -e "\n>>>>" $cmm "\n"
eval $cmm

cmm="cat $outRelatedFile|xargs -I % sh -c  \"sed -i '/%/d' $PLINKPREFIXFILE-kinship.ped\""
echo -e "\n>>>>" $cmm "\n"
eval $cmm

## Convert plink file to binary
#cmm="plink --bfile $PLINKPREFIXFILE-kinship --recode tab --out $PLINKPREFIXFILE-kinship"
#echo -e "\n>>>>" $cmm "\n"
#eval $cmm
