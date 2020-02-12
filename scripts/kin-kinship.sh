#!/bin/bash

# Calculate the related individuals until 2nd degree relationships
# It uses the kin software

PLINKPREFIXFILE=$1
OUTPREFIX=$2

# Convert plink file to binary
plink --file $PLINKPREFIXFILE --make-bed --out $PLINKPREFIXFILE

# Call kin to search for related individuals
king -b $PLINKPREFIXFILE.bed --unrelated --degree 2 --prefix $OUTPREFIX-

# Out filenames
outRelatedFile=$OUTPREFIX-unrelated_toberemoved.txt
outUnRelatedFile=$OUTPREFIX-unrelated.txt

# Remove lines from PED file
cat $outRelatedFile|xargs -I % sh -c  "sed -i '/%/d' $PLINKPEDFILE"

# Convert plink file to binary
plink --bfile $PLINKPREFIXFILE --recode tab --out $PLINKPREFIXFILE
