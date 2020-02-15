
# For diplo
#SHEsis --input geno-diplo.shesis --ploidy 2 --assoc --qtl --snpname-file snpnames.shesis --report-txt --adjust  --output out-shesis-adj-AGs 

GENOPHENO=$1
NAMES=$2
# For tetra
cmm="shesis --input $GENOPHENO  --ploidy 4 --assoc --qtl --snpname-file $NAMES --report-txt --adjust  --output out-shesis-naive-tetra.tbl"
echo $cmm
eval $cmm
