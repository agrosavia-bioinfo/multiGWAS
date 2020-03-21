#!/usr/bin/python

# For genotype transformed to VCF
# Takes as input a tassel phenotype and convert names
# as family_individual by duplicating name

import sys

args = sys.argv

pheno = args [1]

individuals = open (pheno).readlines ()
print individuals[0].strip()
for line in individuals[1:]:
    name  = line.split()[0]
    value = line.split()[1]
    newName = name+"_"+name
    newLine = newName +"\t" + value
    print newLine
