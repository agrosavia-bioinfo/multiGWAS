#!/usr/bin/python

import os, sys

phenotypesFile = "papa-agrosavia-fenotipoGota.txt"
annotationFile = "potato_8303SNPs_potato_dm_v4.03.gff3"
shipFile      = "potato_infinium_8303_map_context_DM_v3_superscaffolds.txt"

#----------------------------------------------------------
#----------------------------------------------------------
def main ():
    gotaPhenotypeDic    = getDicFromPhenotypes (phenotypesFile)
    genomeAnnotationDic = getDicFromAnnotations (annotationFile)
    shipArraySNPsDic    = getDicFromShip (shipFile)

    createPEDFile (gotaPhenotypeDic, genomeAnnotationDic, shipArraySNPsDic)
    
#----------------------------------------------------------
#----------------------------------------------------------

#----------------------------------------------------------
#----------------------------------------------------------
def getDicFromShip (shipFile):
    linesList = open (shipFile).readlines()[1:]
    dic = {}

    outFileErrors = open ("agrosavia-ship-annotations.errors", "w")
    for line in linesList:
        fields      = line.split ()
        snpId       = fields [0]
        position    = fields [2]
        sequence    = fields [3]
        alleles     = getAlleles (sequence)

        repeated      = dic.get (snpId)
        if repeated != None:
            outFileErrors.write (line)
            outFileErrors.write ("%s\n" % repeated)

        dic [snpId] = [position, alleles]

    outFile = open ("agrosavia-ship-annotations.tbl", "w")
    outFile.write ("# snpID, Position, Alleles\n")
    for key in sorted (dic.keys()):
        outFile.write ("%s\t%s\t%s\n" % (key, dic [key][0], dic [key][1]))

def getAlleles (sequence):
    allele1 = sequence.split ("[")[1].split("/")[0]
    allele2 = sequence.split ("/")[1].split("]")[0]

    return (allele1+allele2)


#----------------------------------------------------------
#----------------------------------------------------------
def getDicFromAnnotations (annotationFile):
    linesList = open (annotationFile).readlines()
    dic = {}

    outFileErrors = open ("agrosavia-genome-annotations.errors", "w")
    for line in linesList:
        fields        = line.split ()
        chromosome    = fields [0]
        position      = fields [3]

        snpName       = fields[-1].split("=")[2].strip(";")
        repeated      = dic.get (snpName)
        if repeated != None:
            outFileErrors.write (line)
            outFileErrors.write ("%s\n" % repeated)

        dic [snpName] = [chromosome,position]

    outFile = open ("agrosavia-genome-annotations.tbl", "w")
    outFile.write ("# Name, Chromosome, Position\n")
    for key in sorted (dic.keys()):
        outFile.write ("%s\t%s\t%s\n" % (key, dic [key][0], dic[key][1]))
    outFile.close ()

    return dic

#----------------------------------------------------------
#----------------------------------------------------------
def getDicFromPhenotypes (phenotypesFile):
    dic = {}
    outFile       = open ("agrosavia-fenotipo-gota.tbl", "w")
    outFileErrors = open ("agrosavia-fenotipo-gota.errors", "w")
    

    for line in open(phenotypesFile).readlines()[1:]:
        try:
            fields         = line.split ("\t")
            sampleName     = fields[0]
            sampleId       = int (sampleName.split ("_")[1])
            gotaValue      = fields[5] 
            dic [sampleId] = [sampleName, gotaValue]
        except:
            outFileErrors.write (line)

    outFile = open ("agrosavia-fenotipo-gota.tbl", "w")
    outFile.write ("# SampletId, GotaValue\n")
    for key in sorted (dic.keys()):
        outFile.write ("%s\t%s\n" % (dic [key][0], dic[key][1]))
    outFile.close ()

    return dic

#----------------------------------------------------------
#----------------------------------------------------------
main ()



