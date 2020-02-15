#!/usr/bin/Rscript

library (stringr)
library (dplyr)
options (width=300)
options(scipen=999)

model="Naive"

files =  list.files(".", pattern="^(.*(Naive).*(tbl)[^$]*)$")


files = c("out-Gwasp4-Naive-significativeQTLs.tbl")
#files = c("out-Gwasp4-Naive-significativeQTLs.tbl", "out-Plink-Naive-assoc.linear.adjusted.tbl", "out-Tassel-Naive-GLM_Stats_geno+pheno.tbl")
summTable = data.frame ()
for (f in files) {
    data = read.table (file=f, header=T)
    print (f)
    if (str_detect(f, "Gwasp4")) {
		tool   = "Gwasp4"
		snps   = data$Marker
		pVal   = round (10^(-data$Score),10)
		chrom  = data$Chrom
		pos    = data$Position
		df = data.frame (TOOL=tool, MODEL=model, SNPs=snps, P = pVal, CHR=chrom, POS=pos,pTHR=pThreshold)
    }else if (str_detect (f, "Gwasp2")) {
		tool   = "Gwasp2"
		snps   = data$Marker
		pVal   = round (10^(-data$Score),10)
		chrom  = data$Chrom
		pos    = data$Position
	}else if (str_detect (f, "Plink")) {
		tool   = "Plink"
		snps   = data$SNP
		pVal   = data$BONF
		chrom  = data$CHR
		pos    = NA
	}else if (str_detect (f, "Tassel")) {
		tool   = "Tassel"
		snps   = data$Marker
		pVal   = data$p
		chrom  = data$Chr
		pos    = data$Pos
	}
    	pThreshold = round (0.05/(578*2068), 10)
	#df = data.frame (TOOL=tool, MODEL=model, SNPs=snps, P = pVal, CHR=chrom, POS=pos,pTHR=pThreshold)
	summTable = rbind (summTable, df)
	#summTableSorted = summTable %>% arrange (SNPs, P)
	summTableSorted = summTable %>% add_count (SNPs, sort=T) 
	write.table (file="summary-gwas.tbl", summTable, row.names=F,quote=F, sep="\t")
	write.table (file="summary-gwas-sorted.tbl", summTableSorted, row.names=F,quote=F, sep="\t")
}
