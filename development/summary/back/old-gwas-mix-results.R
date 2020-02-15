#!/usr/bin/Rscript
options (stringsAsFactor=F, width=300)

library ("dplyr")

h <- function (x,n=5) {head (x)[,1:n]}

# Read data
pl = read.table ("out-plink.tbl", header=T)
ts = read.table ("out-tassel.tbl", header=T, sep="\t")
gd = read.table ("out-gwaspd.tbl", header=T, sep="\t")
gt = read.table ("out-gwaspt.tbl", header=T, sep="\t")

# Extract common columns
dgt = gt %>% select (Marker, Score, Chr) %>% mutate (TOOL="Gwasp2", MULTI="g2", N=1,pAVR=Score) %>% 
	rename (SNP=Marker, pVAL=Score, CHR=Chr)
dgd = gd %>% select (Marker, Score, Chr) %>% mutate (TOOL="Gwasp4", MULTI="g4", N=1,pAVR=Score) %>%
	rename (SNP=Marker, pVAL=Score, CHR=Chr)
dpl = pl %>% select (SNP, BONF, CHR) %>% mutate (TOOL="Plink", MULTI="pl", N=1,pAVR=BONF) %>%
	rename (SNP=SNP, pVAL=BONF, CHR=CHR)
dts = ts %>% select (Marker, p, Chr) %>% mutate (TOOL="Tassel", MULTI="ts", N=1,pAVR=p) %>%
	rename (SNP=Marker, pVAL=p, CHR=Chr)


dtt= rbind (dgt,dgd,dpl,dts)
print(head (dtt))
