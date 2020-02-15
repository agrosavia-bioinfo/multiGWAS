kt = read.table ("tmp-plinkg.king.cutoff.out.id")
pl = read.csv (file="filtered-plink-genotype.ped", sep="\t", header=F)
library (dplyr)
markers=kt[,1]
plf = pl %>% filter (V2 %in% markers)
pl[1:10,1:100]
plf = pl %>% filter (V2 %in% markers)
write.table (file="out-king.tbl", plf, quote=F, sep="\t", row.names=F, col.names=F)
