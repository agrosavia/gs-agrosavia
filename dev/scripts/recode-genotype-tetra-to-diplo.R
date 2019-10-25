#!/usr/bin/Rscript
# Recode a tetraploid genotipe (0,1,2,3) to diploid (0,1,2)

genotype = read.csv ("genotype-checked.tbl", header=T)
dossages = genotype [,-c(1,2,3)]

dossages [dossages==2] = 1
dossages [dossages==3] = 1
dossages [dossages==4] = 2

newGenotype = cbind (genotype [,c(1,2,3)], dossages)

write.csv (file="genotype-diplo.tbl", newGenotype, quote=F,row.names=F)

