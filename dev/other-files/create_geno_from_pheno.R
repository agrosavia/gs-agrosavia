#!/usr/bin/Rscript
# Create genotype by checking SNPs in the phenotype file
outGenotypeGwaspoly = "agrosavia-genotype-gwaspoly-checked.tbl"
outInvalidSNPs      = "agrosavia-genotype-gwaspoly-checked.errors"

ph = read.table ("agrosavia-phenotype-gwaspoly.tbl", header=T, sep=",")
snps = as.vector (ph [,1])

gn = read.table ("agrosavia-genotype-gwaspoly.tbl", header=T, sep=",")
gnSnps = colnames (gn)[-(1:3)] 

gnTbl = data.frame (gn [,1:3])
errorsTbl = c ("Invalid SNPs: in Phenotype but not in Genotype")
for (sn in snps) {
  if (sn %in% gnSnps) {
    gnCol = data.frame (gn [, sn])
    colnames (gnCol) = sn
    gnTbl = cbind (gnTbl, gnCol)
  }else {
    print (sprintf (">>> Invalid Snp %s ", sn))
    errorsTbl = append (errorsTbl, sn)
  }
  
}

write.table (file=outGenotypeGwaspoly, gnTbl, quote=F, sep=",", row.names=F)
write.table (file=outInvalidSNPs, errorsTbl, col.names=F,row.names=F)
