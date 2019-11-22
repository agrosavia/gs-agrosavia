#!/usr/bin/Rscript

# Get a table with SNP ids and Reference/Alternate alleles
# Remove duplicated SNPs

args = commandArgs (trailingOnly=T)

args = c ("potato_infinium_8303_map_context_DM_v3_superscaffolds.txt")

inFile = args [1]
name = strsplit (inFile, split="\\.")[[1]][1]
#outFile = "snps-ref-alt-nucleotides.tbl"

snpsAll   = read.table (inFile, header=T)
snps      = snpsAll [!duplicated (snpsAll$Name),]
ids       = snps [1]
seq       = snps [4]
alleles   = unlist(strsplit (substr (seq[,1], 52,54), split="/"))
allelesDF = data.frame (matrix (alleles, ncol=2, byrow=T))
refs      = as.character (allelesDF [,1])
alts      = as.character (allelesDF [,2])
AAs       = paste0 (refs, refs)
ABs       = paste0 (refs, alts)
BBs       = paste0 (alts, alts)


out = cbind (snp=ids, ref=refs, alt=alts, AA=AAs, AB=ABs, BB=BBs)
outFilename = sprintf ("%s-refs-alts-snps.tbl", name)
write.table (file=outFilename, out, quote=F, row.names=F, sep="\t")


