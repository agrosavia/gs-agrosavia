#!/usr/bin/Rscript
library (GWASpoly)

#-------------------------------------------------------------
# Print a log message with the parameter
#-------------------------------------------------------------
msg <- function (...) {
  messages = unlist (list (...))
  cat (">>>>", messages, "\n")
}
#-------------------------------------------------------------
  
args = commandArgs(trailingOnly = TRUE)
print (args)

phenotypeFile = args [1]
genotypeFile  = args [2]

#phenotypeFile = "agrosavia-phenotype-gwaspoly.tbl"
#genotypeFile  = "agrosavia-genotype-gwaspoly-checked.tbl"

#phex = read.table ("example/phenotype-TableS2.tbl", header=T, sep=",")
#gnex = read.table ("example/genotype-TableS1.tbl", header=T, sep=",")

#AnnShip = read.table ("other-files/potato_infinium_8303_map_context_DM_v3_superscaffolds.txt", header=T, sep="\t")
#AnnGeno = read.table ("other-files/potato_8303SNPs_potato_dm_v4.03.gff3", header=F, sep="\t")

ph = read.table (phenotypeFile, header=T, sep=",")
gn = read.table (genotypeFile, header=T, sep=",")


# Read input genotype and genotype (format: "numeric" or "ACGT")
message (">>> Reading data...")
data = read.GWASpoly (ploidy = 4, delim=",", format = "ACGT", n.traits = 1, 
                      pheno.file = phenotypeFile, geno.file = genotypeFile)

message (">>> Calculating kinship...")
# Populations structure by kinship
data2 <- set.K(data)

# Used to include population structure covariates
#params <- set.params(n.PC=10)
#params <- set.params(fixed=c("Grp1","Grp2","Grp3","Grp4"),
#                     fixed.type=rep("numeric",4))
params = NULL
# GWAS execution
message (">>> Running GWASpoly...")
data3 = GWASpoly(data2, models=c("general","additive","1-dom", 
                                 "2-dom"),traits=c("gota"))
# QQ-plot Output
message ("Ploting results...")
pdf (file="plots-qq-gwaspoly.pdf")
  par(mfrow=c(2,3)) #specifies a 2 x 3 panel
  models <- c("additive","general","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")
  for (i in 1:6) {
    qq.plot(data3,trait="gota",model=models[i])
  }
dev.off()

# QTL Detection
data4 = set.threshold (data3, method="Bonferroni",level=0.05)
get.QTL (data4)

# Manhattan plot Output
pdf (file="plots-manhattan-gwaspoly.pdf")
  par(mfrow=c(2,3)) #specifies a 1 x 3 panel
  for (i in 1:6) {
    manhattan.plot (data4, trait="gota", model=models [i])
  }  
dev.off ()


