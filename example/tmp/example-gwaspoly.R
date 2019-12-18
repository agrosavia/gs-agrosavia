#!/usr/bin/Rscript

#----------------------------------------------------------
# Execute GWASPoly with example pheno and geno files
#----------------------------------------------------------

library (GWASpoly)

args = commandArgs(trailingOnly = TRUE)
args = c ("phenotype-TableS2.tbl", "genotype-TableS1.tbl")
print (args)
phenotypeFile = args [1]
genotypeFile  = args [2]

data = read.GWASpoly (ploidy = 4, pheno.file = phenotypeFile, geno.file = genotypeFile, format = "ACGT", n.traits = 13, delim=",")

ph = read.table (phenotypeFile, header=T, row.names = 1, sep=",")
gn = read.table (genotypeFile, header=T, row.names = 1, sep=",")

# Populations structure by kinship
#data2 <- set.K(data)
data2  <- set.K(data,K=NULL)
params <- set.params(fixed=NULL, fixed.type=NULL) # Naive  
#params <- set.params(fixed=c("Grp1","Grp2","Grp3","Grp4"), fixed.type=rep("numeric",4))
# params <- set.params(n.PC=10, fixed=NULL, fixed.type=NULL) # PCs

# GWAS execution
data3 = GWASpoly(data2, models=c("general","additive","1-dom", "2-dom"),traits=c("tuber_shape","tuber_eye_depth"), params=params)

# QQ-plot Output
for (TRAIT in c("tuber_shape", "tuber_eye_depth")) {
	pdf (paste0("out-", TRAIT,".pdf"))
		message (">>>> QQ-plot...")
		par(mfrow=c(2,3)) #specifies a 2 x 3 panel
		models <- c("additive","general","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")
		for (i in 1:6) {
		  qq.plot(data3,trait=TRAIT,model=models[i])
		}  

		# QTL Detection
		message (">>>> QTL Detection ...")
		data4 = set.threshold (data3, method="Bonferroni",level=0.05)
		get.QTL (data4)

		message (">>>> Manhattan plot...")
		# Manhattan plot Output
		par(mfrow=c(2,3)) #specifies a 2 x 3 panel
		models <- c("additive","general","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")
		for (i in 1:6) {
			manhattan.plot (data4, trait=TRAIT, model=models[i])
			write.GWASpoly (data4, TRAIT, paste0("out-",TRAIT,".scores"), what="scores", "delim"="\t")
		}  
	dev.off()
}



