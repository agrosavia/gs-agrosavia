library (GWASpoly)
data = read.GWASpoly (ploidy = 4, pheno.file = "phenotype-TableS2.tbl", 
                      geno.file = "genotype-TableS1.tbl", format = "ACGT", n.traits = 13, delim=",")

ph = read.table ("phenotype-TableS2.tbl", header=T, row.names = 1, sep=",")
gn = read.table ("genotype-TableS1.tbl", header=T, row.names = 1, sep=",")

# Populations structure by kinship
data2 <- set.K(data)

params <- set.params(fixed=c("Grp1","Grp2","Grp3","Grp4"),
                     fixed.type=rep("numeric",4))

# GWAS execution
data3 = GWASpoly(data2, models=c("general","additive","1-dom", 
                                 "2-dom"),traits=c("tuber_shape","tuber_eye_depth"),
                                 params=params)
# QQ-plot Output
par(mfrow=c(2,3)) #specifies a 2 x 3 panel
models <- c("additive","general","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")
for (i in 1:6) {
  qq.plot(data3,trait="tuber_shape",model=models[i])
}  

# QTL Detection
data4 = set.threshold (data3, method="Bonferroni",level=0.05)
get.QTL (data4)

# Manhattan plot Output
par(mfrow=c(1,2)) #specifies a 1 x 3 panel
manhattan.plot (data4, trait="tuber_shape", model="additive")
manhattan.plot (data4, trait="tuber_eye_depth", model="additive")
  



