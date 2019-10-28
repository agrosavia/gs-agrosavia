#!/usr/bin/Rscript
library (parallel)

models = c("Naive", "Structure", "Kinship", "Kinship+Structure")

cmmList = list()
for (mdl in models) {
	cmm = sprintf ("gwas-polypipeline.R --pheno phenotype.tbl --geno genotype-diplo.tbl --struct structure.tbl --snps snps-annotations.tbl --ploidy 2 --model %s", mdl)
	message (cmm)			   
	cmmList = append (cmmList, cmm)
}
#mclapply (cmmList, system, mc.cores=3)
sapply (cmmList, system)
