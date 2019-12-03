#!/usr/bin/Rscript

## Create gwaspoly files (geno and pheno) from plink files (.ped,.map, pheno)

library (dplyr)

options (width=300, stringsAsFactors=T)
args = commandArgs(trailingOnly=T)

#----------------------------------------------------------
# Plink separated alleles (e.g. A C A G A A) are joined (AC AG AA)
#----------------------------------------------------------
joinAlleles <- function (alleles) {
	nCols = ncol (alleles)
	mat = matrix (nrow = nrow (alleles), ncol=nCols/2)
	j = 0
	for (i in seq(1,nCols,2)){
		j=j+1
 		mat [,j] = cbind (paste0 (alleles [,i], alleles [,i+1]))
	}
	return (mat)
}

#----------------------------------------------------------
# Create gwaspoly genotype from plink  .ped and .map files
#----------------------------------------------------------
createGwaspolyGenotype <- function (plinkFile) {
	pedFile   = paste0(plinkFile,".ped")
	mapFile   = paste0(plinkFile,".map")

	outName = paste0 (plinkFile, "-FLT-gwasp.tbl")

	# Process genotype
	ped = read.table (pedFile, header=F,stringsAsFactors=F, colClasses=c("character"))
	samples = ped [,2]
	alleles = ped [,-c(1,2,3,4,5,6)]
	allelesJoined = t(joinAlleles (alleles))

	map     = read.table (mapFile, header=F)
	chrs    = map [,1]
	markers = as.character (map [,2])
	pos     = map [,4]

	message (">>> PED files:")
	print (ped[1:10,1:20])

	message ("\n>>> MAP file:")
	print (map[1:10,])

	#rownames (allelesJoined) = markers
	colnames (allelesJoined) = samples

	genotype = cbind (Markers=markers,Chrom=chrs,Position=pos,allelesJoined)
	gn = genotype 
	#allelesDF = cbind (map [,c(1,2,3,4)], allelesJoined)

	message ("\n>>> Gwaspoly genotype file:")
	print (genotype [1:10,1:10])
	write.table (file=outName, genotype, row.names=F,quote=F, sep=",")
	return (list(samples=as.numeric(samples), markers=markers))
}

#----------------------------------------------------------
#----------------------------------------------------------
createGwaspolyPhenotype <- function (phenoFile, sampleNames) {
	phenoGwaspFile = gsub (".tbl", "-FLT-gwasp.tbl", phenoFile)

	pheno      <<- read.table (phenoFile, header=T)
	#phTmp   <<-  pheno %>% filter (sampleNames %in% IID)
	#phenoDF = pheno %>% arrange (sampleNames %in% IID) %>% as.data.frame()
	#pheno %>% filter (sampleNames %in% IID)

	phenoDF = pheno %>% filter (IID %in% sampleNames) %>% arrange (match (IID, sampleNames)) %>% data.frame()

	phenoGwasp <<- cbind (Samples=phenoDF[,2],BLIGHT=phenoDF[,3])

	message ("\n>>> Gwaspoly penotype file:")
	print (head (phenoGwasp))
	write.table (file=phenoGwaspFile, phenoGwasp, row.names=F,quote=F,sep=",")

}

#----------------------------------------------------------
# Filter the tetra gwaspoly genotype using the plink filtered genotype
#----------------------------------------------------------
filterTetraGwaspolyWithPlinkGeno <- function (gwaspGenoTetraFile, samples, markers) {
	gwp = read.csv (file=gwaspGenoTetraFile, header=T, check.names=F)
	gwp [,2] = gwp [,2]+1

	gwpDF = gwp %>% select (c(1,2,3), as.character(samples)) %>% filter (Markers %in% markers) %>% 
			arrange (match(Markers, markers)) %>% data.frame(check.names=F)

	outName = gsub (".tbl", "-FLT.tbl", gwaspGenoTetraFile)
	write.csv (file=outName, gwpDF, quote=F,row.names=F)
}

#----------------------------------------------------------
# Main
#----------------------------------------------------------
args = c("geno-6", "pheno.tbl", "agrosavia-genotype-tetra-NUM.tbl", "agrosavia-structure.tbl")

plinkFile  = args [1]
phenoFile  = args [2]
gwaspGenoTetraFile = args [3]
structFile = args [4]


genotype = createGwaspolyGenotype (plinkFile)
samples = genotype$samples
markers = genotype$markers

createGwaspolyPhenotype (phenoFile, genotype$samples)

filterTetraGwaspolyWithPlinkGeno (gwaspGenoTetraFile, samples, markers)

struct = read.csv (file=structFile, header=T, check.names=F)

structDF = struct %>% filter (Samples %in% samples) %>% arrange (match (Samples, samples)) %>% data.frame (check.names=F)
outName = gsub (".tbl", "-FLT.tbl", structFile)
write.csv (file=outName, structDF, quote=F,row.names=F)



