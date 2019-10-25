#!/usr/bin/Rscript
# r2.0: Support for multiple command line arguments (default only)
# r1.5: Added SNPs annotations from Ship
# r1.4: Working with data from agrosavia, with improved graphics
# r1.2: With basic parameters handling and best graphics
# r0.1: Functions + write to RData cache

USAGE="USAGE: gwas-pipeline-poly.R <phenotype> <genotype> [structure]"

library (GWASpoly)
args = commandArgs(trailingOnly = TRUE)

setClass ("GWASpolyStruct", slots=c(params="list"), contains="GWASpoly.K")

options (width=300)

#gwasModelsTypes = c("Naive", "Kinship", "Structure", "PCs", "Kinship+Structure", "Kinship+PCs", "Structure+PCs", "Kinship+Structure+PCs")
snpModels       = c("general","additive","1-dom", "2-dom")
testModels      = c("general", "additive","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")
##gwasModelsTypes = c("Kinship")
##snpModels       = c("general")
##testModels      = c("general")
##multipleTestingModel = "FDR"
##multipleTestingModel = "permute"
multipleTestingModel = "Bonferroni" 
#snpModels            = c("general","additive","1-dom")
#testModels           = c("general", "additive","1-dom-alt","1-dom-ref")
#snpModels            = c("general")
#testModels           = c("general")
#-------------------------------------------------------------
# Global configs
#-------------------------------------------------------------
#genotypeFormati = "numeric"|"AB"|"ACGT"
genotypeFormat = "numeric"
#gwasModel = "Naive"|"Structure"|"Kinship"|"PCs"
#gwasModel      = "Kinship"
PLOIDY         = 4
#-------------------------------------------------------------
data = data1 = data2 = data3 = data4 = NULL
phenotype = genotype = structure = NULL

main <- function (args) {
	#args = c ("phenotype-gwaspoly-tuber_shape.tbl", "genotype-TableS1.tbl", "structure-gwaspoly.tbl")
	args = c ("--pheno", "phenotype-checked.tbl", "--geno", "genotype-checked.tbl", 
			  "--struct", "structure-checked.tbl", "--snps", "snps-annotations.tbl", "--model", "Naive")
	#args = c ("phenotype-checked.tbl", "genotype-checked-diplo.tbl", "structure-checked.tbl", "snps-annotations.tbl")

	# Read and check command line arguments
	params = readCheckCommandLineArguments (args)

	msg (">>>>>>>>>>>>", params$gwasModel, "<<<<<<<<<<<") 

	# Load cache data
	if (file.exists ("gwas.RData")) load (file="gwas.RData")

	# Read and check matchs in params (phenotype and genotype) and write new params
	data = readData (params$phenotypeFile, params$genotypeFile, params$structureFile)

	# Read input genotype and genotype (format: "numeric", "AB", or "ACGT")
	data1 <- initGWAS (params$phenotypeFile, params$genotypeFile, data1)

	# Set the kinship
	data2 <- setKinship (data1, params$gwasModel, data2)

	# Populations structure and kinship
	data3 <- setPopulationStructure (data2, params$gwasModel, data$phenotype, 
									data$structure, data3)

	# GWAS execution
	data4 <- runGwaspoly (data3, params$gwasModel, snpModels, data$traits, data4)

	# Plot results
	showResults (data4, testModels, data$traits, params$gwasModel, params$phenotypeFile, params$snpsFile)

	save(data, data1, data2, data3, data4, file="gwas.RData") 
}
#-------------------------------------------------------------
# Read and check command line arguments
# Currently: phenotye, genotype, structure
#-------------------------------------------------------------
readCheckCommandLineArguments <- function (args) {
	msg();msg("Reading Files...")

	i = match ("--pheno", args)
	phenotypeFile = args [i+1]

	i = match ("--geno", args)
	genotypeFile = args [i+1]

	i = match ("--struct", args)
	structureFile = args [i+1]

	i = match ("--snps", args)
	snpsFile = args [i+1]

	i = match ("--model", args)
	gwasModel = args [i+1]

	return (list(phenotypeFile=phenotypeFile, 
				 genotypeFile=genotypeFile,
				 structureFile=structureFile, 
				 snpsFile=snpsFile,
				 gwasModel=gwasModel))
}

#-------------------------------------------------------------
# Read the genotype and phenotype
#-------------------------------------------------------------
readData <- function (phenotypeFile, genotypeFile, structureFile) {
	msg();msg ("Reading data...");msg()

	if (!is.null (c(phenotype, genotype))) {
		msg (">>>> Loading phenotype and genotype"); 
	}else{
		msg (">>>> Reading phenotype and genotype"); 
		phenotype   <- read.table (phenotypeFile, header=T, sep=",")
		genotype    <- read.table (genotypeFile, header=T, sep=",")
		if (!is.null(structureFile))		
			structure   <- read.csv   (structureFile, header=T)
		else
			structure = NULL

		traits  <- c(colnames (phenotype)[2])
		print (traits)
	}

	return (list (phenotype=phenotype, genotype=genotype,structure=structure, traits=traits))
}

#-------------------------------------------------------------
#-------------------------------------------------------------
readCheckFiles <- function (phenotypeFile, genotypeFile, structureFile) {
	msg ("Reading files and checking colnames and rownames...")
	phenotypeAll = read.table (phenotypeFile, header=T, sep=",")
	genotypeAll  = read.table (genotypeFile, header=T, sep=",")
	structureAll = read.table (structureFile, header=T)

	samplesPheno  = phenotypeAll$sample
	samplesGeno   = colnames (genotypeAll)
	samplesStruct = structureAll$sample

	commonMarkers = Reduce (intersect, list (samplesGeno, samplesPheno, samplesStruct))

	phenotype  = phenotypeAll [phenotypeAll$sample %in% commonMarkers,]
	genoColums = c("Markers","Chrom","Position", commonMarkers)
	genotype   = genotypeAll  [,colnames(genotypeAll) %in% genoColums]
	structure  = structureAll [structureAll[,1] %in% commonMarkers,]

	write.table (file="phenotype-checked.tbl", phenotype, row.names=F, quote=F, sep=",")
	write.table (file="genotype-checked.tbl", genotype, row.names=F, quote=F, sep=",")
	write.table (file="structure-checked.tbl", structure, row.names=F, quote=F, sep=",")

	return (list (phenotype, genotype, structure))

}

#-------------------------------------------------------------
# Set the kinship matrix using the realized relationship matrix  
#-------------------------------------------------------------
setKinship <- function (data1,  gwasModel, data2) {
	msg();msg("Setting kinship...");msg()
 	# Load data instead calculate it
	if (!is.null (data2)) {msg ("Loading kinship..."); return (data2) }

	#kinshipMatrix = NULL
	if (gwasModel %in% c("Kinship", "Kinship|Structure")) {
		msg (">>>> With default kinship... ")
		kinshipMatrix = NULL
	}else { 
		msg (">>>> Without kinship...") 
		markerNames   = data1@pheno [,1]
		n             = length (markerNames)
		kinshipMatrix = matrix (diag (n), n, n, dimnames=list (markerNames, markerNames))
	}		
	data2  = set.K (data1, K=kinshipMatrix)
	return (data2)
}
		
#-------------------------------------------------------------
# Fix Populations structure and kinship
#-------------------------------------------------------------
setPopulationStructure <- function (data2, gwasModel, phenotype, structure, data3) {
	msg();msg ("Setting population structure...");msg()
	#setClass ("GWASpolyStruct", slots=c(params="list"), contains="GWASpoly.K")
	data3 = new ("GWASpolyStruct", data2)
	
	nPCs=0
	structNames = structTypes = NULL

	if (gwasModel %in% "PCs") {
		nPCs=20
		msg (">>>> With nPCS=", nPCs, "...")
	}else if (gwasModel %in% c("Structure", "Kinship|Structure")) {
		msg (">>>> With population structure...")
		structNames = colnames (structure[-1])
		structTypes = rep ("numeric", length (structNames))

		#data3@pheno = phenoStruct
		data3@pheno = phenotype
		data3@fixed = structure [,-1]
	}else 
		msg (">>>> Without populations structure")

	data3@params = set.params(n.PC=nPCs, fixed=structNames, fixed.type=structTypes)

	return (data3)
}
#-------------------------------------------------------------
# GWAS execution
#-------------------------------------------------------------
runGwaspoly <- function (data3, gwasModel, snpModels, traits, data4) {
	msg();msg("Running GWASpoly...");msg()

	if (!is.null (data4)) { msg (">>>> Loading GWASpoly..."); return (data4) }

 	if (gwasModel %in% c("Naive","Kinship")) {
		msg (">>>> Without params")
		data4 = GWASpoly(data3, models=snpModels, traits=traits, params=NULL)
	}else {
		msg (">>>> With params")
		data4 = GWASpoly(data3, models=snpModels, traits=traits, params=data3@params)
	}
	
	return (data4)
}

#-------------------------------------------------------------
# Plot results
#-------------------------------------------------------------
showResults <- function (data4, testModels, traits, gwasModel, phenotypeFile, snpsAnnotationsFile) {
	msg();msg ("Plotting results...");msg()
	phenoName = strsplit (phenotypeFile, split=".tbl")[[1]][1]
	plotName = sprintf("out-qqPlot-%s-%s.pdf", gwasModel, traits [1])

	# QTL Detection
	data5 = set.threshold (data4, method=multipleTestingModel,level=0.05,n.permute=100,n.core=3)
	#significativeQTLs = get.QTL (data5)
	significativeQTLs = getQTL (data5, snpsAnnotationsFile)

	msg (">>>>", "Writing results...")
	outFile = sprintf ("out-significativeQTLs-%s-%s.tbl", gwasModel, traits[1]) 
	write.table (file=outFile, significativeQTLs, quote=F, sep="\t", row.names=F)

	print (significativeQTLs)


	pdf (file=plotName, width=11, height=7)
	n = length (testModels)
	#op <- par(mfrow = c(2,n),   # 2xn layout
          #oma = c(0,3,0,2) + 0.1,   # zero rows of text at the outer left and bottom margin
          #mar = c(0,3,1.5,1) + 0.1, # space for three row of text at ticks and to separate plots
		  #mgp = c(1.8,0.6,0)+0.1,   # axis label at 1.8 rows distance, tick labels at 0.6 row
		  #	xpd = F,            # allow content to protrude into outer margin (and beyond)
		  #cex.main=1)
	
	#par(mfrow=c(2,6), mgp=c(2,0.7,0), oma=c(0,0,0,2), mar=c(0,4,1.5,0), xpd=F)
	
	# QQ-plot Output
	op <- par(mfrow = c(2,n), oma=c(0,0,3,0))
	for (i in 1:length(testModels)) {
		#par (cex.main=0.5, cex.lab=0.5, cex.axis=0.5, ann=T)
		qqPlot(data4,trait=traits [1], model=testModels[i], cex=0.3)
	}

	# Manhattan plot Output
	for (i in 1:length(testModels)) {
		#par (cex=1.5)
		manhattan.plot (y.max=20,data5, trait=traits[1], model=testModels [i])
	}
	plotTitle = sprintf ("GWAS %s-ploidy with %s for %s trait", PLOIDY, gwasModel, traits[1])  
	mtext(plotTitle, outer=T,  cex=1.5,  line=0)
	par(op)
	dev.off()
}

#-------------------------------------------------------------
# Extracts significant QTL
#-------------------------------------------------------------
getQTL <- function(data,snpsAnnotationsFile, traits=NULL,models=NULL) {
	stopifnot(inherits(data,"GWASpoly.thresh"))
	if (is.null(traits)) {
		traits <- names(data@scores)
	} else {
		stopifnot(is.element(traits,names(data@scores)))
	}
	if (is.null(models)) {
		models <- colnames(data@scores[[1]])
	} else {
		stopifnot(is.element(models,colnames(data@scores[[1]])))
	}

	snpsAnnotations <<- read.csv (file=snpsAnnotationsFile, header=T)

	n.model <- length(models)
	n.trait <- length(traits)
	output <- data.frame(NULL)
	for (i in 1:n.trait) {
		for (j in 1:n.model) {
			ix <- which(data@scores[[traits[i]]][,models[j]] > data@threshold[traits[i],models[j]])
			markers <-  data.frame (SNP=data@map[ix,c("Marker")])
			print (markers)
			snpAnn  <- merge (markers, snpsAnnotations, by.x="SNP",by.y="SNP_id", sort=F)[,c(2,7)]
			print (snpAnn)

			n.ix <- length(ix)
			output <- rbind(output,data.frame(Trait=rep(traits[i],n.ix),
											  Model=rep(models[j],n.ix),
											  Threshold=round(rep(data@threshold[traits[i],models[j]],n.ix),2),
											  data@map[ix,],
											  Score=round(data@scores[[traits[i]]][ix,models[j]],2),
											  Effect=round(data@effects[[traits[i]]][ix,models[j]],2),
											  snpAnn,
											  stringsAsFactors=F,check.names=F))
		}
	}
	return(output)
}
#-------------------------------------------------------------
# QQ plot
#-------------------------------------------------------------
qqPlot <- function(data,trait,model,cex=1,filename=NULL) {
	stopifnot(inherits(data,"GWASpoly.fitted"))
	traits <- names(data@scores)
	stopifnot(is.element(trait,traits))
	models <- colnames(data@scores[[trait]])
	stopifnot(is.element(model,models))
	scores <- data@scores[[trait]][,model]
	remove <- which(is.na(scores))
	if (length(remove)>0) {
		x <- sort(scores[-remove],decreasing=TRUE)
	} else {
		x <- sort(scores,decreasing=TRUE)
	}
	n <- length(x)
	ifactor = inflationFactor (x)
	unif.p <- -log10(ppoints(n))
	if (!is.null(filename)) {postscript(file=filename,horizontal=FALSE)}
	par(pty="s")
	plot(unif.p,x,pch=16,cex=cex,
		 xlab=expression(paste("Expected -log"[10],"(p)",sep="")),
		 ylab=expression(paste("Observed -log"[10],"(p)",sep="")),
		 main=paste(trait," (",model,") ",sep=""))

	mtext (bquote(lambda[GC] == .(ifactor)), side=3, line=-2, cex=0.7)

	lines(c(0,max(unif.p)),c(0,max(unif.p)),lty=2)
	if (!is.null(filename)) {dev.off()}
	return(NULL)
}

#-------------------------------------------------------------
# Calculate the inflation factor from -log10 values
#-------------------------------------------------------------
inflationFactor <- function (scores){
	pvalues = 10^-scores
	chisq <- na.omit (qchisq(1-pvalues,1))
	delta  = round (median(chisq)/qchisq(0.5,1), 3)

	return (delta)
}

#-------------------------------------------------------------
# Read input genotype and genotype (format: "numeric", "AB", or "ACGT")
#-------------------------------------------------------------
initGWAS <- function (phenotypeFile, genotypeFile, data1) {
	msg();msg ("Initializing GWAS...");msg()
	# When data is previously loaded
	if (!is.null (data)) {msg(">>>> Loading GWAS data..."); return (data)}

	data1 = read.GWASpoly (ploidy = PLOIDY, pheno.file = phenotypeFile, geno.file = genotypeFile, 
						  format = genotypeFormat, n.traits = 1, delim=",")

	return (data1)
}

#-------------------------------------------------------------
# Print a log message with the parameter
#-------------------------------------------------------------
msg <- function (...) {
  messages = unlist (list (...))
  cat (">>>>", messages, "\n")
}

#create a new binary pipe operator
`%notin%` <- function (x, table)
    is.na(match(x, table, nomatch = NA_integer_))

#-------------------------------------------------------------
# Define a set of initial models, traits, and gwas analyses
#-------------------------------------------------------------
old_setInitials <- function () {
	msg ("Setting initials...")
	#testModels  = snpModels = c ("general")
	snpModels = c("general","additive","1-dom", "2-dom")
	testModels <- c ("additive","general","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")
	#testModels  = snpModels = c ("additive", "general")
	testTraits <- c ("gota")
	msg (">>>> Models: ", testModels)
	msg (">>>> Traits: ", testTraits)

	return (list (testModels, snpModels, testTraits))
}


#-------------------------------------------------------------
# Define a set of initial models, traits, and gwas analyses
#-------------------------------------------------------------
old_setInitials <- function () {
	msg ("Setting initials...")
	#testModels  = snpModels = c ("general")
	snpModels = c("general","additive","1-dom", "2-dom")
	testModels <- c ("additive","general","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")
	#testModels  = snpModels = c ("additive", "general")
	testTraits <- c ("gota")
	msg (">>>> Models: ", testModels)
	msg (">>>> Traits: ", testTraits)

	return (list (testModels, snpModels, testTraits))
}


#-------------------------------------------------------------
# Call main 
#-------------------------------------------------------------
main (args)


