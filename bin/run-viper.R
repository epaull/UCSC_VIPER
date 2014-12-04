#!/usr/bin/env	Rscript

# Get command line options:
library('getopt')

# the real network score --
opt = getopt(matrix(c(
    'expression', 'e', 1, "character",
    'phenotypes', 'p', 1, "character",
    'output', 'o', 1, "character",
    'regulon', 'n', 1, "character",
    'test_phenotype', 't', 1, "character",
    'reference_phenotype', 'r', 1, "character",
    'num_results', 'y', 1, "character",
    'viper_null', 'x', 0, "logical"
    ),ncol=4,byrow=TRUE));

library(mixtools)
#library(bcellViper)
library(viper)
library(Biobase)

# Input: 
#   - A TAB separated file: the first line is column header (sample) info. 
#       The first column has feature names 
#
# Returns:
#
#       - A data matrix with the 'colnames' and 'rownames' attributes set, the first column
#       a 'character' vector with the description/HUGO id's, and the following columns are
#       numeric data. 
parse.tab <- function (file) {

  # skip the first two lines (GCT metadata)
  # The first two columns are the IDS and names, respectively. Note, the gene names 
  # might not be unique in the training data...
  m <- read.table(file, sep="\t", row.names=1, header=TRUE, quote="", check.names=FALSE)
  return (as.matrix(m))
}

parse.phenotypes <- function (file) {
	m <- read.table(file, sep="\t", row.names=1, header=TRUE, quote="", check.names=FALSE)
	# hack to fix sample IDs and presevere the 'names' attribute. This needs to be fixed on
	# their side, most likely. 
	sampleID = m$sampleID
	attr(sampleID, 'names') <- rownames(m)
	Tissue = m$Tissue
	attr(Tissue, 'names') <- rownames(m)
	description = m$description
	attr(description, 'names') <- rownames(m)
	m = data.frame(I(sampleID), I(Tissue), I(description))
	metadata <- data.frame(labelDescription = c("TCGA Sample ID", "Tissue Type", "Sample Type"), row.names = c("sampleID", "Tissue", "description"))
	phenoData <- new("AnnotatedDataFrame", data = m, varMetadata = metadata)
	return (phenoData)
}

#
# Input:  
#       - A data matrix with columns as samples, rows are genes
#       - phenotype data should include sample annotations (tumor/normal classes)
#
# Returns:
#
#       A data frame of predictors, filtered down the the top X
#
makeExpressionSet <- function(data.matrix, phenoData) {

  exprSet <- new("ExpressionSet", exprs = data.matrix, phenoData = phenoData, 
		 annotation = "viper_input")

	return (exprSet)
}

run.viper.MR <- function (exp.obj, regulon, set1.label, set2.label) {

	# get set 1 indexes
	set1.idx <- which(colnames(exp.obj) %in% names(c(phenoData(exp.obj)$sampleID[which(phenoData(exp.obj)$description == set1.label)])) )
	set2.idx <- which(colnames(exp.obj) %in% names(c(phenoData(exp.obj)$sampleID[which(phenoData(exp.obj)$description == set2.label)])) )

	#set1.idx <- c(phenoData(exp.obj)$sampleID[which(phenoData(exp.obj)$description == set1.label)])
	#set2.idx <- c(phenoData(exp.obj)$sampleID[which(phenoData(exp.obj)$description == set2.label)])

	print (paste("Comparing ", length(set1.idx), " test samples with ", length(set2.idx), " reference samples..." ))

	# get just the data matrix from this object
	data.matrix = exprs(exp.obj)
	# generate the signature and normalize to z-scores
	signature <- rowTtest( data.matrix[,set1.idx], data.matrix[,set2.idx] )
	signature <- (qnorm(signature$p.value/2, lower.tail=F) * sign(signature$statistic))[, 1]
	# create the null model signatures
	nullmodel <- ttestNull(data.matrix[,set1.idx], data.matrix[,set2.idx], per=100, repos=T)
	# compute MARINa scores based on the null model
	mrs <- msviper(signature, regulon, nullmodel)
	print (summary(mrs))
	png(file=paste(opt$output, "/", "masterRegulators.png", sep=""))
	plot(mrs, cex=0.7)
	# create a background viper signature based on relative levels, then compute the final scores
	vpres <- NULL
	if (!is.null(opt$viper_null)) {
		print ("Constructing Viper Signature")
		vpsig <- viperSignature(data.matrix[,-set2.idx], data.matrix[,set2.idx], method="ttest", verbose=T)
		print ("Constructing Viper Inferences")
		vpres <- viper(vpsig, regulon)
	} else {
		print ("Constructing Viper Inferences")
		vpres <- viper(exp.obj, regulon)
	}

	return (list(mrs, vpres))
}

##
## Read-in expression and phenotype data
##

exprs = parse.tab(opt$expression)
pdata = parse.phenotypes(opt$phenotypes)

# check phenotype data matches
#
if (!all(rownames(pdata) == colnames(exprs))) {
	print ("Error: phenotype and expression sample lists must be identical!")
	q();
}
#
expset.obj = makeExpressionSet(exprs, pdata)
print (expset.obj)
# parse the adj file to get a regulon object:
# note: all candidate regulators and genes must be in the 
# dataset you parsed
regulon <- opt$regulon 
regul <- aracne2regulon(regulon, expset.obj)


##
## The vignette example 
##
#data(bcellViper)
#naiveBcell <- which(colnames(dset) %in% names(c(phenoData(dset)$sampleID[which(phenoData(dset)$description == "N")])) )
#GCBcell <- which(colnames(dset) %in% c(names(c(phenoData(dset)$sampleID[which(phenoData(dset)$description == "CB")])), names(c(phenoData(dset)$sampleID[which(phenoData(dset)$description == "CC")]))))

#adjfile <- file.path(find.package("bcellViper"), "aracne", "bcellaracne.adj")
#regul <- aracne2regulon(adjfile, dset)
#
#result = run.viper.MR(dset, regul, "CB", "N")
#print (result)


num_results <- as.numeric(opt$num_results)
if (is.null(opt$num_results)) {
num_results <- 100
}

result = run.viper.MR(expset.obj, regul, opt$test_phenotype, opt$reference_phenotype)
mr.result = summary(result[[1]], mrs=num_results)
viper.result = result[[2]]

options(digits=2)
write.table(mr.result, file=paste(opt$output, "/", "masterRegulators.txt", sep=""), col.names = NA, sep="\t", quote=F)
write.table(format(viper.result, digits=2), file=paste(opt$output, "/", "viperScores.txt", sep=""), col.names = NA, sep="\t", quote=F)

png(file=paste(opt$output, "/", "masterRegulators.png", sep=""))
plot(result[[1]], cex=0.7)

save(mr.result, file=paste(opt$output, "/", "master-reg.RData", sep=""))
q();

#d1 = exprs(dset)
#signature <- rowTtest( d1[,GCBcell], d1[,naiveBcell] )
#signature <- (qnorm(signature$p.value/2, lower.tail=F) * sign(signature$statistic))[, 1]
##nullmodel <- ttestNull(d1[, GCBcell], d1[, naiveBcell], per=1000, repos=T)
#mrs <- msviper(signature, regulon, nullmodel)
#vpres <- viper(d1[,-naiveBcell], d1[,naiveBcell])
#
#summary(mrs)
#plot(mrs, cex=.7)

