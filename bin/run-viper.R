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
    'num_results', 'y', 2, "integer",
    'regulon_minsize', 'i', 2, "integer",
    'permutations', 'j', 2, "integer",
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

  m <- read.table(file, sep="\t", row.names=1, header=TRUE, quote="", check.names=FALSE)
  return (as.matrix(m))
}

parse.phenotypes <- function (file) {
	m <- read.table(file, sep="\t", row.names=1, header=TRUE, quote="", check.names=FALSE)
	phenoData <- new("AnnotatedDataFrame", data = m)
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

run.viper.MR <- function (exp.obj, regulon, set1.label, set2.label, max.results, regul.minsize, num.permutations) {

	# get set 1 indexes
	set1.idx <- which(colnames(exp.obj) %in% rownames(pData(exp.obj))[which(pData(exp.obj)[,1] == set1.label)])
	set2.idx <- which(colnames(exp.obj) %in% rownames(pData(exp.obj))[which(pData(exp.obj)[,1] == set2.label)])

	print (paste("Comparing ", length(set1.idx), " test samples with ", length(set2.idx), " reference samples..." ))

	# get just the data matrix from this object
	data.matrix = exprs(exp.obj)
	# generate the signature and normalize to z-scores
	signature <- rowTtest( data.matrix[,set1.idx], data.matrix[,set2.idx] )
	signature <- (qnorm(signature$p.value/2, lower.tail=F) * sign(signature$statistic))[, 1]
	# create the null model signatures
	nullmodel <- ttestNull(data.matrix[,set1.idx], data.matrix[,set2.idx], per=num.permutations, repos=T)
	# compute MARINa scores based on the null model
	mrs <- msviper(signature, regulon, nullmodel, minsize=regul.minsize)
	mr.summary <- summary(mrs, max.results)
	write.table(mr.summary, file=paste(opt$output, "/", "masterRegulators.txt", sep=""), col.names = NA, sep="\t", quote=F)
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
pheno = parse.phenotypes(opt$phenotypes)

# check phenotype data matches
#
if (length(setdiff(rownames(pData(pheno)),colnames(exprs)))!=0) {
	print ("Error: phenotype and expression sample lists must be identical!")
	q();
}

expset.obj = makeExpressionSet(exprs, pheno)
print (expset.obj)
# parse the adj file to get a regulon object:
# note: all candidate regulators and genes must be in the 
# dataset you parsed
regulon <- opt$regulon 
if (grepl('.rda', regulon)) {
	# can't determine the variable name before this step, but it should
	# be named 'regul'
	load(regulon)
} else if (grepl('.adj', regulon)) {
	regul <- aracne2regulon(regulon, expset.obj)
} else {
	print("Unrecognized regulon file type!")
	q();
}

# display the top X master regulators
num_results <- as.numeric(opt$num_results)
if (is.null(opt$num_results)) {
num_results <- 25
}
regulon_minsize <- as.numeric(opt$regulon_minsize)
if (is.null(opt$regulon_minsize)) {
regulon_minsize <- 25
}
num_permutations <- as.numeric(opt$permutations)
if (is.null(opt$permutations)) {
num_permutations <- 1000
}

result = run.viper.MR(expset.obj, regul, opt$test_phenotype, opt$reference_phenotype, num_results, regulon_minsize, num_permutations)
mr.result = result[[1]]
viper.result = result[[2]]

mr.summary = summary(mr.result, num_results)
write.table(mr.summary, file=paste(opt$output, "/", "masterRegulators.txt", sep=""), col.names = NA, sep="\t", quote=F)
write.table(viper.result, file=paste(opt$output, "/", "viperScores.txt", sep=""), col.names = NA, sep="\t", quote=F)

pdf(file=paste(opt$output, "/", "masterRegulators.pdf", sep=""))
plot(mr.result,num_results,cex=0.7)
dev.off()

save.image(file=paste(opt$output, "/", "master-reg.RData", sep=""))
q();


