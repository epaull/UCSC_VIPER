#!/usr/bin/env	Rscript

# Get command line options:
library('getopt')

opt = getopt(matrix(c(
    'expression', 'e', 1, "character",
    'output', 'o', 1, "character",
    'regulon', 'n', 1, "character",
    'regulon_minsize', 'i', 2, "integer"
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

run.viper.unsupervised <- function (exp.obj, regulon, regul.minsize) {
	vpres <- viper(exp.obj, regulon, minsize=regul.minsize)
	return (vpres)
}


##
## Read-in expression and phenotype data
##

exprs = parse.tab(opt$expression)

# parse the adj file to get a regulon object:
# note: all candidate regulators and genes must be in the 
# dataset you parsed
regulon <- opt$regulon 
if (grepl('.rda', regulon)) {
	# can't determine the variable name before this step, but it should
	# be named 'regul'
	load(regulon)
} else if (grepl('.adj', regulon)) {
	regul <- aracne2regulon(regulon, exprs)
} else {
	print("Unrecognized regulon file type!")
	q();
}

regulon_minsize <- as.numeric(opt$regulon_minsize)
if (is.null(opt$regulon_minsize)) {
regulon_minsize <- 25
}


##
## Just run unsupervised VIPER (i.e. median centered data)
##

print ("Running unsupervised VIPER inference...")
vpres <- run.viper.unsupervised(exprs, regul, regulon_minsize)
print ("Done!")
print ("Writing result..")
viper.result <- t(vpres)
write.table(cbind(id=rownames(viper.result),viper.result), file=paste(opt$output, "/", "viperScores.txt", sep=""), row.names=FALSE, sep="\t", quote=F)



save.image(file=paste(opt$output, "/", "master-reg.RData", sep=""))
q();


