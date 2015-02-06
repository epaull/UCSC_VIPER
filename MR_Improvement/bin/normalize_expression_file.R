#!/usr/bin/env Rscript

# Get command line options:
library('getopt')

opt = getopt(matrix(c(
	'exprData', 'e', 1, "character",
	'phenotypes', 'p', 1, "character",
	'normalization', 'n', 1, "character",	# either ('normal' or 'variance')
	'output', 'o', 1, "character"
	),ncol=4,byrow=TRUE));

# Input:
#	- An expression dataset 'e' (NOTE: requires that expression data be in RAW count format. NO PRIOR NORMALIZATION
#		SHOULD BE MADE PRIOR TO RUNNING THIS SCRIPT).
#	- A phenotype annotation file 'p'.
#	- A normalization strategy 'n'
#		If normalization is set using the parameter 'normal' the program will normalize the count data via
#			DESeq normalization.
#		If normalization is set using the parameter 'variance' the program will normalize the count data via
#			DESeq normalization along with variance stabilizations.
#	- The name of the output file for which to write the normalized expression data to 'o'.
#
# Returns: 
#	- An expression dataset file with its expression values normalized via which DESeq normalization strategy
#		was specified.

library("DESeq")

readCountTable = read.table(opt$exprData, header=TRUE, row.names=1, check.names=FALSE)
phenotypeTable = read.table(opt$phenotypes, header=TRUE, row.names=1, check.names=FALSE)

exprDesign = data.frame(
	row.names = colnames( readCountTable ),
	condition = as.vector(phenotypeTable[,1]),
	libType = rep('paired-end', ncol(readCountTable)) )

condition = factor(as.vector(phenotypeTable[,1]))

cds = newCountDataSet( readCountTable, condition )
cds = estimateSizeFactors( cds )

if (opt$normalization == 'normal'){
	outTable = counts(cds, normalized=TRUE)
} else if (opt$normalization == 'variance'){
	cdsBlind = estimateDispersions( cds, method="blind" )
	outTable = getVarianceStabilizedData(cdsBlind)
}
write.table(outTable, file=opt$output, col.names = NA, sep="\t", quote=F)

q()