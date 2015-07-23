#!/usr/bin/env Rscript

options(warn = -1)
library('getopt')

opt = getopt(matrix(c(
    'regulator' , 'r', 1, "character",
    'class_1' , '1', 1, "character",
    'class_2' , '2', 1, "character"
	),ncol=4,byrow=TRUE));

if (!is.null(opt$help) || is.null(opt$regulator) || is.null(opt$class_1)) {
	self = commandArgs()[1];
	#print a friendly message and exit with a non-zero error code
	cat(paste("Usage: ",self,"  --class_1 <viper scores> --class_2 <viper scores> --regulator <gene name> \n"))
	q();
}

# get all scores for this regulator
m_1 = read.delim(opt$class_1, header=T, row.names=1, as.is=T, sep='\t')
m_2 = read.delim(opt$class_2, header=T, row.names=1, as.is=T, sep='\t')
m1_scores = as.numeric(m_1[which(rownames(m_1) == opt$regulator),])
m2_scores = as.numeric(m_2[which(rownames(m_2) == opt$regulator),])

# plot distributions of two classes
png(paste(opt$regulator, '.png', sep=''))
plot(density(m1_scores), col="red", lwd=3, xlim=c(-10,10))
# blue == MAPK
lines(density(m2_scores), col="blue", lwd=3)
q();
