#!/usr/bin/env	Rscript

options(warn = -1)
library('getopt')

library(ggplot2)

opt = getopt(matrix(c(
    'scores' , 's', 1, "character",
    'output' , 'o', 1, "character"
	),ncol=4,byrow=TRUE));

getDensity <- function(score, fit) {

	epsilon = 0.0001
	val <- fit(score)
	if (is.na(val)) {
		val <- epsilon
	}

	return (val)
}

getFDRs <- function(scores, background) {

	# modify the density kernel function here
	#fit <- density(scores)
	#approx.fit <- approxfun(fit)

	#bg.fit <- density(background)
	#approx.bg.fit <- approxfun(bg.fit)

	#pvals <- sapply(scores, function(x) {1 - (approx.fit(x) - getDensity(x, approx.bg.fit))})
	#fdrs <- sapply(scores, function(x) {1 - (approx.fit(x) - getDensity(x, approx.bg.fit))})

	pvals <- sapply(scores, function(x) { (length(which(background > x))+1)/(length(background)+1)})
	fdrs <- p.adjust(pvals, method="BH")

	return(fdrs)
}

#folder <- opt$folder
folder <- '.'
output_folder <- paste(folder, "analysis", sep='/')
if (!file.exists(output_folder)) {
	dir.create(file.path(output_folder))
}
report_file <- paste(output_folder, opt$output, sep='/')

score_file = opt$scores
scores = read.delim(score_file, header=F, sep="\t")

test_vals = list()

for ( subtype in unique(scores[,1]) ) {
	print(subtype)
	subtype.idx <- which(scores[,1] == subtype)
	ref.idx <- which(scores[2] == 'reference')

	ref.index <- intersect(ref.idx, subtype.idx)	
	ref.vals <- as.numeric(scores[ref.index, 4])
	png(paste(output_folder, '/', subtype, '.png', sep=""))
	plot(density(ref.vals), col="black", xlim=c(0.0, 1.0))

	test.idx <- which(scores[2] == 'test')
	test.index <- intersect(test.idx, subtype.idx)	
	test.vals <- as.numeric(scores[test.index, 4])
	for ( val in test.vals ) {
		abline(v=val, col="red")
	}

	print (summary(test.vals))
	test_vals[[subtype]] = test.vals
}

print (ks.test(test_vals[['MEK']], test_vals[['MTOR']]))

