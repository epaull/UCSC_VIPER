#!/usr/bin/env	Rscript

options(warn = -1)
library('getopt')

opt = getopt(matrix(c(
    'tissue' , 't', 1, "character",
    'folder' , 'h', 1, "character"
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

folder <- opt$folder
print(opt$folder)
score_file = paste(opt$folder, "scores.txt", sep='/')
bg_file = paste(folder, "background_models", "swapped_events_scores_global.txt", sep='/')
bg_specific_file = paste(folder, "background_models", "swapped_events_scores_specific.txt", sep='/')

sample.scores = as.matrix(read.table(score_file, sep="\t", header=FALSE))
real.scores = as.numeric(sample.scores[,2])
permuted.data = as.matrix(read.table(bg_file, sep="\t", header=FALSE))
permuted.scores = as.numeric(permuted.data[,2])
pss = as.numeric(as.matrix(read.table(bg_specific_file, sep="\t", header=FALSE))[,2])

print("Real:")
summary(real.scores)
print("Permuted Global:")
summary(permuted.scores)

#ks.test(real.scores, permuted.scores)
wilcox.test(real.scores, permuted.scores)

if (is.null(opt$tissue)) {
	opt <- data.frame(tissue="tiedie_backgroundModel")
}

#for (drug in sample.scores[,1]) {
#	background <- permuted.scores[which(permuted.data[,1] == drug)]
#	real <- real.scores[which(sample.scores[,1] == drug)]
#	print (paste(drug, (real - mean(background)), sep=" "))
#}


png(paste(folder, paste(opt$tissue, ".png", sep=""), sep="/"))
plot(density(real.scores), col="black", lwd=2, main=paste(opt$tissue, " real (black) vs permuted (red) distributions"), xlab="Spearman Correlation to Centroid", xlim=c(0.8, 1))
lines(density(permuted.scores), col="red", lwd=2)

png(paste(folder, paste(opt$tissue, "_specific.png", sep=""), sep="/"))
plot(density(real.scores), col="black", lwd=2, main=paste(opt$tissue, " real (black) vs permuted (red) distributions"), xlab="Spearman Correlation to Centroid")
lines(density(pss), col="red", lwd=2)

# calculate p-vals
fdrs <- getFDRs(real.scores, permuted.scores)
all.scores <- cbind(sample.scores, fdrs)
write.table(all.scores, file=paste(folder, 'fdrs.txt', sep='/'), sep="\t", quote=FALSE)

q();
