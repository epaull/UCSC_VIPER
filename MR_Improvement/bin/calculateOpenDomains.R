#!/usr/bin/env Rscript

# Get command line options:
library('getopt')

opt = getopt(matrix(c(
	'hicIntersect', 'i', 1, "character",
	'network', 'n', 1, "character",
	'output', 'o', 1, "character",
	'significance', 's', 1, "double"
	),ncol=4,byrow=TRUE));

## Test Options ##
#opt <- list()
#opt$hicIntersect <- 'dnaseCoverage.sorted.bed'
#opt$network <- 'brca_rnaseq851_multinet.adj'
#opt$output <- 'test.adj'
#opt$significance <- 0.05

# Input:
#	- A bam file representing post 'bamtools coverage' run against a Hi-C domain bedfile.
#		The first 3 columns are the chromosome the starting position and end position of
#		a Hi-C experiment domain. The last 4 columns represent the number of DNase regions
#		that overlapped the Hi-C domain interval, the number of bases in the Hi-C domain that
#		had non-zero coverage from intervals in a DNase bed file, the length of the Hi-C domain,
#		and the fraction of bases in the Hi-C domain that had non-zero coverage from DNase intervals.
#
# Returns: 
#	- A file formatted in the same way as the inputted Hi-C domain bedfile except containing only
#		the intervals of Hi-C domains that are determined to be open based off of the
#		ttest p-value of 0.05 on the fractional coverage data that is found on the 8th column of
#		the input file.

parse.tab <- function(tab_file){

  m <- read.table(tab_file, sep="\t", header=FALSE, quote="", check.names=FALSE)
  return (as.matrix(m))
}

dnaseCoverage <- parse.tab(opt$hicIntersect)
coverageFractionList <- as.double(dnaseCoverage[,8])

sigValue <- qnorm(opt$significance, mean=mean(coverageFractionList), sd=sd(coverageFractionList))

# Obtain gene list of genes that lie in closed Hi-C domains.
#	'closed' domains are defined as domains with a DNase coverage fraction that
#	is less than the cutoff value 'sigValue'.
closedGeneList <- c()
for (index in 1:length(coverageFractionList)){

	coverageFraction <- dnaseCoverage[index,8]

	# find and collect genes that lie in Hi-C 'closed domains' into a list
	if (coverageFraction < sigValue){
		closedGeneList <- c(closedGeneList, unlist(strsplit(dnaseCoverage[index,4], ":")))
	}
}


# parse out master regulators that are contained within the closed gene list then
#	parse out all targets of MRs that are contained within the closed gene list
#	of the '.adj' gene expression regulatory network file.
regNetwork <- file(opt$network, open = 'r')

finalRegNetwork <- list()
counter <- 1
while (length(oneLine <- readLines(regNetwork, n=1, warn = FALSE)) > 0){
	print(counter)
	targetList <- c()
	targetList <- unlist(strsplit(oneLine, "\t"))
	finalTargetList <- c()
	if (!(targetList[1] %in% closedGeneList)){
		index <- 2
		
		while (index < length(targetList)) {
			if (!(targetList[index] %in% closedGeneList)){
				finalTargetList <- c(finalTargetList, targetList[index:(index+1)])
				index <- index + 2
			} else{
				index <- index + 2
				next
			}
		}
	} else {
		next
	}
	finalRegNetwork[[counter]] <- c(targetList[1], finalTargetList)
	counter <- counter + 1
}

close(regNetwork)

lapply(finalRegNetwork, write, file=opt$output, append=T, ncolumns=5000, sep="\t" )


q();

#DNase_Coverage_Fraction = coverageFractionList
#
#library(mixtools)
#
## plot fitted normal distribution
#
#pdf(file='hist.pdf')
#
#h<-hist(DNase_Coverage_Fraction, breaks=100, col="red", main="Histogram of DNase Coverage", xlab="DNase Coverage", ylab="Frequency")
#xfit<-seq(min(DNase_Coverage_Fraction),max(DNase_Coverage_Fraction),length=length(DNase_Coverage_Fraction)) 
#yfit<-dnorm(xfit,mean=mean(DNase_Coverage_Fraction),sd=sd(DNase_Coverage_Fraction))
#yfit <- yfit*diff(h$mids[1:2])*length(DNase_Coverage_Fraction)
#lines(xfit, yfit, col="blue", lwd=2)
#
#dev.off()
#
#
#sigValue <- qnorm(0.05, mean=mean(DNase_Coverage_Fraction), sd=sd(DNase_Coverage_Fraction))
#abline(v=sigValue, col='green')
#
#
## plot fitted bimodal distribution with kernel density
#
#pdf(file='kd_Hist.pdf')
#
#bimodalDist <- normalmixEM(DNase_Coverage_Fraction, k=2, mu=c(0.01,0.19), sigma=c(0.005, 0.05), lambda=c(0.05,0.95)) 
#plot(bimodalDist, which=2)
#lines(density(DNase_Coverage_Fraction), lty=2, lwd=2)
#
#
#dev.off()
