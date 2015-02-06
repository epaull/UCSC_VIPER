#!/usr/bin/env Rscript

library('getopt')

opt = getopt(matrix(c(
	'input', 'i', 1, "character",
	'output', 'o', 1, "character"
	),ncol=4,byrow=TRUE));

library(org.Hs.eg.db)

# Input:
#	- An R object file containing a 'regul' object where the object contains entrez gene ids
#		that designate transcription factors as an attribute. Each of these has an attribute
#		'tfmode' and 'likelihood'. 'tfmode' lists all of the targets of the transcription factor
#		with its corresponding mutual information (MI) value. 'likelihood' lists all of the
#		targets of the transcription factor with its corresponding MI value likelihood.
#
# Returns:
#	- A file formatted in a way that can be converted into a '.adj' file format. It outputs
#		a TAB separated file: The first line is a header listed as 'Transcriptional Regulator'
#		'Predicted Target' and 'Mutual Information'. The following lines list 
#		the Transcriptional regulator in gene symbol format, the target the TF
#		regulates also in gene symbol format and then the mutual information score attached
#		to the relationship. 

# Set up entrez gene id to gene symbol converting object as a dictionary 'xx'
x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])

load(opt$input)

printStr = paste('Transcriptional Regulator', 'Predicted Target', 'Mutual Information', sep='	')
cat(printStr, file=opt$output, sep="\n")

for (entrezID in names(regul)) {

	regulator <- xx[[entrezID]]

	for (targetEntrezID in names(regul[[entrezID]]$tfmode)){
		target <- xx[[targetEntrezID]]
		MI <- regul[[entrezID]]$tfmode[[targetEntrezID]]

		printStr = paste(regulator, target, toString(MI), sep = '	')

		cat(printStr, file=opt$output, sep="\n")
	}
}