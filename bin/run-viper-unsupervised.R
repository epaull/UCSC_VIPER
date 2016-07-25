#!/usr/bin/env	Rscript

# Get command line options:
library('getopt')

opt = getopt(matrix(c(
    'expression', 'e', 1, "character",
    'output', 'o', 1, "character",
    'regulon', 'n', 1, "character",
    'regulon_minsize', 'i', 2, "integer",
    'sample_signature_method', 's', 2, "character"
    ),ncol=4,byrow=TRUE));

## source the library, always in this relative location
calling_directory = dirname(get_Rscript_filename())
source(paste(calling_directory, "../lib", "viper-tools.R", sep="/"))


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


# Alana adding for paired ttest
sample_signature_method <- opt$sample_signature_method
if (is.null(sample_signature_method)) {
sample_signature_method <- "scale"
print(sample_signature_method)
}



##
## Just run unsupervised VIPER (i.e. median centered data)
##

print ("Running unsupervised VIPER inference...")
vpres <- run.viper.unsupervised(exprs, regul, regulon_minsize)
print ("Done!")
print ("Writing result..")
viper.result <- vpres #it was originally t(vpres). Alana un-transposed b/c TFs by samples is easier to work with downstream
write.table(cbind(id=rownames(viper.result),viper.result), file=paste(opt$output, "/", "viperScores.txt", sep=""), row.names=FALSE, sep="\t", quote=F)

# Alana added this so we save the .Rdata:
save.image(file=paste(opt$output, "/", "master-reg.RData", sep=""))

q();
