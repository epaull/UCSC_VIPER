#!/usr/bin/env	Rscript

# Get command line options:
library('getopt')

opt = getopt(matrix(c(
    'expression', 'e', 1, "character",
    'phenotypes', 'p', 1, "character",
    'output', 'o', 1, "character",
    'regulon', 'n', 1, "character",
    'test_phenotype', 't', 1, "character",
    'reference_phenotype', 'r', 1, "character",
    'num_results', 'y', 2, "double",
    'permutations','z',1000,"integer",
    'num_combin','c',2,"double",
    'min_size', 'm',2, "double",
    'viper_null', 'x', 0, "logical"
    ),ncol=4,byrow=TRUE));

## source the library, always in this relative location
calling_directory = dirname(get_Rscript_filename())
source(paste(calling_directory, "../lib", "viper-tools.R", sep="/"))

#options(error = quote({dump.frames(to.file = TRUE); q(status = 1)}))
#source("./shadow_combin.R")

##
## Read-in expression and phenotype data, perform basic error checking
##

exprs = parse.tab(opt$expression)
pheno = parse.phenotypes(opt$phenotypes)

# check phenotype data matches
#
if (length(setdiff(rownames(pData(pheno)),colnames(exprs)))!=0) {
	print ("Error: phenotype and expression sample lists must be identical!")
	q();
}

if ( is.null( opt$num_results) ) { opt$num_results<- 20 }
if ( is.null( opt$num_combin) ) { opt$num_combin<- 50 }
if ( is.null( opt$min_size) ) { opt$min_size<- 25 }
if ( is.null( opt$permutations) ) { opt$min_size<- 1000 }
expset.obj = makeExpressionSet(exprs, pheno)
print (expset.obj)
# parse the adj file to get a regulon object:
# note: all candidate regulators and genes must be in the 
# dataset you parsed
regulon <- opt$regulon 
regul <- NULL
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

##
## Do VIPER analysis and write out results
##

viper.result = run.viper.supervised(expset.obj, regul, opt$test_phenotype, opt$reference_phenotype, min.size=opt$min_size)
write.table(cbind(id=rownames(viper.result),viper.result), file=paste(opt$output, "/", "viperScores.txt", sep=""), row.names=FALSE,sep="\t", quote=F)

q();

