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
## Do MARINa Master Regulator Analysis, and write out results.
##
#mrs = run.marina(expset.obj, regul, opt$test_phenotype, opt$reference_phenotype, regul.minsize=opt$min_size, num.permutations=opt$permutations)
#mr.summary = summary(mrs,mrs=length(mrs$regulon))
#
#write.table(mr.summary, file=paste(opt$output, "/", "masterRegulators.txt", sep=""),row.names=FALSE,sep="\t", quote=F)
#
#pdf(file=paste(opt$output, "/", "masterRegulators.pdf", sep=""))
#plot(mr.summary, ,mrs=min(opt$num_results,length(mrs$regulon)),cex=0.7)
#dev.off()

#leading edge
#mr.le <- lapply(names(mrs$regulon),function(x) mrs$ledge[[x]])
#names(mr.le) <- names(mrs$regulon)
#writeSetList(mr.le,out.file=paste(opt$output, "/", "masterRegulatorsLedge.listt", sep=""))
#
##
## Do VIPER analysis and write out results
##

viper.result = run.viper.supervised(expset.obj, regul, opt$test_phenotype, opt$reference_phenotype, min.size=opt$min_size)
write.table(cbind(id=rownames(viper.result),viper.result), file=paste(opt$output, "/", "viperScores.txt", sep=""), row.names=FALSE,sep="\t", quote=F)
save.image(file=paste(opt$output, "/", "master-reg.RData", sep=""))

q();

##
## Everything else: combinatorial analysis, shadow, synnergy, etc. 
##

#combinatorial analysis
mrs.combin <- msviperCombinatorial(mr.summary,regulators = opt$num_combin,level=5,minsize=opt$min_size)
mrs.combin <- ledge(mrs.combin)

combin.summary=summary(mrs.combin,mrs=length(mrs.combin$regulon))
write.table(combin.summary, file=paste(opt$output, "/", "masterRegulatorsCombinatorial.txt", sep=""),row.names=FALSE, sep="\t", quote=F)

pdf(file=paste(opt$output, "/", "masterRegulatorsCombinatorial.pdf", sep=""))
plot(mrs.combin,mrs=min(opt$num_results,length(mrs.combin$regulon)),cex=0.7)
dev.off()

#combin regulons leading edge
#leading edge
mr.combin.le <- lapply(names(mrs.combin$regulon),function(x) mrs.combin$ledge[[x]])
names(mr.combin.le) <- names(mrs.combin$regulon)

writeSetList(mr.combin.le,out.file=paste(opt$output, "/", "masterRegulatorsCombinatorialLedge.listt", sep=""))

#synergy between MRs
#only compute if combinatorial analysis produces any synergy regulons
if(length(mrs.combin$regulon)>length(mr.summary$regulon)){
mrs.synergy <- msviperSynergy(mrs.combin)
synergy.summary <- summary(mrs.synergy,mrs=length(mrs.synergy$regulon))
write.table(synergy.summary, file=paste(opt$output, "/", "masterRegulatorsSynergy.txt", sep=""), row.names=FALSE, sep="\t", quote=F)

pdf(file=paste(opt$output, "/", "masterRegulatorsSynergy.pdf", sep=""))
plot(mrs.synergy,mrs=min(opt$num_results,length(mrs.synergy$regulon)),cex=0.7)
dev.off()
}

#shadow analysis
mrs.shadow.001 <- shadow(mr.summary,regulators=0.1,shadow=0.01)
mrs.shadow.001.shadow.tfs <- lapply(unique(mrs.shadow.001$shadow[,2]),function(x) mrs.shadow.001$shadow[mrs.shadow.001$shadow[,2]==x,1])
names(mrs.shadow.001.shadow.tfs) <- unique(mrs.shadow.001$shadow[,2])

if(length(mrs.shadow.001.shadow.tfs)>0){
mrs.shadow.001.shadow.tfs <- mrs.shadow.001.shadow.tfs[names(sort(sapply(mrs.shadow.001.shadow.tfs,length),decreasing = T))]
writeSetList(mrs.shadow.001.shadow.tfs,out.file=paste(opt$output, "/", "shadow_single_tfs_pval_01.listt", sep=""))
}

## mrs.combin.shadow.001 <- shadow.combin(mrs.combin,regulators=0.1,shadow=0.01)
## mrs.combin.shadow.001.shadow.tfs <- lapply(unique(mrs.combin.shadow.001$shadow[,2]),function(x) mrs.combin.shadow.001$shadow[mrs.combin.shadow.001$shadow[,2]==x,1])
## names(mrs.combin.shadow.001.shadow.tfs) <- unique(mrs.combin.shadow.001$shadow[,2])
## mrs.combin.shadow.001.shadow.tfs <- mrs.combin.shadow.001.shadow.tfs[names(sort(sapply(mrs.combin.shadow.001.shadow.tfs,length),decreasing = T))]
## writeSetList(mrs.combin.shadow.001.shadow.tfs,out.file=paste(opt$output, "/", "shadow_combin_tfs_pval_01.listt", sep=""))

#store the run
save.image(file=paste(opt$output, "/", "master-reg.RData", sep=""))
q();


