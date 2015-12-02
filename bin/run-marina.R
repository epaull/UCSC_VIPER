#!/usr/bin/env	Rscript

# Get command line options:
library('getopt')

# the real network score --
opt = getopt(matrix(c(
    'expression', 'e', 1, "character",
    'phenotypes', 'p', 1, "character",
    'output', 'o', 1, "character",
    'regulon', 'n', 1, "character",
    'test_phenotype', 't', 1, "character",
    'reference_phenotype', 'r', 1, "character",
    'num_results', 'y', 2, "integer",
    'regulon_minsize', 'i', 2, "integer",
    'num_combin','c',2,"double",
    'permutations', 'j', 2, "integer"
    ),ncol=4, byrow=TRUE));

## source the library, always in this relative location
calling_directory = dirname(get_Rscript_filename())
source(paste(calling_directory, "../lib", "viper-tools.R", sep="/"))

##
## Read-in expression and phenotype data
##

exprs = parse.tab(opt$expression)
if (is.null(opt$phenotypes)) {
	print("Error: must supply phenotype file!")
	q();
} else {
	pheno = parse.phenotypes(opt$phenotypes)
}

# check phenotype data matches
#
if (length(setdiff(rownames(pData(pheno)),colnames(exprs)))!=0) {
	print ("Error: phenotype and expression sample lists must be identical!")
	q();
}

if (is.null(opt$num_combin)) {
	# default to top 25 regulators for synnergy analysis
	opt$num_combin <- 25
}

expset.obj = makeExpressionSet(exprs, pheno)
print (expset.obj)
# parse the adj file to get a regulon object:
# note: all candidate regulators and genes must be in the 
# dataset you parsed
regulon <- opt$regulon 
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


# display the top X master regulators
num_results <- as.numeric(opt$num_results)
if (is.null(opt$num_results)) {
num_results <- 25
}
regulon_minsize <- as.numeric(opt$regulon_minsize)
if (is.null(opt$regulon_minsize)) {
regulon_minsize <- 25
}
num_permutations <- as.numeric(opt$permutations)
if (is.null(opt$permutations)) {
num_permutations <- 1000
}

mrs = run.marina(expset.obj, regul, opt$test_phenotype, opt$reference_phenotype, regul.minsize=regulon_minsize, num.permutations=num_permutations)
mr.summary = summary(mrs,mrs=length(mrs$regulon))

write.table(mr.summary, file=paste(opt$output, "/", "masterRegulators.txt", sep=""), col.names = NA, sep="\t", quote=F)

pdf(file=paste(opt$output, "/", "masterRegulators.pdf", sep=""))
plot(mrs,num_results,cex=0.7)
dev.off()

save.image(file=paste(opt$output, "/", "master-reg.RData", sep=""))

#synergy between MRs
#only compute if combinatorial analysis produces any synergy regulons
mrs.combin <- msviperCombinatorial(mrs,regulators = opt$num_combin,level=5,minsize=opt$regulon_minsize)
mrs.combin <- ledge(mrs.combin)

if(length(mrs.combin$regulon)>length(mrs$regulon)){
	mrs.synergy <- msviperSynergy(mrs.combin)
	synergy.summary <- summary(mrs.synergy,mrs=length(mrs.synergy$regulon))
	write.table(synergy.summary, file=paste(opt$output, "/", "masterRegulatorsSynergy.txt", sep=""), row.names=FALSE, sep="\t", quote=F)

	pdf(file=paste(opt$output, "/", "masterRegulatorsSynergy.pdf", sep=""))
	plot(mrs.synergy,mrs=min(opt$num_results,length(mrs.synergy$regulon)),cex=0.7)
	dev.off()
}

#leading edge
mr.le <- lapply(names(mrs$regulon),function(x) mrs$ledge[[x]])
names(mr.le) <- names(mrs$regulon)
writeSetList(mr.le,out.file=paste(opt$output, "/", "masterRegulatorsLedge.listt", sep=""))

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

#shadow analysis
mrs.shadow.001 <- shadow(mrs,regulators=0.1,shadow=0.01)
mrs.shadow.001.shadow.tfs <- lapply(unique(mrs.shadow.001$shadow[,2]),function(x) mrs.shadow.001$shadow[mrs.shadow.001$shadow[,2]==x,1])
names(mrs.shadow.001.shadow.tfs) <- unique(mrs.shadow.001$shadow[,2])

if(length(mrs.shadow.001.shadow.tfs)>0){
mrs.shadow.001.shadow.tfs <- mrs.shadow.001.shadow.tfs[names(sort(sapply(mrs.shadow.001.shadow.tfs,length),decreasing = T))]
writeSetList(mrs.shadow.001.shadow.tfs,out.file=paste(opt$output, "/", "shadow_single_tfs_pval_01.listt", sep=""))
}


save.image(file=paste(opt$output, "/", "master-reg.RData", sep=""))
q();


