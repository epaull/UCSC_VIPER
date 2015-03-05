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
    'num_combin','c',2,"double",
    'min_size', 'm',2, "double",
    'viper_null', 'x', 0, "logical"
    ),ncol=4,byrow=TRUE));

#options(error = quote({dump.frames(to.file = TRUE); q(status = 1)}))
library(mixtools)
#library(bcellViper)
library(viper)
library(Biobase)
#source("./shadow_combin.R")
# Input: 
#   - A TAB separated file: the first line is column header (sample) info. 
#       The first column has feature names 
#
# Returns:
#
#       - A data matrix with the 'colnames' and 'rownames' attributes set, the first column
#       a 'character' vector with the description/HUGO id's, and the following columns are
#       numeric data. 
parse.tab <- function (file) {

  m <- read.table(file, sep="\t", row.names=1, header=TRUE, quote="", check.names=FALSE)
  return (as.matrix(m))
}

parse.phenotypes <- function (file) {
	m <- read.table(file, sep="\t", row.names=1, header=TRUE, quote="", check.names=FALSE)
	phenoData <- new("AnnotatedDataFrame", data = m)
	return (phenoData)
}

#
# Input:  
#       - A data matrix with columns as samples, rows are genes
#       - phenotype data should include sample annotations (tumor/normal classes)
#
# Returns:
#
#       A data frame of predictors, filtered down the the top X
#
makeExpressionSet <- function(data.matrix, phenoData) {

  exprSet <- new("ExpressionSet", exprs = data.matrix, phenoData = phenoData, 
		 annotation = "viper_input")

	return (exprSet)
}


writeSetList <- function(setList,out.file,delim="\t"){
        #since we are appending to a file,
        #we need to check if it exists and remove it
        if(out.file!=stdout() && file.exists(out.file)){
                    file.remove(out.file)
                }
            for (i in 1:length(setList)){
                        write(paste(names(setList)[i],paste(setList[[i]],collapse=delim),sep=delim),file=out.file,append=TRUE)
                    }
    }



run.viper.MR <- function (exp.obj, regulon, set1.label, set2.label,min.size=25,perm.count=300) {

	# get set 1 indexes
	set1.idx <- which(colnames(exp.obj) %in% rownames(pData(exp.obj))[which(pData(exp.obj)[,1] == set1.label)])
	set2.idx <- which(colnames(exp.obj) %in% rownames(pData(exp.obj))[which(pData(exp.obj)[,1] == set2.label)])

	print (paste("Comparing ", length(set1.idx), " test samples with ", length(set2.idx), " reference samples..." ))

	# get just the data matrix from this object
	data.matrix = exprs(exp.obj)
	# generate the signature and normalize to z-scores
	signature <- rowTtest( data.matrix[,set1.idx], data.matrix[,set2.idx] )
    
    signature <- (qnorm(signature$p.value/2, lower.tail=F) * sign(signature$statistic))[, 1]
	# create the null model signatures
    
	nullmodel <- ttestNull(data.matrix[,set1.idx], data.matrix[,set2.idx], per=perm.count, repos=T)
    
	# compute MARINa scores based on the null model
	mrs <- msviper(signature, regulon, nullmodel,minsize=min.size)
    mrs <- ledge(mrs)
	print (summary(mrs))

	# create a background viper signature based on relative levels, then compute the final scores
	vpres <- NULL
	if (!is.null(opt$viper_null)) {
		print ("Constructing Viper Signature")
		vpsig <- viperSignature(data.matrix[,-set2.idx], data.matrix[,set2.idx], method="ttest", verbose=T)
		print ("Constructing Viper Inferences")
		vpres <- viper(vpsig, regulon,minsize=min.size)
	} else {
		print ("Constructing Viper Inferences")
		vpres <- viper(exp.obj, regulon,minsize=min.size)
	}

	return (list(mrs, vpres))
}

##
## Read-in expression and phenotype data
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
expset.obj = makeExpressionSet(exprs, pheno)
print (expset.obj)
# parse the adj file to get a regulon object:
# note: all candidate regulators and genes must be in the 
# dataset you parsed
regulon <- opt$regulon 
regul <- aracne2regulon(regulon, expset.obj)

save.image(file=paste(opt$output, "/", "master-reg.RData", sep=""))
result = run.viper.MR(expset.obj, regul, opt$test_phenotype, opt$reference_phenotype,min.size=opt$min_size)
mr.result = summary(result[[1]],mrs=length(result[[1]]$regulon))


viper.result = t(exprs(result[[2]]))

write.table(mr.result, file=paste(opt$output, "/", "masterRegulators.txt", sep=""),row.names=FALSE,sep="\t", quote=F)
write.table(cbind(id=rownames(viper.result),viper.result), file=paste(opt$output, "/", "viperScores.txt", sep=""), row.names=FALSE,sep="\t", quote=F)

pdf(file=paste(opt$output, "/", "masterRegulators.pdf", sep=""))
plot(result[[1]],mrs=min(opt$num_results,length(result[[1]]$regulon)),cex=0.7)
dev.off()

#leading edge
mr.le <- lapply(names(result[[1]]$regulon),function(x) result[[1]]$ledge[[x]])
names(mr.le) <- names(result[[1]]$regulon)

writeSetList(mr.le,out.file=paste(opt$output, "/", "masterRegulatorsLedge.listt", sep=""))

#combinatorial analysis
mrs.combin <- msviperCombinatorial(result[[1]],regulators = opt$num_combin,level=5,minsize=opt$min_size)
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
if(length(mrs.combin$regulon)>length(result[[1]]$regulon)){
mrs.synergy <- msviperSynergy(mrs.combin)
synergy.summary <- summary(mrs.synergy,mrs=length(mrs.synergy$regulon))
write.table(synergy.summary, file=paste(opt$output, "/", "masterRegulatorsSynergy.txt", sep=""), row.names=FALSE, sep="\t", quote=F)

pdf(file=paste(opt$output, "/", "masterRegulatorsSynergy.pdf", sep=""))
plot(mrs.synergy,mrs=min(opt$num_results,length(mrs.synergy$regulon)),cex=0.7)
dev.off()
}

#shadow analysis
mrs.shadow.001 <- shadow(result[[1]],regulators=0.1,shadow=0.01)
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


