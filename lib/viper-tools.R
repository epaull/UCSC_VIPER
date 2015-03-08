
# viper/mixtools import 
library(mixtools)
library(viper)
library(Biobase)
library(GenomicFeatures)
library(mygene)


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


HugoGene_to_Entrez <- function (GeneList) {
  GenEnt <- queryMany(GeneList, scopes="symbol", fields="entrezgenes", species="human", size=1)
  GenEnt <- GenEnt[,"_id"]
}

Entrez_to_HugoGene <- function (GeneList) {
  GenSym <- queryMany(GeneList, fields="symbol", species="human", size=1)
  GenSym <- GenSym[,"symbol"]
}

write.df <- function(df,row.names.id='',out.file){
  output <- cbind(rownames(df),df)
  colnames(output)[1] <-row.names.id
  write.table(output ,file=out.file,quote=FALSE,append=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
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


run.marina <- function (exp.obj, regulon, set1.label, set2.label, regul.minsize, num.permutations) {

	##
	## Run MARINa based on the expression set object and regulon object, using the provided dichotomy. 
	##
	## Inputs:
	##	exp.obj - expressionSet object
	##	regulon - regulon object type
	##	set1.label - text label for the first phenotype, stored in exp.obj 
	##	set2.label - text for the second phenotype
	##	regul.minsize - minimum number of downstream genes to consider a regulator in this test 
	##	num.permutations - how many permutations to do for the background model
	##
	## Returns:
	##	a msviper MARINa object, output by 'msviper' function

	# get set 1 indexes
	set1.idx <- which(colnames(exp.obj) %in% rownames(pData(exp.obj))[which(pData(exp.obj)[,1] == set1.label)])
	set2.idx <- which(colnames(exp.obj) %in% rownames(pData(exp.obj))[which(pData(exp.obj)[,1] == set2.label)])

	print (paste("Comparing ", length(set1.idx), " test samples with ", length(set2.idx), " reference samples..." ))
	print (paste("Running ", as.numeric(num.permutations), " permutations for null model..."))

	# get just the data matrix from this object
	data.matrix = exprs(exp.obj)
	# generate the signature and normalize to z-scores
	signature <- rowTtest( data.matrix[,set1.idx], data.matrix[,set2.idx] )
	signature <- (qnorm(signature$p.value/2, lower.tail=F) * sign(signature$statistic))[, 1]
	# create the null model signatures
	nullmodel <- ttestNull(data.matrix[,set1.idx], data.matrix[,set2.idx], per=as.numeric(num.permutations), repos=T)
	# compute MARINa scores based on the null model
	mrs <- msviper(signature, regulon, nullmodel, minsize=regul.minsize)
    mrs <- ledge(mrs)

	return (mrs)
}

run.viper.unsupervised <- function (exp.obj, regulon, regul.minsize) {
	vpres <- viper(exp.obj, regulon, minsize=regul.minsize)
	return (vpres)
}

run.viper.supervised <- function (exp.obj, regulon, set1.label, set2.label,	min.size=25, perm.count=1000) {

	##
	## Run VIPER supervised 
	##

	# get set 1 indexes
	set1.idx <- which(colnames(exp.obj) %in% rownames(pData(exp.obj))[which(pData(exp.obj)[,1] == set1.label)])
	set2.idx <- which(colnames(exp.obj) %in% rownames(pData(exp.obj))[which(pData(exp.obj)[,1] == set2.label)])

	print (paste("Comparing ", length(set1.idx), " test samples with ", length(set2.idx), " reference samples..." ))
	# get just the data matrix from this object
	data.matrix = exprs(exp.obj)

	# create a background viper signature based on relative levels, then compute the final scores
	vpres <- NULL
	print ("Constructing Viper Signature")
	vpsig <- viperSignature(data.matrix[,-set2.idx], data.matrix[,set2.idx], method="ttest", verbose=T)
	print ("Constructing Viper Inferences")
	vpres <- viper(vpsig, regulon,minsize=min.size)

	return (vpres)
}

