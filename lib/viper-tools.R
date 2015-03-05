

# viper/mixtools import 
library(mixtools)
library(viper)
library(Biobase)

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

