#!/usr/bin/env Rscript

library('getopt')

opt = getopt(matrix(c(
	'input', 'i', 1, "character",
	'output', 'o', 1, "character"
	),ncol=4,byrow=TRUE));

# Input:
#	- A reference gene set tab-delimited file where the first column is a list of gene identifiers and the 
#		second column is the database for which that gene id is annotated from. The first row must be labels of
#		these columns of data as this script skips the first row of data.
#
#	- The name of the output file that the script writes the converted gene ids to 'o'.
#	
#	This script detects gene ids from 'Ensembl', 'Entrez Gene', and 'Uniprot-TrEMBL' databases and converts
#		the gene ids to gene symbols.
#	
#	NOTE: this script cannot convert PubChem-compound and Chemspider ids to gene symbol ids
#	as there exists no framework currently available to support those annotation conversions.
#
# Returns: 
#	- A reference gene set file formatted the same way as the input file except the gene ids are converted to
#		their respective gene symbol id.

library(org.Hs.eg.db)
library(biomaRt)

parse.tab <- function (file) {

  m <- read.table(file, sep="\t", header=TRUE, quote="", check.names=FALSE)
  return (as.matrix(m))
}

mart=useMart("ensembl", dataset="hsapiens_gene_ensembl")
mart2 = useMart("unimart", dataset="uniprot")

geneList = parse.tab(opt$input)

printStr = paste('Identifier', 'Database', sep='	')
cat(printStr, file=opt$output, sep="\n")

for (geneIdx in 1:length(geneList[,1])) {

	geneID = geneList[geneIdx,1]
	geneID_type = geneList[geneIdx,2]

	#print(paste(geneID, geneID_type, sep=' '))

	g = ''

	if (geneID_type == "Ensembl"){
		g = getGene( id = geneID, type = "ensembl_gene_id", mart = mart)[2]
	} else if (geneID_type == "Entrez Gene"){
		g = getGene( id = geneID, type = "entrezgene", mart = mart)[2]
	} else if (geneID_type == "Uniprot-TrEMBL"){
		g = getGene( id = geneID, type = "uniprot_sptrembl", mart = mart)[2]
		if (nrow(g) == 0){
			g = getGene( id = geneID, type = "uniprot_swissprot", mart = mart)[2]
			if (nrow(g) == 0){
				g = getBM(attributes=c("gene_name"),filter="accession",values=geneID ,mart=mart2)
			}
		}
	}
	g = unique(g)
	if ((nrow(g) > 0) && !(g == '') && !is.na(g)){
		printStr = paste(g, 'gene symbol', sep='\t')
		cat(printStr, file=opt$output, sep="\n", append=TRUE)
	}
}
