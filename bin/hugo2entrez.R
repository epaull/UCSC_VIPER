#!/usr/bin/env	Rscript

# convert the names
library(GenomicFeatures)
library(mygene)

library('getopt')

# the real network score --
opt = getopt(matrix(c(
    'data', 'd', 1, "character"), ncol=4, byrow=TRUE));

## source the library, always in this relative location
calling_directory = dirname(get_Rscript_filename())
source(paste(calling_directory, "../lib", "viper-tools.R", sep="/"))

data <- parse.tab(opt$data)
entrez.ids <- HugoGene_to_Entrez(rownames(data))

map <- cbind(entrez.ids, rownames(data))
write.table(map, file="entrez2hugo.txt", sep='\t', quote=FALSE)

rownames(data) <- entrez.ids

write.table(data, file="data.entrez.tab", sep="\t", quote=FALSE)
