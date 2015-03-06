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
hugo.ids <- Entrez_to_HugoGene(rownames(data))

## map <- cbind(entrez.ids, rownames(data))
## write.table(map, file="entrez2hugo.txt", sep='\t', quote=FALSE)

rownames(data) <- hugo.ids

write.df(data, out.file="data.hugo.tab")
